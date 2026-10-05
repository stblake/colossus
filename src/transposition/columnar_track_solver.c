#include "columnar_track_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Columnar with within-column track permutation L (TRANSCOL_L)
// =====================================================================
//
//  Unlike the plain columnar solver (which recovers only the column order), this
//  solver also recovers a uniform within-column row permutation L -- equivalently a
//  permutation of the R grid rows, the jarl / "dave transposition" scheme.
//
//  Search strategy (from ciphers/W168/gridlib.py): the engine anneals the COLUMN
//  ORDER under the cheap IDENTITY-L reading (read the grid rows top-to-bottom). Even
//  when the true L is not identity this still discriminates the column order, because
//  the within-row n-grams -- which dominate the score -- are already correct at the
//  true column order; only the R-1 row seams are wrong. Nesting the exact best-L
//  inside every score eval is both far slower AND worse for the search (it gives every
//  wrong column order its best-case row reading, flattening the contrast). So best-L
//  is applied ONCE, at report time, to recover L exactly via the Held-Karp seam
//  decomposition (seam_best_row_order, trans_common.c) on the winning column order.
//
//  The exact best-L needs a complete grid (len % K == 0); a ragged grid degrades
//  gracefully to a plain columnar (identity L). The column count K is swept like the
//  single-columnar solver, optionally crossed with the column read direction
//  (-readdir) and the row read direction (-readrowdir).

// Word-coverage weights (Rec 2). The space-preserving jarl recipe: reward column
// orders that produce word-complete rows during the search, and break ties in the
// best-L row ordering with dictionary coverage. CT_WORD_SEARCH is normalized by
// length (so it is comparable to the mean n-gram score); CT_WORD_SEAM is a raw
// per-letter reward in the additive seam objective.
#define CT_WORD_SEARCH 0.6
#define CT_WORD_SEAM   2.0
static const WordSet *ct_ws = NULL;   // built once per solve (NULL => no dictionary)

// Per-decrypt scratch (single-threaded program; same idiom as the engine).
static int    ct_grid[MAX_CIPHER_LENGTH];
static int    ct_rowbuf[MAX_CIPHER_LENGTH];
static int   *ct_rowptr[MAX_CIPHER_LENGTH];
static int    ct_rowlen[MAX_CIPHER_LENGTH];
static double ct_indiv[HELD_KARP_MAX_NODES];
static double ct_delta[HELD_KARP_MAX_NODES * HELD_KARP_MAX_NODES];

// Best-L decode (REPORT only): refill the grid from the K columns (read per
// col_dir), build the row strings (reversed if row_rev), recover the best L exactly
// (Held-Karp seam), and emit the best-L plaintext into out[]. L_out (length R)
// receives the recovered row order. Returns R. For a ragged grid L is identity and
// out is the plain columnar decrypt.
static int coltrack_decode(const int *cipher, int len, int K, int col_dir, int row_rev,
                           const int *order, const float *ngram_data, int ngram_size,
                           int *L_out, int *out) {
    if (K <= 1 || len % K != 0) {
        decrypt_columnar((int *) cipher, len, K, (int *) order, col_dir, out);
        int R = (len + K - 1) / K;
        for (int i = 0; i < R; i++) L_out[i] = i;
        return R;
    }
    int R = len / K;
    int pos = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j];
        if (col_dir == COL_READ_BT) for (int r = R - 1; r >= 0; r--) ct_grid[r * K + c] = cipher[pos++];
        else                        for (int r = 0; r < R; r++)     ct_grid[r * K + c] = cipher[pos++];
    }
    for (int r = 0; r < R; r++) {
        for (int c = 0; c < K; c++)
            ct_rowbuf[r * K + c] = ct_grid[r * K + (row_rev ? (K - 1 - c) : c)];
        ct_rowptr[r] = &ct_rowbuf[r * K];
        ct_rowlen[r] = K;
    }
    seam_best_row_order(R, ct_rowptr, ct_rowlen, ngram_data, ngram_size,
                        ct_ws, CT_WORD_SEAM, ct_indiv, ct_delta, L_out);
    int o = 0;
    for (int i = 0; i < R; i++) {
        int r = L_out[i];
        for (int c = 0; c < K; c++) out[o++] = ct_rowbuf[r * K + c];
    }
    return R;
}

// =====================================================================
//  Structural crib-anchored columnar (-cribanchored)  (Rec 3)
// =====================================================================
//
//  For a complete-grid keyed columnar (L = identity) the cipher's R-char blocks ARE
//  the grid columns in a permuted order. A crib fixes some plaintext cells, hence
//  some characters of each grid column; a cipher block can fill grid column gc only
//  if it agrees with the crib at every known row of that column. We build that
//  block<->column compatibility, then backtrack over a perfect matching
//  (most-constrained column first), scoring each complete assignment by the n-gram +
//  word-coverage objective and keeping the best. A crib of even a few rows collapses
//  the (otherwise intractable) shallow K-column search to a unique (or tiny) set --
//  the technique that cracked jarl. Assumes L = identity (the standard keyed
//  columnar); the row read direction is not applied here.
#define CBA_LEAF_CAP 3000000L

static int cba_K, cba_R, cba_len, cba_dir;
static const int *cba_cipher;
static int cba_known[MAX_CIPHER_LENGTH];           // known[r*K+gc] = crib char index, or -1
static int cba_compat[MAX_COLS][MAX_COLS], cba_ncompat[MAX_COLS];
static int cba_colorder[MAX_COLS], cba_assign[MAX_COLS];
static bool cba_used[MAX_COLS];
static const float *cba_ng; static int cba_ns; static const WordSet *cba_wset;
static double cba_best; static int cba_best_order[MAX_COLS]; static long cba_leaves;

static void cba_eval_leaf(void) {
    int order[MAX_COLS];
    for (int gc = 0; gc < cba_K; gc++) order[cba_assign[gc]] = gc;  // block -> grid column
    static int out[MAX_CIPHER_LENGTH];
    decrypt_columnar((int *) cba_cipher, cba_len, cba_K, order, cba_dir, out);
    double s = ngram_sum_raw(out, cba_len, cba_ng, cba_ns);
    if (cba_wset) s += CT_WORD_SEAM * word_coverage(out, cba_len, cba_wset);
    if (s > cba_best) { cba_best = s; for (int j = 0; j < cba_K; j++) cba_best_order[j] = order[j]; }
}

static void cba_recurse(int idx) {
    if (cba_leaves > CBA_LEAF_CAP) return;
    if (idx == cba_K) { cba_eval_leaf(); cba_leaves++; return; }
    int gc = cba_colorder[idx];
    for (int t = 0; t < cba_ncompat[gc]; t++) {
        int b = cba_compat[gc][t];
        if (cba_used[b]) continue;
        cba_used[b] = true; cba_assign[gc] = b;
        cba_recurse(idx + 1);
        cba_used[b] = false;
    }
}

// Run the matcher for one (K, dir). Returns the best objective and fills order[0..K-1]
// (decrypt_columnar convention: cipher block j -> grid column order[j]); -1e300 if no
// crib-consistent perfect matching exists.
static double cba_solve_one(const int *cipher, int len, int K, int dir,
                            const int *crib_indices, const int *crib_positions, int n_cribs,
                            const float *ng, int ns, const WordSet *wset, int *order_out) {
    if (K < 2 || len % K != 0) return -1e300;
    int R = len / K;
    cba_K = K; cba_R = R; cba_len = len; cba_dir = dir; cba_cipher = cipher;
    cba_ng = ng; cba_ns = ns; cba_wset = wset;

    for (int i = 0; i < R * K; i++) cba_known[i] = -1;
    for (int i = 0; i < n_cribs; i++)
        if (crib_positions[i] >= 0 && crib_positions[i] < R * K)
            cba_known[crib_positions[i]] = crib_indices[i];

    // Build block<->column compatibility. Grid cell (r,gc) takes block b's element
    // b[r] (dir TB) or b[R-1-r] (dir BT); compatible if it matches every known cell.
    for (int gc = 0; gc < K; gc++) {
        cba_ncompat[gc] = 0;
        for (int b = 0; b < K; b++) {
            int ok = 1;
            for (int r = 0; r < R && ok; r++) {
                int k = cba_known[r * K + gc];
                if (k < 0) continue;
                int br = (dir == COL_READ_BT) ? (R - 1 - r) : r;
                if (cipher[b * R + br] != k) ok = 0;
            }
            if (ok) cba_compat[gc][cba_ncompat[gc]++] = b;
        }
        if (cba_ncompat[gc] == 0) return -1e300;      // a column no block can fill
    }
    // Most-constrained-first column order (simple selection sort on ncompat).
    for (int i = 0; i < K; i++) cba_colorder[i] = i;
    for (int i = 0; i < K; i++) {
        int m = i;
        for (int j = i + 1; j < K; j++)
            if (cba_ncompat[cba_colorder[j]] < cba_ncompat[cba_colorder[m]]) m = j;
        int t = cba_colorder[i]; cba_colorder[i] = cba_colorder[m]; cba_colorder[m] = t;
    }
    for (int b = 0; b < K; b++) cba_used[b] = false;
    cba_best = -1e300; cba_leaves = 0;
    cba_recurse(0);
    if (cba_best <= -1e300) return -1e300;
    for (int j = 0; j < K; j++) order_out[j] = cba_best_order[j];
    return cba_best;
}

// Entry: sweep the column count (and read direction) for the structural crib match,
// report the best. Used when -cribanchored is set with cribs present.
static void solve_columnar_track_crib(ColossusConfig *cfg, SharedData *shared,
    int cipher[], int cipher_len, int crib_indices[], int crib_positions[], int n_cribs,
    char *cribtext_str) {

    const WordSet *wset = trans_word_set(shared);
    int lo = cfg->min_cols, hi = cfg->max_cols;
    int cap = cipher_len / 2;
    if (lo < 2) lo = 2;
    if (hi > cap) hi = cap;
    if (hi > MAX_COLS) hi = MAX_COLS;

    int dirs[2], nd = 0;
    if (cfg->read_direction == COL_READ_BOTH) { dirs[nd++] = COL_READ_TB; dirs[nd++] = COL_READ_BT; }
    else dirs[nd++] = cfg->read_direction;

    double best = -1e300; int best_K = 0, best_dir = COL_READ_TB, best_order[MAX_COLS];
    for (int K = lo; K <= hi; K++) {
        if (cipher_len % K != 0) continue;
        for (int di = 0; di < nd; di++) {
            int order[MAX_COLS];
            double s = cba_solve_one(cipher, cipher_len, K, dirs[di],
                crib_indices, crib_positions, n_cribs, shared->ngram_data, cfg->ngram_size, wset, order);
            if (s > best) { best = s; best_K = K; best_dir = dirs[di];
                            for (int j = 0; j < K; j++) best_order[j] = order[j]; }
        }
    }

    if (best <= -1e300) {
        printf("\ntranscol-L (crib-anchored): no crib-consistent columnar found "
               "(check -mincols/-maxcols, the crib, and that L = identity).\n\n");
        return;
    }

    int out[MAX_CIPHER_LENGTH];
    decrypt_columnar(cipher, cipher_len, best_K, best_order, best_dir, out);
    double score = state_score(out, cipher_len, crib_indices, crib_positions, n_cribs,
        shared->ngram_data, cfg->ngram_size, cfg->weight_ngram, cfg->weight_crib,
        cfg->weight_ioc, cfg->weight_entropy);
    char params[64];
    snprintf(params, sizeof params, "K=%d dir=%s crib-anchored", best_K,
        best_dir == COL_READ_BT ? "bt" : "tb");
    printf("\ntranscol-L (crib-anchored): %d columns (%s), L = identity\n",
        best_K, best_dir == COL_READ_BT ? "bottom-to-top" : "top-to-bottom");
    printf("column order:");
    for (int j = 0; j < best_K; j++) printf(" %d", best_order[j]);
    printf("\n");
    report_transposition(cfg, shared, cipher, cipher_len, out, score,
        cribtext_str, n_cribs, params);
}

// ---- engine model ----------------------------------------------------------
// State: key[0..K-1] = column order. Config: period = K, aux[0] = col read dir,
// aux[1] = row read dir (reversal flag).

static void ct_krange(const SolverCtx *ctx, int *lo, int *hi) {
    int l = ctx->cfg->min_cols, h = ctx->cfg->max_cols;
    int cap = ctx->cipher_len / 2;
    if (cap < 2) cap = 2;
    if (cap > MAX_COLS) cap = MAX_COLS;
    if (l < 2) l = 2;
    if (h > cap) h = cap;
    if (l > h) l = h;
    *lo = l; *hi = h;
}

static int ct_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    int lo, hi;
    ct_krange(ctx, &lo, &hi);
    int cdirs[2], ncd = 0, rdirs[2], nrd = 0;
    if (ctx->cfg->read_direction == COL_READ_BOTH) { cdirs[ncd++] = COL_READ_TB; cdirs[ncd++] = COL_READ_BT; }
    else cdirs[ncd++] = ctx->cfg->read_direction;
    if (ctx->cfg->read_row_direction == ROW_READ_BOTH) { rdirs[nrd++] = ROW_READ_LR; rdirs[nrd++] = ROW_READ_RL; }
    else rdirs[nrd++] = ctx->cfg->read_row_direction;

    int n = 0;
    for (int K = lo; K <= hi; K++)
        for (int ci = 0; ci < ncd; ci++)
            for (int ri = 0; ri < nrd && n < cap; ri++) {
                out[n].period = K; out[n].aux[0] = cdirs[ci]; out[n].aux[1] = rdirs[ri];
                out[n].j = 0; out[n].k = 0;
                n++;
            }
    return n;
}

static void ct_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    perm_seed(st->key, cc->period);
}
static void ct_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st, bool *fp) {
    (void) ctx; (void) fp;
    perm_move(st->key, cc->period);
}
static void ct_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}
// SEARCH decrypt (per eval): cheap IDENTITY-L reading -- read the grid rows
// top-to-bottom (within-row reversed if row_rev). Discriminates the column order
// without the cost / contrast-flattening of nesting best-L. Ragged grids fall back
// to plain columnar.
static int ct_id_L[MAX_CIPHER_LENGTH];
static void ct_decrypt(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                       int *out, double *adj) {
    int len = ctx->cipher_len, K = cc->period;
    if (K <= 1 || len % K != 0) {
        decrypt_columnar((int *) ctx->cipher, len, K, st->key, cc->aux[0], out);
    } else {
        int R = len / K;
        for (int i = 0; i < R; i++) ct_id_L[i] = i;
        decrypt_columnar_tracked((int *) ctx->cipher, len, K, st->key, cc->aux[0],
                                 ct_id_L, cc->aux[1], out);
    }
    // Reward column orders that yield word-complete rows (the jarl recipe). Robust
    // to row order, since words sit within rows and identity-L preserves each row.
    if (ct_ws && len > 0) {
        double w = (ctx->cfg->weight_word > 0.0) ? ctx->cfg->weight_word : CT_WORD_SEARCH;
        *adj += w * word_coverage(out, len, ct_ws) / (double) len;
    }
}

static void ct_param_summary(const SolverConfig *cc, char *buf, size_t n) {
    snprintf(buf, n, "K=%d dir=%s row=%s", cc->period,
        cc->aux[0] == COL_READ_BT ? "bt" : "tb",
        cc->aux[1] == ROW_READ_RL ? "rl" : "lr");
}

static void ct_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) st;
    char params[64]; ct_param_summary(cc, params, sizeof params);
    report_transposition_verbose(ctx, score, decrypted, stats, params);
}

static void ct_report(const SolverCtx *ctx, const SolverConfig *cc, const SolverState *st,
                      double score, int *decrypted) {
    (void) score; (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    char params[64]; ct_param_summary(cc, params, sizeof params);
    // Apply the exact best-L to the winning column order to recover L + final plaintext.
    int L[MAX_CIPHER_LENGTH], out[MAX_CIPHER_LENGTH];
    int R = coltrack_decode(ctx->cipher, ctx->cipher_len, cc->period, cc->aux[0], cc->aux[1],
                            st->key, ctx->ngram_data, cfg->ngram_size, L, out);
    double best_l_score = state_score(out, ctx->cipher_len,
        ctx->crib_indices, ctx->crib_positions, ctx->n_cribs, ctx->ngram_data, cfg->ngram_size,
        cfg->weight_ngram, cfg->weight_crib, cfg->weight_ioc, cfg->weight_entropy);
    printf("\ntranscol-L: %d columns (%s), within-column track L over %d rows (row %s)\n",
        cc->period, cc->aux[0] == COL_READ_BT ? "bottom-to-top" : "top-to-bottom", R,
        cc->aux[1] == ROW_READ_RL ? "right-to-left" : "left-to-right");
    printf("column order:");
    for (int c = 0; c < cc->period; c++) printf(" %d", st->key[c]);
    printf("\ntrack order L:");
    for (int i = 0; i < R; i++) printf(" %d", L[i]);
    printf("\n");
    report_transposition(cfg, ctx->shared, ctx->cipher, ctx->cipher_len, out,
        best_l_score, ctx->cribtext, ctx->n_cribs, params);
}

static const CipherModel COLUMNAR_TRACK_MODEL = {
    .name = "columnar-track", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = ct_enumerate, .key_len = NULL,
    .seed = ct_seed, .perturb = ct_perturb, .copy_state = ct_copy,
    .decrypt = ct_decrypt, .report = ct_report, .report_verbose = ct_report_verbose,
};

void solve_columnar_track(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a transcol-L solve.\n\n");
        return;
    }
    ct_ws = trans_word_set(shared);     // dictionary coverage signal (NULL if no dict)

    // -cribanchored: use the crib as a STRUCTURAL constraint (block<->column matching)
    // rather than the stochastic anneal -- the only reliable attack on a shallow
    // many-column keyed columnar (the jarl collapse). Requires a crib.
    if (cfg->crib_anchored && n_cribs > 0) {
        solve_columnar_track_crib(cfg, shared, cipher_indices, cipher_len,
            crib_indices, crib_positions, n_cribs, cribtext_str);
        return;
    }

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    run_solver(&COLUMNAR_TRACK_MODEL, &ctx);
}

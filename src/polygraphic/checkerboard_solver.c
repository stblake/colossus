#include "checkerboard_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Checkerboard solver (TYPE checkerboard)
// =====================================================================
//
// THE CIPHER. A keyed 5x5 Polybius square (25 letters, J->I). Each plaintext letter -> a DIGRAPH
// (row label, column label). SIMPLE case: one label per row/column, so the digraph is unique --
// a bijection between 25 codes and 25 letters. COMPLEX case: two labels per row/column, so a
// letter has 2x2 codes chosen freely -> HOMOPHONIC (100 codes -> 25 letters).
//
// KEY REDUCTION. Label ORDER (which observed label denotes line 0 vs line 1) is NOT identifiable
// ciphertext-only: it is absorbed by a row/column permutation of the recovered square (the same
// fact nihilist_sub_solver.c documents for its keyed labels). So the ACA label keywords
// (BLACK/WHITE/...) are unrecoverable; only the square, and -- in the complex case -- the label
// PAIRING (which two labels share a line), are identifiable. Applied per axis:
//   SIMPLE  axis: the observed labels, sorted, ARE the lines (0..side-1). Nothing to search.
//   COMPLEX axis: the two labels of a line have identical conditional other-axis distributions,
//     so a chi-squared homogeneity distance + minimum-weight pairing over the observed labels
//     recovers the pairing INDEPENDENTLY of the other axis and of the square (a decoupling reward).
//
// SEARCH. Once the per-axis grouping (label letter -> line 0..side-1) is fixed, each digraph (a,b)
// collapses to a single MERGED CODE rowgroup[a]*side + colgroup[b] in 0..24, and the decrypt is
// out[i] = square[merged[i]] -- exactly the Aristocrat/homophonic kernel over 25 symbols. The
// square is a free 25-permutation annealed on the homophonic INCREMENTAL fast path (sync_caches /
// score_neighbor / commit_neighbor + per-thread scratch_clone), a swap of two codes' plaintext
// images per move. The pre-pass keeps the top-K pairings per axis and CROSSES them into engine
// configs (one merged stream each); the simple case is a single config.
//
// CAPABILITY. The simple case is a 25-letter Aristocrat: reliable from ~150-200 letters, partial
// and high-variance in the ACA 60-90 range. The complex case sits BELOW the ACA range: the pairing
// statistic has an O(N) bias favouring wrong pairings against an O(N^2) signal, so it needs
// ~400-600+ plaintext letters (800-1200 ciphertext) before the true pairing ranks first. This is a
// property of the cipher (documented, like Tridigital / CM-Bifid even periods), not a solver bug.

// ---------------------------------------------------------------------
//  Tunables and file-static problem buffers (written once before the search, then read-only)
// ---------------------------------------------------------------------

#define CB_KEEP_DEFAULT   4     // top-K pairings kept per complex axis (raise with -nperiods)
#define CB_MAX_KEEP       6
#define CB_MAX_CONFIGS    (CB_MAX_KEEP * CB_MAX_KEEP)
#define CB_MAX_PARTITIONS 2048  // >= 945 (the 10-into-5-pairs count); safety cap on the enumerator

static int g_cb_cipher[MAX_CIPHER_LENGTH];                 // letters-only digraph stream (0..24)
static int g_cb_merged[CB_MAX_CONFIGS * MAX_CIPHER_LENGTH]; // per-config merged 25-code streams

typedef struct {
    int    N;                  // plaintext length == scoring length == clen/2
    int    side, alpha;        // 5, 25
    int    row_complex, col_complex;
    int    nR, nC;             // distinct observed row / column labels
    int    n_configs;
    int    rlabels[CHECKERBOARD_MAX_LABELS];   // observed row-label letters, sorted
    int    clabels[CHECKERBOARD_MAX_LABELS];   // observed col-label letters, sorted
    int    rowgrp[CB_MAX_CONFIGS][CHECKERBOARD_MAX_LABELS]; // per config: rlabels position -> line
    int    colgrp[CB_MAX_CONFIGS][CHECKERBOARD_MAX_LABELS];

    int    eng_rank[CHECKERBOARD_MAX_GRID];    // English letters by descending monogram weight

    // --- incremental-scoring caches (homophonic fast path; simple case) ---
    double scale;              // n-gram scale factor (matches ngram_score())
    double ngsum;              // running raw sum of ngram_data[] over all windows
    int   *pos;                // length N: cipher positions grouped by merged code
    int   *pos_off;            // length alpha + 1
    char  *win_mark;           // length N (touched-window set scratch; stays all-zero)
    int   *win_list;           // length N
    double pend_ngsum;
    int    pend_nsym;
    int   *pend_sym;           // length alpha: changed merged codes
    int   *pend_newc;          // length alpha: their new plaintext letters
} CheckerboardScratch;

// ---------------------------------------------------------------------
//  Pre-pass: per-axis chi-squared homogeneity pairing (complex case)
// ---------------------------------------------------------------------

typedef struct { double dist; int grp[CHECKERBOARD_MAX_LABELS]; } CbPartition;
static CbPartition g_cb_parts[CB_MAX_PARTITIONS];

// Symmetric 2xC chi-squared homogeneity statistic between two count profiles a[], b[] (length C).
// Small => the two labels have the same other-axis distribution => they belong on the same line.
static double cb_profile_distance(const int a[], const int b[], int C) {
    long na = 0, nb = 0;
    for (int c = 0; c < C; c++) { na += a[c]; nb += b[c]; }
    if (na == 0 || nb == 0) return 0.0;    // an unseen label pairs freely (no evidence)
    double chi2 = 0.0;
    for (int c = 0; c < C; c++) {
        double tot = (double) a[c] + b[c];
        if (tot <= 0.0) continue;
        // expected a[c] under homogeneity = na*tot/(na+nb); (obs-exp)^2/exp summed over both rows.
        double e_a = (double) na * tot / (double) (na + nb);
        double e_b = (double) nb * tot / (double) (na + nb);
        chi2 += (a[c] - e_a) * (a[c] - e_a) / e_a + (b[c] - e_b) * (b[c] - e_b) / e_b;
    }
    return chi2;
}

// Recursive enumerator: partition the m observed labels into groups of size 1 or 2, at most
// max_groups groups, accumulating the within-pair distance. Appends each complete partition (with
// group ids 0..) to out[]; caps at `cap`.
typedef struct {
    int m, max_groups, cap, count;
    const double (*dist)[CHECKERBOARD_MAX_LABELS];
    CbPartition *out;
    int grp[CHECKERBOARD_MAX_LABELS];
    char used[CHECKERBOARD_MAX_LABELS];
    double acc;
} CbEnum;

static void cb_enum_rec(CbEnum *e, int ngroups) {
    int u = -1;
    for (int i = 0; i < e->m; i++) if (!e->used[i]) { u = i; break; }
    if (u < 0) {                                   // all labels grouped -> record
        if (e->count < e->cap) {
            e->out[e->count].dist = e->acc;
            for (int i = 0; i < e->m; i++) e->out[e->count].grp[i] = e->grp[i];
        }
        e->count++;
        return;
    }
    if (ngroups >= e->max_groups) return;          // no room for another group
    // u as a singleton line.
    e->used[u] = 1; e->grp[u] = ngroups;
    cb_enum_rec(e, ngroups + 1);
    e->used[u] = 0;
    // u paired with a later label v.
    for (int v = u + 1; v < e->m; v++) {
        if (e->used[v]) continue;
        e->used[u] = e->used[v] = 1; e->grp[u] = e->grp[v] = ngroups;
        e->acc += e->dist[u][v];
        cb_enum_rec(e, ngroups + 1);
        e->acc -= e->dist[u][v];
        e->used[u] = e->used[v] = 0;
    }
}

// Rank the pairings of `m` observed labels (profiles prof[m][C]) and write the best `keep`
// groupings (label position -> line 0..side-1) into grp_out[keep][m]. Returns the count kept.
static int cb_rank_pairings(int m, int side, const int prof[][CHECKERBOARD_MAX_LABELS], int C,
                            int keep, int grp_out[][CHECKERBOARD_MAX_LABELS], bool verbose) {
    double dist[CHECKERBOARD_MAX_LABELS][CHECKERBOARD_MAX_LABELS];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            dist[i][j] = (i == j) ? 0.0 : cb_profile_distance(prof[i], prof[j], C);

    CbEnum e;
    e.m = m; e.max_groups = side; e.cap = CB_MAX_PARTITIONS; e.count = 0;
    e.dist = dist; e.out = g_cb_parts; e.acc = 0.0;
    for (int i = 0; i < m; i++) { e.used[i] = 0; e.grp[i] = 0; }
    cb_enum_rec(&e, 0);

    int total = (e.count < e.cap) ? e.count : e.cap;
    if (total < 1) return 0;
    // Partial selection sort of the top `keep` by ascending distance.
    if (keep > total) keep = total;
    for (int k = 0; k < keep; k++) {
        int best = k;
        for (int j = k + 1; j < total; j++)
            if (g_cb_parts[j].dist < g_cb_parts[best].dist) best = j;
        CbPartition tmp = g_cb_parts[k]; g_cb_parts[k] = g_cb_parts[best]; g_cb_parts[best] = tmp;
        for (int i = 0; i < m; i++) grp_out[k][i] = g_cb_parts[k].grp[i];
    }
    if (verbose)
        printf("  pairing pre-pass: %d partitions, best chi2=%.2f, kept %d\n",
            total, g_cb_parts[0].dist, keep);
    return keep;
}

// ---------------------------------------------------------------------
//  CipherModel hooks (shared by the simple and complex models)
// ---------------------------------------------------------------------

// n-gram window packing, identical to ngram_score() (see aristocrat_solver.c).
#define CB_WINDOW_INDEX(w, ng, letter) ({                   \
    int _idx = 0;                                           \
    for (int _j = 0; _j < (ng); _j++)                       \
        _idx = _idx * g_alpha + letter((w) + _j);           \
    _idx; })

// Swap two merged codes' plaintext images: keeps key[] a permutation of 0..alpha-1.
static void cb_move(int *key, int alpha) {
    int a = rand_int(0, alpha), b = rand_int(0, alpha);
    int t = key[a]; key[a] = key[b]; key[b] = t;
}

static int cb_enumerate_simple(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = g_alpha; out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

static int cb_enumerate_complex(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const CheckerboardScratch *s = (const CheckerboardScratch *) ctx->model_scratch;
    int n = s->n_configs; if (n > cap) n = cap;
    for (int c = 0; c < n; c++) {
        out[c].period = g_alpha; out[c].j = 0; out[c].k = 0;
        out[c].aux[0] = c;        // config id -> which merged stream to decrypt
        out[c].aux[1] = 0;
    }
    return n;
}

// Frequency rank-match warm start over the config's merged codes + a few random swaps.
static void cb_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const CheckerboardScratch *s = (const CheckerboardScratch *) ctx->model_scratch;
    const int *merged = &g_cb_merged[(size_t) cc->aux[0] * s->N];
    int hist[CHECKERBOARD_MAX_GRID], rank[CHECKERBOARD_MAX_GRID];
    for (int c = 0; c < g_alpha; c++) { hist[c] = 0; rank[c] = c; }
    for (int i = 0; i < s->N; i++) hist[merged[i]]++;
    for (int x = 0; x < g_alpha; x++)
        for (int y = x + 1; y < g_alpha; y++)
            if (hist[rank[y]] > hist[rank[x]]) { int t = rank[x]; rank[x] = rank[y]; rank[y] = t; }
    for (int r = 0; r < g_alpha; r++) st->key[rank[r]] = s->eng_rank[r];
    st->key_len = g_alpha;
    int kicks = rand_int(0, 9);
    for (int k = 0; k < kicks; k++) cb_move(st->key, g_alpha);
}

static void cb_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    cb_move(st->key, g_alpha);
    if (frand() < 0.10) cb_move(st->key, g_alpha);
}

static void cb_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < g_alpha; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void cb_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, int *out, double *score_adjust) {
    const CheckerboardScratch *s = (const CheckerboardScratch *) ctx->model_scratch;
    const int *merged = &g_cb_merged[(size_t) cc->aux[0] * s->N];
    for (int i = 0; i < s->N; i++) out[i] = st->key[merged[i]];
    *score_adjust = 0.0;
}

// ---------------------------------------------------------------------
//  Incremental fast path (simple case; homophonic/aristocrat twin over the merged codes)
// ---------------------------------------------------------------------

static void cb_sync_caches(const SolverCtx *ctx, const SolverConfig *cc, const int *dec) {
    (void) cc;
    CheckerboardScratch *s = (CheckerboardScratch *) ctx->model_scratch;
    int len = ctx->cipher_len, ng = ctx->cfg->ngram_size;
    s->ngsum = 0.0;
    int n_windows = len - ng + 1;
    #define dec_at(q) (dec[q])
    for (int w = 0; w < n_windows; w++)
        s->ngsum += ctx->ngram_data[CB_WINDOW_INDEX(w, ng, dec_at)];
    #undef dec_at
}

static double cb_score_neighbor(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *cur, const SolverState *loc, const int *cur_dec, double cur_score) {
    (void) cur_score;
    CheckerboardScratch *s = (CheckerboardScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int len = ctx->cipher_len, ng = cfg->ngram_size, n = cc->period;
    const int *cipher = ctx->cipher;                 // the merged-code stream (config 0)

    s->pend_nsym = 0;
    for (int sym = 0; sym < n; sym++) {
        if (loc->key[sym] == cur->key[sym]) continue;
        s->pend_sym[s->pend_nsym]  = sym;
        s->pend_newc[s->pend_nsym] = loc->key[sym];
        s->pend_nsym++;
    }

    double dsum = 0.0;
    int n_windows = len - ng + 1;
    if (cfg->weight_ngram > 1.e-4 && n_windows > 0) {
        int nlist = 0;
        for (int k = 0; k < s->pend_nsym; k++) {
            int sym = s->pend_sym[k];
            for (int pi = s->pos_off[sym]; pi < s->pos_off[sym + 1]; pi++) {
                int p = s->pos[pi];
                int lo = p - ng + 1; if (lo < 0) lo = 0;
                int hi = p;          if (hi > n_windows - 1) hi = n_windows - 1;
                for (int w = lo; w <= hi; w++)
                    if (!s->win_mark[w]) { s->win_mark[w] = 1; s->win_list[nlist++] = w; }
            }
        }
        #define dec_at(q) (cur_dec[q])
        #define loc_at(q) (loc->key[cipher[q]])
        for (int i = 0; i < nlist; i++) {
            int w = s->win_list[i];
            int old_idx = CB_WINDOW_INDEX(w, ng, dec_at);
            int new_idx = CB_WINDOW_INDEX(w, ng, loc_at);
            dsum += ctx->ngram_data[new_idx] - ctx->ngram_data[old_idx];
            s->win_mark[w] = 0;
        }
        #undef loc_at
        #undef dec_at
    }
    s->pend_ngsum = s->ngsum + dsum;

    if (cfg->weight_ngram > 1.e-4 && n_windows > 0)
        return s->scale * s->pend_ngsum / (len - ng);
    return 0.0;
}

static void cb_commit_neighbor(const SolverCtx *ctx, const SolverConfig *cc, int *cur_dec) {
    (void) cc;
    CheckerboardScratch *s = (CheckerboardScratch *) ctx->model_scratch;
    for (int k = 0; k < s->pend_nsym; k++) {
        int sym = s->pend_sym[k], c = s->pend_newc[k];
        for (int pi = s->pos_off[sym]; pi < s->pos_off[sym + 1]; pi++)
            cur_dec[s->pos[pi]] = c;
    }
    s->ngsum = s->pend_ngsum;
}

static void *cb_scratch_clone(const SolverCtx *ctx) {
    const CheckerboardScratch *src = (const CheckerboardScratch *) ctx->model_scratch;
    int len = src->N;
    CheckerboardScratch *s = malloc(sizeof *s);
    if (!s) return NULL;
    memcpy(s, src, sizeof *s);                       // share read-only fields (pos/pos_off, ranks)
    s->win_mark  = calloc(len > 0 ? len : 1, sizeof(char));
    s->win_list  = malloc(sizeof(int) * (len > 0 ? len : 1));
    s->pend_sym  = malloc(sizeof(int) * CHECKERBOARD_GRID);
    s->pend_newc = malloc(sizeof(int) * CHECKERBOARD_GRID);
    if (!s->win_mark || !s->win_list || !s->pend_sym || !s->pend_newc) {
        free(s->win_mark); free(s->win_list); free(s->pend_sym); free(s->pend_newc);
        free(s);
        return NULL;
    }
    return s;
}

static void cb_scratch_free(void *scratch) {
    CheckerboardScratch *s = (CheckerboardScratch *) scratch;
    if (!s) return;
    free(s->win_mark); free(s->win_list); free(s->pend_sym); free(s->pend_newc);
    free(s);
}

// ---------------------------------------------------------------------
//  Reporting
// ---------------------------------------------------------------------

static void cb_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) st;
    report_transposition_verbose(ctx, score, decrypted, stats, "checkerboard");
}

// Print the recovered square: cell (line rg, line cg) holds key[rg*side+cg].
static void cb_print_square(const int key[], int side) {
    for (int rg = 0; rg < side; rg++) {
        printf("    ");
        for (int cg = 0; cg < side; cg++) printf("%c ", index_to_char(key[rg * side + cg]));
        printf("\n");
    }
}

static void cb_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const CheckerboardScratch *s = (const CheckerboardScratch *) ctx->model_scratch;
    int N = s->N, side = s->side, cid = cc->aux[0];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < N; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[N] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    const char *rc = s->row_complex ? "complex" : "simple";
    const char *cc_case = s->col_complex ? "complex" : "simple";
    printf("\nResult Score: %.2f | Words: %d | rows=%s(%d) cols=%s(%d)\n",
        score, n_words_found, rc, s->nR, cc_case, s->nC);

    print_text(decrypted, N);
    printf("\n%s\n", ctx->cribtext);

    printf("\nrecovered 5x5 square (plaintext letters by line; label ORDER is not identifiable):\n");
    cb_print_square(st->key, side);

    // Show the observed labels grouped onto lines for the winning config.
    printf("  row lines:");
    for (int line = 0; line < side; line++) {
        printf(" [");
        int first = 1;
        for (int p = 0; p < s->nR; p++)
            if (s->rowgrp[cid][p] == line) { printf("%s%c", first ? "" : "", index_to_char(s->rlabels[p])); first = 0; }
        printf("]");
    }
    printf("\n  col lines:");
    for (int line = 0; line < side; line++) {
        printf(" [");
        int first = 1;
        for (int p = 0; p < s->nC; p++)
            if (s->colgrp[cid][p] == line) { printf("%s%c", first ? "" : "", index_to_char(s->clabels[p])); first = 0; }
        printf("]");
    }
    printf("\n");

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(decrypted, ctx->result->decrypted, N);
        ctx->result->decrypted_len = N;
    }

    // One-liner: >>> score, [words,] type, rows=, cols=, file, PLAINTEXT (letters-only, so the
    // accuracy suite's whitespace-stripped compare matches the .solution).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, rows=%s, cols=%s, ", score, n_words_found, cfg->cipher_type, rc, cc_case);
    else
        printf(">>> %.2f, %d, rows=%s, cols=%s, ", score, cfg->cipher_type, rc, cc_case);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(decrypted, N);
    printf("\n");
}

static const CipherModel CHECKERBOARD_SIMPLE_MODEL = {
    .name = "checkerboard", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = cb_enumerate_simple, .key_len = NULL,
    .seed = cb_seed, .perturb = cb_perturb, .copy_state = cb_copy,
    .decrypt = cb_decrypt_hook, .report = cb_report, .report_verbose = cb_report_verbose,
    .sync_caches = cb_sync_caches, .score_neighbor = cb_score_neighbor,
    .commit_neighbor = cb_commit_neighbor,
    .scratch_clone = cb_scratch_clone, .scratch_free = cb_scratch_free,
};

static const CipherModel CHECKERBOARD_COMPLEX_MODEL = {
    .name = "checkerboard", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = cb_enumerate_complex, .key_len = NULL,
    .seed = cb_seed, .perturb = cb_perturb, .copy_state = cb_copy,
    .decrypt = cb_decrypt_hook, .report = cb_report, .report_verbose = cb_report_verbose,
    // No incremental hooks: K configs each carry a distinct merged stream, so the generic
    // per-move re-decrypt is used (correct under multi-config + -nthreads; the complex case is a
    // documented low-yield limitation where the incremental speedup does not matter).
};

// ---------------------------------------------------------------------
//  Entry point
// ---------------------------------------------------------------------

// Sorted distinct labels seen at digraph position `parity` (0 = row, 1 = col); returns the count.
static int cb_observed_labels(const int cipher[], int clen, int parity, int labels[]) {
    char seen[MAX_ALPHABET_SIZE] = {0};
    int n = 0, ndig = clen / 2;
    for (int i = 0; i < ndig; i++) {
        int v = cipher[2 * i + parity];
        if (v >= 0 && v < MAX_ALPHABET_SIZE && !seen[v]) { seen[v] = 1; labels[n++] = v; }
    }
    for (int x = 0; x < n; x++)                       // insertion sort (n <= 10)
        for (int y = x + 1; y < n; y++)
            if (labels[y] < labels[x]) { int t = labels[x]; labels[x] = labels[y]; labels[y] = t; }
    return n;
}

void solve_checkerboard(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;
    int side = CHECKERBOARD_SIDE, alpha = side * side;
    if (g_alpha != alpha) {
        printf("\n\nERROR: Checkerboard needs a %d-letter alphabet (got %d). "
               "Run -type checkerboard so the alphabet is forced (J->I).\n\n", alpha, g_alpha);
        return;
    }

    // Letters-only digraph stream (drop any spaces/punctuation sentinels).
    int L = 0;
    for (int i = 0; i < cipher_len && L < MAX_CIPHER_LENGTH; i++)
        if (cipher_indices[i] >= 0) g_cb_cipher[L++] = cipher_indices[i];
    int N = L / 2;                                    // plaintext length (drop a dangling letter)
    if (N < 4) {
        printf("\n\nERROR: Checkerboard ciphertext too short (%d digraphs); need >= 4.\n\n", N);
        return;
    }

    CheckerboardScratch scratch;
    memset(&scratch, 0, sizeof scratch);
    scratch.N = N; scratch.side = side; scratch.alpha = alpha;

    scratch.nR = cb_observed_labels(g_cb_cipher, L, 0, scratch.rlabels);
    scratch.nC = cb_observed_labels(g_cb_cipher, L, 1, scratch.clabels);
    scratch.row_complex = (scratch.nR > side);
    scratch.col_complex = (scratch.nC > side);

    // Inverse: label letter -> its position in the sorted observed-label list (else -1).
    int rpos[MAX_ALPHABET_SIZE], cpos[MAX_ALPHABET_SIZE];
    for (int i = 0; i < MAX_ALPHABET_SIZE; i++) { rpos[i] = -1; cpos[i] = -1; }
    for (int p = 0; p < scratch.nR; p++) rpos[scratch.rlabels[p]] = p;
    for (int p = 0; p < scratch.nC; p++) cpos[scratch.clabels[p]] = p;

    // Count matrix M[row label][col label] over the observed labels (the pre-pass profiles).
    static int M[CHECKERBOARD_MAX_LABELS][CHECKERBOARD_MAX_LABELS];
    for (int r = 0; r < scratch.nR; r++)
        for (int c = 0; c < scratch.nC; c++) M[r][c] = 0;
    for (int i = 0; i < N; i++) {
        int r = rpos[g_cb_cipher[2 * i]], c = cpos[g_cb_cipher[2 * i + 1]];
        if (r >= 0 && c >= 0) M[r][c]++;
    }

    int keep = (cfg->n_periods > 0) ? cfg->n_periods : CB_KEEP_DEFAULT;
    if (keep > CB_MAX_KEEP) keep = CB_MAX_KEEP;
    if (keep < 1) keep = 1;

    // Row groupings (label position -> line). Simple: the sorted labels ARE the lines. Complex:
    // the top-`keep` pairings from the chi-squared pre-pass over row profiles M[r][*].
    int rowgrp[CB_MAX_KEEP][CHECKERBOARD_MAX_LABELS];
    int colgrp[CB_MAX_KEEP][CHECKERBOARD_MAX_LABELS];
    int keep_r, keep_c;
    if (scratch.row_complex) {
        keep_r = cb_rank_pairings(scratch.nR, side, M, scratch.nC, keep, rowgrp, cfg->verbose);
        if (keep_r < 1) { for (int p = 0; p < scratch.nR; p++) rowgrp[0][p] = p % side; keep_r = 1; }
    } else {
        for (int p = 0; p < scratch.nR; p++) rowgrp[0][p] = p;   // identity: label p -> line p
        keep_r = 1;
    }
    if (scratch.col_complex) {
        // Column profiles are the columns of M; transpose into a scratch matrix.
        static int Mt[CHECKERBOARD_MAX_LABELS][CHECKERBOARD_MAX_LABELS];
        for (int c = 0; c < scratch.nC; c++)
            for (int r = 0; r < scratch.nR; r++) Mt[c][r] = M[r][c];
        keep_c = cb_rank_pairings(scratch.nC, side, Mt, scratch.nR, keep, colgrp, cfg->verbose);
        if (keep_c < 1) { for (int p = 0; p < scratch.nC; p++) colgrp[0][p] = p % side; keep_c = 1; }
    } else {
        for (int p = 0; p < scratch.nC; p++) colgrp[0][p] = p;
        keep_c = 1;
    }

    // Cross the per-axis groupings into configs; build each config's merged 25-code stream.
    int n_configs = keep_r * keep_c;
    if (n_configs > CB_MAX_CONFIGS) n_configs = CB_MAX_CONFIGS;
    scratch.n_configs = n_configs;
    for (int cr = 0; cr < keep_r; cr++)
        for (int cc = 0; cc < keep_c; cc++) {
            int cid = cr * keep_c + cc;
            if (cid >= n_configs) break;
            for (int p = 0; p < scratch.nR; p++) scratch.rowgrp[cid][p] = rowgrp[cr][p];
            for (int p = 0; p < scratch.nC; p++) scratch.colgrp[cid][p] = colgrp[cc][p];
            int *merged = &g_cb_merged[(size_t) cid * N];
            for (int i = 0; i < N; i++) {
                int rp = rpos[g_cb_cipher[2 * i]], cp = cpos[g_cb_cipher[2 * i + 1]];
                int rg = (rp >= 0) ? rowgrp[cr][rp] : 0;
                int cg = (cp >= 0) ? colgrp[cc][cp] : 0;
                merged[i] = rg * side + cg;
            }
        }

    bool complex_case = (scratch.row_complex || scratch.col_complex);
    if (cfg->verbose)
        printf("\ncheckerboard: %d digraphs, rows=%s(%d) cols=%s(%d), %d config(s)%s\n",
            N, scratch.row_complex ? "complex" : "simple", scratch.nR,
            scratch.col_complex ? "complex" : "simple", scratch.nC, n_configs,
            complex_case ? " [complex: recovery needs long text]" : "");

    // English letters by descending monogram weight (warm-start target ranking).
    for (int i = 0; i < g_alpha; i++) scratch.eng_rank[i] = i;
    for (int x = 0; x < g_alpha; x++)
        for (int y = x + 1; y < g_alpha; y++)
            if (g_monograms[scratch.eng_rank[y]] > g_monograms[scratch.eng_rank[x]]) {
                int t = scratch.eng_rank[x]; scratch.eng_rank[x] = scratch.eng_rank[y]; scratch.eng_rank[y] = t;
            }
    scratch.scale = g_ngram_logprob ? 1.0 : pow((double) g_alpha, cfg->ngram_size);

    // Cribs align to per-character ciphertext positions (2x the plaintext), so they do not line up
    // with the merged digraph stream -- ignored, like ADFGVX / nihilist-sub.
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    // The "cipher" the engine sees IS config 0's merged stream (the incremental fast path reads it).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_cb_merged, N, crib_indices, crib_positions, 0);
    ctx.result = result;

    if (complex_case) {
        ctx.model_scratch = &scratch;
        run_solver(&CHECKERBOARD_COMPLEX_MODEL, &ctx);
        return;
    }

    // Simple case: build the position index over the (single) merged stream for the fast path.
    scratch.pos      = malloc(sizeof(int) * (N > 0 ? N : 1));
    scratch.pos_off  = malloc(sizeof(int) * (alpha + 1));
    scratch.win_mark = calloc(N > 0 ? N : 1, sizeof(char));
    scratch.win_list = malloc(sizeof(int) * (N > 0 ? N : 1));
    scratch.pend_sym  = malloc(sizeof(int) * alpha);
    scratch.pend_newc = malloc(sizeof(int) * alpha);
    if (!scratch.pos || !scratch.pos_off || !scratch.win_mark ||
        !scratch.win_list || !scratch.pend_sym || !scratch.pend_newc) {
        printf("\n\nERROR: out of memory in checkerboard solve.\n\n");
        free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
        free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
        return;
    }
    for (int sym = 0; sym <= alpha; sym++) scratch.pos_off[sym] = 0;
    for (int i = 0; i < N; i++) scratch.pos_off[g_cb_merged[i] + 1]++;
    for (int sym = 0; sym < alpha; sym++) scratch.pos_off[sym + 1] += scratch.pos_off[sym];
    {
        int cursor[CHECKERBOARD_GRID];
        for (int sym = 0; sym < alpha; sym++) cursor[sym] = scratch.pos_off[sym];
        for (int i = 0; i < N; i++) scratch.pos[cursor[g_cb_merged[i]]++] = i;
    }

    ctx.model_scratch = &scratch;
    run_solver(&CHECKERBOARD_SIMPLE_MODEL, &ctx);

    free(scratch.pos); free(scratch.pos_off); free(scratch.win_mark);
    free(scratch.win_list); free(scratch.pend_sym); free(scratch.pend_newc);
}

int checkerboard_pairing_rank(const int cipher[], int clen, int axis, const int true_line[]) {
    int side = CHECKERBOARD_SIDE;
    int labels[CHECKERBOARD_MAX_LABELS], olabels[CHECKERBOARD_MAX_LABELS];
    int m  = cb_observed_labels(cipher, clen, axis, labels);
    int om = cb_observed_labels(cipher, clen, 1 - axis, olabels);
    if (m < 2) return 0;

    int pos[MAX_ALPHABET_SIZE], opos[MAX_ALPHABET_SIZE];
    for (int i = 0; i < MAX_ALPHABET_SIZE; i++) { pos[i] = -1; opos[i] = -1; }
    for (int p = 0; p < m; p++)  pos[labels[p]]   = p;
    for (int p = 0; p < om; p++) opos[olabels[p]] = p;

    static int P[CHECKERBOARD_MAX_LABELS][CHECKERBOARD_MAX_LABELS];
    for (int r = 0; r < m; r++) for (int c = 0; c < om; c++) P[r][c] = 0;
    int ndig = clen / 2;
    for (int i = 0; i < ndig; i++) {
        int a = cipher[2 * i + axis], b = cipher[2 * i + (1 - axis)];
        int r = (a >= 0 && a < MAX_ALPHABET_SIZE) ? pos[a] : -1;
        int c = (b >= 0 && b < MAX_ALPHABET_SIZE) ? opos[b] : -1;
        if (r >= 0 && c >= 0) P[r][c]++;
    }

    double dist[CHECKERBOARD_MAX_LABELS][CHECKERBOARD_MAX_LABELS];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            dist[i][j] = (i == j) ? 0.0 : cb_profile_distance(P[i], P[j], om);

    // The true partition's total within-pair distance (labels sharing a true line).
    double true_dist = 0.0;
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
            if (true_line[labels[i]] == true_line[labels[j]]) true_dist += dist[i][j];

    CbEnum e;
    e.m = m; e.max_groups = side; e.cap = CB_MAX_PARTITIONS; e.count = 0;
    e.dist = dist; e.out = g_cb_parts; e.acc = 0.0;
    for (int i = 0; i < m; i++) { e.used[i] = 0; e.grp[i] = 0; }
    cb_enum_rec(&e, 0);

    int total = (e.count < e.cap) ? e.count : e.cap;
    int rank = 0;
    for (int k = 0; k < total; k++)
        if (g_cb_parts[k].dist < true_dist - 1.e-9) rank++;
    return rank;
}

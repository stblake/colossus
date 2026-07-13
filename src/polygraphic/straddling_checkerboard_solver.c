#include "straddling_checkerboard_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Straddling Checkerboard solver -- see straddling_checkerboard_solver.h
// =====================================================================

// Decoded-plaintext + tiled scratch (single-threaded), off the stack.
static _Thread_local int g_sc_decode[MAX_CIPHER_LENGTH];
static _Thread_local int g_sc_tiled[MAX_CIPHER_LENGTH];
// Parsed ciphertext digit stream, filled by solve_straddling_checkerboard.
static int g_sc_digits[MAX_CIPHER_LENGTH];

typedef struct {
    int *digits;                             // ciphertext digit stream (length clen)
    int  clen;                               // digit-stream length (== engine cipher_len)
    int  ngram_size;
    float *ngram;
    int  n_configs;                          // number of (top-K) indicator-pair configs
    int  cfg_r0[STRADDLING_N_CONFIGS];       // config indicator pairs (r0 < r1), pre-pass ranked
    int  cfg_r1[STRADDLING_N_CONFIGS];
    int  warm[STRADDLING_N_CONFIGS][SC_MAP_CELLS];   // mini-solve warm-start map per kept config
} StraddlingScratch;

// ===================================================================
//  Code indexing + map decode (the solver's board-free decrypt)
// ===================================================================
//
// Within a config {r0,r1} the 28 token codes are indexed: the 8 non-indicator single digits
// -> 0..7 (via single_slot[]); the (r0,d2) doubles -> 8+d2; the (r1,d2) doubles -> 18+d2.
// map[code] in {0..25 letter, SC_MAP_FIG, SC_MAP_NULL}. A token's figure-mode digit is its
// own column digit (the single digit, or the double's 2nd digit) -- intrinsic to the code --
// so keyed labels need not be modelled. This matches the primitive's board decrypt exactly
// (the equivalent map is map[code(token)] = board.sym_at_cell[cell(token)]).

static void sc_single_slot(int r0, int r1, int single_slot[10]) {
    int s = 0;
    for (int d = 0; d < 10; d++) single_slot[d] = (d == r0 || d == r1) ? -1 : s++;
}

static int sc_decode_map(const int *digits, int clen, int r0, int r1, const int single_slot[10],
                         const int map[SC_MAP_CELLS], int out[], int filler,
                         int *n_tokens, int *n_valid) {
    int o = 0, nt = 0, nv = 0, mode = 0, i = 0;
    while (i < clen) {
        int g = digits[i], code, coldigit;
        if (g < 0 || g > 9) { i++; continue; }
        if (g == r0 || g == r1) {
            if (i + 1 >= clen) { out[o++] = filler; nt++; i++; continue; }   // truncated
            int g2 = digits[i + 1];
            i += 2;
            if (g2 < 0 || g2 > 9) { out[o++] = filler; nt++; continue; }
            code = ((g == r0) ? 8 : 18) + g2;
            coldigit = g2;
        } else {
            i += 1;
            code = single_slot[g];                                           // 0..7
            coldigit = g;
        }
        int s = map[code];
        if (s == SC_MAP_NULL) { out[o++] = filler; nt++; continue; }         // unused cell
        if (s == SC_MAP_FIG)  { mode ^= 1; nt++; nv++; continue; }           // figure-shift toggle
        nt++; nv++;
        out[o++] = (mode == 0) ? s : (26 + coldigit);                        // letter, or figure digit
    }
    if (n_tokens) *n_tokens = nt;
    if (n_valid)  *n_valid  = nv;
    return o;
}

// Score a decoded map: mean n-gram of the decode tiled to `clen` (length-fair) + validity.
static double sc_score_map(const StraddlingScratch *a, int r0, int r1, const int single_slot[10],
                           const int map[SC_MAP_CELLS]) {
    int nt = 0, nv = 0;
    int m = sc_decode_map(a->digits, a->clen, r0, r1, single_slot, map,
                          g_sc_decode, STRADDLING_FILLER, &nt, &nv);
    if (m <= 0) return -1e9;
    for (int i = 0; i < a->clen; i++) g_sc_tiled[i] = g_sc_decode[i % m];
    double s = ngram_score(g_sc_tiled, a->clen, a->ngram, a->ngram_size);
    s += STRADDLING_VALID_WEIGHT * ((nt > 0) ? (double) nv / nt : 0.0);
    return s;
}

// ===================================================================
//  Config pre-pass: a short substitution hill-climb per indicator pair (the real ranking)
// ===================================================================

// One simulated-annealing swap climb of the code->cell map from a random start; returns the
// best score reached and its map. Metropolis accept-worse (geometric T0 -> Tmin) escapes the
// gibberish local optima a greedy climb sticks in, so a WRONG config cannot masquerade as the
// true one and the true config reliably ranks near-optimal.
static double sc_anneal_map(const StraddlingScratch *a, int r0, int r1, const int single_slot[10],
                            int iters, double t0, double tmin, int best_map_out[SC_MAP_CELLS]) {
    int map[SC_MAP_CELLS];
    for (int i = 0; i < SC_MAP_CELLS; i++) map[i] = i;   // 0..25 letters, FIG, NULL
    shuffle(map, SC_MAP_CELLS);
    double cur = sc_score_map(a, r0, r1, single_slot, map);
    double best = cur;
    for (int i = 0; i < SC_MAP_CELLS; i++) best_map_out[i] = map[i];

    double ratio = (iters > 1) ? pow(tmin / t0, 1.0 / (iters - 1)) : 1.0, temp = t0;
    for (int it = 0; it < iters; it++, temp *= ratio) {
        int i = rand_int(0, SC_MAP_CELLS), j = rand_int(0, SC_MAP_CELLS);
        if (i == j) continue;
        int t = map[i]; map[i] = map[j]; map[j] = t;
        double s = sc_score_map(a, r0, r1, single_slot, map);
        double d = s - cur;
        if (d >= 0 || frand() < exp(d / temp)) {
            cur = s;
            if (s > best) { best = s; for (int k = 0; k < SC_MAP_CELLS; k++) best_map_out[k] = map[k]; }
        } else {
            t = map[i]; map[i] = map[j]; map[j] = t;                          // revert
        }
    }
    return best;
}

// Rank-and-warm-start a config: a few independent SA restarts, keep the best map + score.
static double sc_mini_solve(const StraddlingScratch *a, int r0, int r1, int iters,
                            int best_map_out[SC_MAP_CELLS]) {
    int single_slot[10];
    sc_single_slot(r0, r1, single_slot);
    int nrst = 3, per = iters / nrst;
    double best = -1e18;
    int mp[SC_MAP_CELLS];
    for (int r = 0; r < nrst; r++) {
        double s = sc_anneal_map(a, r0, r1, single_slot, per, 0.30, 0.02, mp);
        if (s > best) { best = s; for (int k = 0; k < SC_MAP_CELLS; k++) best_map_out[k] = mp[k]; }
    }
    return best;
}

// Rank all 45 indicator pairs by their mini-solve score; keep the top `keep` (with warm maps).
static void sc_rank_configs(StraddlingScratch *a, int keep, bool verbose) {
    int r0s[STRADDLING_N_CONFIGS], r1s[STRADDLING_N_CONFIGS];
    double sc[STRADDLING_N_CONFIGS];
    static int wm[STRADDLING_N_CONFIGS][SC_MAP_CELLS];
    int n = 0;
    for (int r0 = 0; r0 < 10; r0++)
        for (int r1 = r0 + 1; r1 < 10; r1++) {
            r0s[n] = r0; r1s[n] = r1;
            sc[n] = sc_mini_solve(a, r0, r1, STRADDLING_PREPASS_ITERS, wm[n]);
            n++;
        }
    if (keep > n) keep = n;
    for (int i = 0; i < keep; i++) {
        int best = i;
        for (int j = i + 1; j < n; j++) if (sc[j] > sc[best]) best = j;
        double ts = sc[i]; sc[i] = sc[best]; sc[best] = ts;
        int t0 = r0s[i]; r0s[i] = r0s[best]; r0s[best] = t0;
        int t1 = r1s[i]; r1s[i] = r1s[best]; r1s[best] = t1;
        for (int k = 0; k < SC_MAP_CELLS; k++) { int t = wm[i][k]; wm[i][k] = wm[best][k]; wm[best][k] = t; }
    }
    a->n_configs = keep;
    for (int i = 0; i < keep; i++) {
        a->cfg_r0[i] = r0s[i]; a->cfg_r1[i] = r1s[i];
        for (int k = 0; k < SC_MAP_CELLS; k++) a->warm[i][k] = wm[i][k];
    }
    if (verbose) {
        printf("\nstraddling pre-pass (mini-solve of all 45 indicator pairs, keep top %d):\n", keep);
        for (int i = 0; i < keep && i < 8; i++)
            printf("  {%d,%d} score=%.4f\n", r0s[i], r1s[i], sc[i]);
    }
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int sc_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const StraddlingScratch *a = (const StraddlingScratch *) ctx->model_scratch;
    int n = a->n_configs;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = i;                          // config index (for the warm-start lookup)
        out[i].j = 0; out[i].k = 0;
        out[i].aux[0] = a->cfg_r0[i];
        out[i].aux[1] = a->cfg_r1[i];
    }
    return n;
}

// Seed from the config's mini-solve warm map, jittered by a few random swaps for diversity.
static void sc_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const StraddlingScratch *a = (const StraddlingScratch *) ctx->model_scratch;
    int idx = cc->period;
    for (int i = 0; i < SC_MAP_CELLS; i++) st->key[i] = a->warm[idx][i];
    int jit = rand_int(0, 7);
    for (int s = 0; s < jit; s++) {
        int i = rand_int(0, SC_MAP_CELLS), j = rand_int(0, SC_MAP_CELLS);
        int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    }
    st->key_len = SC_MAP_CELLS;
}

static void sc_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int i = rand_int(0, SC_MAP_CELLS), j = rand_int(0, SC_MAP_CELLS);
    int t = st->key[i]; st->key[i] = st->key[j]; st->key[j] = t;
    if (frand() < 0.10) {                            // occasional second swap (coarser move)
        int p = rand_int(0, SC_MAP_CELLS), q = rand_int(0, SC_MAP_CELLS);
        t = st->key[p]; st->key[p] = st->key[q]; st->key[q] = t;
    }
}

static void sc_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < SC_MAP_CELLS; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void sc_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, int *out, double *score_adjust) {
    const StraddlingScratch *a = (const StraddlingScratch *) ctx->model_scratch;
    int single_slot[10];
    sc_single_slot(cc->aux[0], cc->aux[1], single_slot);
    int C = a->clen, nt = 0, nv = 0;
    int m = sc_decode_map(a->digits, C, cc->aux[0], cc->aux[1], single_slot, st->key,
                          g_sc_decode, STRADDLING_FILLER, &nt, &nv);
    if (m <= 0) {
        for (int i = 0; i < C; i++) out[i] = STRADDLING_FILLER;
        *score_adjust = 0.0;
        return;
    }
    for (int i = 0; i < C; i++) out[i] = g_sc_decode[i % m];
    *score_adjust = STRADDLING_VALID_WEIGHT * ((nt > 0) ? (double) nv / nt : 0.0);
}

// ===================================================================
//  Reporting
// ===================================================================

static char sc_cell_char(int s) {
    return (s == SC_MAP_FIG) ? '/' : (s == SC_MAP_NULL) ? '.' : ('A' + s);
}

// Render the recovered board (canonical 0..9 labels): top row (8 non-indicator singles),
// then the r0 and r1 rows. Keyed labels are not separately identifiable, so this is the
// equivalent fixed-label board.
static void sc_render_board(int r0, int r1, const int map[SC_MAP_CELLS], char out[]) {
    int single_slot[10];
    sc_single_slot(r0, r1, single_slot);
    int o = 0;
    for (int d = 0; d < 10; d++) if (single_slot[d] >= 0) out[o++] = sc_cell_char(map[single_slot[d]]);
    for (int d = 0; d < 10; d++) out[o++] = sc_cell_char(map[8 + d]);
    for (int d = 0; d < 10; d++) out[o++] = sc_cell_char(map[18 + d]);
    out[o] = '\0';
}

static void sc_print_digits(const int digits[], int n) {
    for (int i = 0; i < n; i++) putchar('0' + digits[i]);
}

static void sc_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    const StraddlingScratch *a = (const StraddlingScratch *) ctx->model_scratch;
    int r0 = cc->aux[0], r1 = cc->aux[1], C = a->clen, nt = 0, nv = 0;
    int single_slot[10];
    sc_single_slot(r0, r1, single_slot);
    int m = sc_decode_map(a->digits, C, r0, r1, single_slot, st->key,
                          g_sc_decode, STRADDLING_FILLER, &nt, &nv);
    if (m > MAX_CIPHER_LENGTH) m = MAX_CIPHER_LENGTH;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < m; i++) plaintext_string[i] = index_to_char(g_sc_decode[i]);
    plaintext_string[m] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char boardstr[SC_MAP_CELLS + 1]; sc_render_board(r0, r1, st->key, boardstr);

    printf("\nResult Score: %.2f | Words: %d | valid=%d/%d | indicators=%d%d | board=%s\n",
        score, n_words_found, nv, nt, r0, r1, boardstr);

    sc_print_digits(a->digits, C);
    printf("\n");
    print_text(g_sc_decode, m);
    printf("\n");
    print_spaces_line(g_spaces_table, g_sc_decode, m);
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(g_sc_decode, ctx->result->decrypted, m);
        ctx->result->decrypted_len = m;
    }

    // >>> score, [words,] type, board=, ind=, valid=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, board=%s, ind=%d%d, valid=%d/%d, ",
            score, n_words_found, cfg->cipher_type, boardstr, r0, r1, nv, nt);
    else
        printf(">>> %.2f, %d, board=%s, ind=%d%d, valid=%d/%d, ",
            score, cfg->cipher_type, boardstr, r0, r1, nv, nt);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    sc_print_digits(a->digits, C);
    printf(", ");
    print_text(g_sc_decode, m);
    printf("\n");
}

static const CipherModel STRADDLING_MODEL = {
    .name = "straddling-checkerboard", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = sc_enumerate, .key_len = NULL,
    .seed = sc_seed, .perturb = sc_perturb, .copy_state = sc_copy,
    .decrypt = sc_decrypt_hook, .report = sc_report,
};

// ===================================================================
//  Entry point
// ===================================================================

static int sc_parse_digits(const char *s, int out[], int cap) {
    int n = 0;
    for (int i = 0; s[i] && n < cap; i++)
        if (s[i] >= '0' && s[i] <= '9') out[n++] = s[i] - '0';
    return n;
}

void solve_straddling_checkerboard(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_len; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != MAX_ALPHABET_SIZE) {
        printf("\n\nERROR: Straddling Checkerboard needs the 36-symbol alphabet (A..Z + 0..9, "
               "got %d). Run -type sc so the alphabet is forced.\n\n", g_alpha);
        return;
    }

    int clen = sc_parse_digits(ciphertext_str, g_sc_digits, MAX_CIPHER_LENGTH);
    if (clen < 4) {
        printf("\n\nERROR: parsed only %d ciphertext digits; need >= 4. The Straddling "
               "Checkerboard ciphertext must be a stream of digits 0-9.\n\n", clen);
        return;
    }

    static StraddlingScratch scratch;
    scratch.digits = g_sc_digits;
    scratch.clen = clen;
    scratch.ngram_size = cfg->ngram_size;
    scratch.ngram = shared->ngram_data;

    int keep = (cfg->n_primers > 0) ? cfg->n_primers : STRADDLING_KEEP;
    sc_rank_configs(&scratch, keep, cfg->verbose);

    if (cfg->verbose)
        printf("\nstraddling-checkerboard: %d ciphertext digits, %d-symbol alphabet %s, %d config(s)\n",
            clen, g_alpha, g_idx_to_char_arr, scratch.n_configs);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, clen, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&STRADDLING_MODEL, &ctx);
}

#include "route_chain_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Route + column-key two-stage chain (TRANSROUTECOL)
// =====================================================================
//
//  Decrypt of one (route, shape, column-key) candidate:
//    1. SCATTER: place cipher[k] at the k-th cell of the route over the R x C grid
//       (the inverse of "read the grid along the route"). Complete grid: R*C == len.
//    2. COLUMN KEY: row r becomes  g[r][order[0]], g[r][order[1]], ...
//    3. SEAM BEST-L: read the R permuted rows in the exact best order (Held-Karp).
//  The engine anneals the column key; the route + grid shape + row direction are the
//  swept config. Reuses route_cells() (the 6 colossus routes) and the shared seam
//  best-L. See ciphers/W168/chains.py for the validated model.

#define RC_WORD_SEARCH 0.6
#define RC_WORD_SEAM   2.0
static const WordSet *rc_ws = NULL;

static int    rc_cells[MAX_CIPHER_LENGTH];
static int    rc_grid[MAX_CIPHER_LENGTH];
static int    rc_rowbuf[MAX_CIPHER_LENGTH];
static int   *rc_rowptr[MAX_CIPHER_LENGTH];
static int    rc_rowlen[MAX_CIPHER_LENGTH];
static double rc_indiv[HELD_KARP_MAX_NODES];
static double rc_delta[HELD_KARP_MAX_NODES * HELD_KARP_MAX_NODES];

// Decode one candidate into out[]; fills L_out (length R) with the row order and
// returns R. Requires a complete grid (R*C == len). use_best_l selects the exact
// Held-Karp row ordering (report) vs the cheap identity-L reading (search).
static int rc_decode(const int *cipher, int len, int R, int C, int route_id, int row_rev,
                     const int *order, const float *ngram_data, int ngram_size,
                     int use_best_l, int *L_out, int *out) {
    int n = route_cells(R, C, len, route_id, rc_cells);
    if (n != len) {
        for (int i = 0; i < len; i++) out[i] = cipher[i];
        for (int r = 0; r < R; r++) L_out[r] = r;
        return R;
    }
    for (int k = 0; k < len; k++) rc_grid[rc_cells[k]] = cipher[k];   // scatter (row-major grid)
    for (int r = 0; r < R; r++) {
        for (int c = 0; c < C; c++) {
            int src_c = order[row_rev ? (C - 1 - c) : c];
            rc_rowbuf[r * C + c] = rc_grid[r * C + src_c];
        }
        rc_rowptr[r] = &rc_rowbuf[r * C];
        rc_rowlen[r] = C;
    }
    if (use_best_l)
        seam_best_row_order(R, rc_rowptr, rc_rowlen, ngram_data, ngram_size,
                            rc_ws, RC_WORD_SEAM, rc_indiv, rc_delta, L_out);
    else
        for (int r = 0; r < R; r++) L_out[r] = r;       // identity-L reading (search)
    int o = 0;
    for (int i = 0; i < R; i++) { int r = L_out[i];
        for (int c = 0; c < C; c++) out[o++] = rc_rowbuf[r * C + c]; }
    return R;
}

// ---- engine model: state key[0..C-1] = column key; config period=C,
//      aux[0]=route_id, aux[1]=row_rev (R = len/C). ----
static int rc_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    int len = ctx->cipher_len;
    int rdirs[2], nrd = 0;
    if (ctx->cfg->read_row_direction == ROW_READ_BOTH) { rdirs[nrd++] = ROW_READ_LR; rdirs[nrd++] = ROW_READ_RL; }
    else rdirs[nrd++] = ctx->cfg->read_row_direction;
    int n = 0;
    for (int C = 2; C <= len / 2; C++) {
        if (len % C != 0) continue;                  // complete rectangles only
        int R = len / C;
        if (R < 2 || C > MAX_COLS) continue;
        for (int route_id = 0; route_id < N_ROUTES; route_id++)
            for (int ri = 0; ri < nrd && n < cap; ri++) {
                out[n].period = C; out[n].aux[0] = route_id; out[n].aux[1] = rdirs[ri];
                out[n].j = 0; out[n].k = 0; n++;
            }
    }
    return n;
}
static void rc_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx; perm_seed(st->key, cc->period);
}
static void rc_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st, bool *fp) {
    (void) ctx; (void) fp; perm_move(st->key, cc->period);
}
static void rc_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}
static void rc_decrypt(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                       int *out, double *adj) {
    int L[MAX_CIPHER_LENGTH];
    int len = ctx->cipher_len, C = cc->period, R = len / C;
    rc_decode(ctx->cipher, len, R, C, cc->aux[0], cc->aux[1],
              st->key, ctx->ngram_data, ctx->cfg->ngram_size, 0 /* identity-L search */, L, out);
    if (rc_ws && len > 0) {
        double w = (ctx->cfg->weight_word > 0.0) ? ctx->cfg->weight_word : RC_WORD_SEARCH;
        *adj += w * word_coverage(out, len, rc_ws) / (double) len;
    }
}
static const char *rc_route_names[N_ROUTES] = {
    "rows-snake", "cols-snake", "spiral-cw", "spiral-ccw", "diag-snake", "diag" };
static void rc_param_summary(const SolverCtx *ctx, const SolverConfig *cc, char *buf, size_t n) {
    int C = cc->period, R = ctx->cipher_len / C;
    snprintf(buf, n, "%dx%d route=%d row=%s", R, C, cc->aux[0],
        cc->aux[1] == ROW_READ_RL ? "rl" : "lr");
}
static void rc_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) st;
    char params[80]; rc_param_summary(ctx, cc, params, sizeof params);
    report_transposition_verbose(ctx, score, decrypted, stats, params);
}
static void rc_report(const SolverCtx *ctx, const SolverConfig *cc, const SolverState *st,
                      double score, int *decrypted) {
    (void) score; (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    char params[80]; rc_param_summary(ctx, cc, params, sizeof params);
    int L[MAX_CIPHER_LENGTH], out[MAX_CIPHER_LENGTH];
    int C = cc->period, R = ctx->cipher_len / C;
    rc_decode(ctx->cipher, ctx->cipher_len, R, C, cc->aux[0], cc->aux[1],
              st->key, ctx->ngram_data, cfg->ngram_size, 1 /* best-L */, L, out);
    double best_l_score = state_score(out, ctx->cipher_len,
        ctx->crib_indices, ctx->crib_positions, ctx->n_cribs, ctx->ngram_data, cfg->ngram_size,
        cfg->weight_ngram, cfg->weight_crib, cfg->weight_ioc, cfg->weight_entropy);
    printf("\ntransroutecol: %d x %d grid, route %d (%s), column key + seam best-L (row %s)\n",
        R, C, cc->aux[0], rc_route_names[cc->aux[0]], cc->aux[1] == ROW_READ_RL ? "rl" : "lr");
    printf("column key:");
    for (int c = 0; c < C; c++) printf(" %d", st->key[c]);
    printf("\ntrack order L:");
    for (int i = 0; i < R; i++) printf(" %d", L[i]);
    printf("\n");
    report_transposition(cfg, ctx->shared, ctx->cipher, ctx->cipher_len, out,
        best_l_score, ctx->cribtext, ctx->n_cribs, params);
}

static const CipherModel ROUTE_CHAIN_MODEL = {
    .name = "route-chain", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = rc_enumerate, .key_len = NULL,
    .seed = rc_seed, .perturb = rc_perturb, .copy_state = rc_copy,
    .decrypt = rc_decrypt, .report = rc_report, .report_verbose = rc_report_verbose,
};

void solve_route_chain(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str;
    if (cipher_len < 6) {
        printf("\n\nERROR: ciphertext too short for a transroutecol solve.\n\n");
        return;
    }
    rc_ws = trans_word_set(shared);
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    run_solver(&ROUTE_CHAIN_MODEL, &ctx);
}

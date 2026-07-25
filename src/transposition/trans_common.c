#include "trans_common.h"
#include "scoring.h"

// =====================================================================
//  Shared reporting for the dedicated transposition solvers
// =====================================================================
//
// Every transposition solver recovers a single best plaintext plus a short,
// type-specific parameter description; this helper prints the common
// human-readable block and the ">>> ..." one-line CSV summary so the output shape
// (and, crucially, the recovered plaintext as the final CSV field that the
// regression suite scrapes) is identical across types. `param_summary` is one
// already-formatted CSV field describing the recovered key (e.g. "rails=5 off=0").
void report_transposition(ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len, int best_decrypted[],
    double best_score, char *cribtext_str, int n_cribs,
    const char *param_summary) {

    int n_words_found = 0;

    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < cipher_len; i++) plaintext_string[i] = index_to_char(best_decrypted[i]);
    plaintext_string[cipher_len] = '\0';

    if (cfg->dictionary_present && shared->dict != NULL) {
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);
    }

    printf("\nResult Score: %.2f | Words: %d | %s\n", best_score, n_words_found, param_summary);

    print_text(cipher_indices, cipher_len);
    printf("\n");
    print_text(best_decrypted, cipher_len);
    printf("\n");
    print_solution_check(best_decrypted, cipher_len);
    print_spaces_line(g_spaces_table, best_decrypted, cipher_len);
    printf("%s\n", cribtext_str);

    if (PARTIAL_CRIB_MATCH && n_cribs > 0) {
        for (int i = 0; i < cipher_len; i++) {
            if (cribtext_str[i] == '_') {
                printf("_");
            } else {
                int diff = abs(best_decrypted[i] - (g_char_to_idx[toupper((unsigned char)cribtext_str[i]) & 127]));
                if (diff < 10) printf("%d", diff); else printf("*");
            }
        }
        printf("\n");
    }

    // One-liner summary: >>> score, [words,] type, <params>, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present) {
        printf(">>> %.2f, %d, %d, ", best_score, n_words_found, cfg->cipher_type);
    } else {
        printf(">>> %.2f, %d, ", best_score, cfg->cipher_type);
    }
    printf("%s, ", param_summary);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(cipher_indices, cipher_len);
    printf(", ");
    print_text(best_decrypted, cipher_len);
    printf("\n");
}

// Live (-verbose) progress block for the optimization-based transposition models, in
// the same shape as permutation_report_verbose / columnar_model_report_verbose: timing
// and search counters, a type-specific param line, and the current best plaintext.
// best_decrypted is the just-accepted best decrypt the engine passes to report_verbose.
void report_transposition_verbose(const SolverCtx *ctx, double best_score,
    int best_decrypted[], const EngineStats *stats, const char *param_summary) {

    double elapsed = engine_elapsed_sec(stats);
    double n_iter_per_sec = (elapsed > 0.) ? ((double) stats->n_iterations)/elapsed : 0.;
    printf("\n%.2f\t[sec]\n", elapsed);
    printf("%.0fK\t[it/sec/thread]\n", 1.e-3*n_iter_per_sec);
    printf("%d\t[restarts]\n", stats->n_restarts);
    printf("%d\t[backtracks]\n", stats->n_backtracks);
    printf("%d\t[slips]\n", stats->n_slips);
    printf("%.4f\t[entropy]\n", entropy(best_decrypted, ctx->cipher_len));
    printf("%.2f\t[score]\n", best_score);
    printf("%s\t[params]\n\n", param_summary);
    print_text(best_decrypted, ctx->cipher_len); printf("\n");
    print_solution_check(best_decrypted, ctx->cipher_len);
    fflush(stdout);
}


// =====================================================================
//  Shared key-anneal hooks for the permutation-style transposition models
// =====================================================================
//
// Amsco, Myszkowski, Redefence, Cadenus, Nihilist, Swagman and Grille all reduce
// to annealing a short integer key (st->key) through the generic engine. Only the
// key->plaintext decrypt, the neighbour move, the restart seed, the outer parameter
// sweep and the report differ per type. The move and seed are supplied through a
// TransKeyOps descriptor (placed in ctx->model_scratch by the solve_ entry point);
// the decrypt, enumerate and report are per-type hooks. cc->period carries the key
// length; cc->aux[0..1] carry the fixed per-config parameters (start / offset /
// readmode / N), with -variant read straight from cfg.


void tkey_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    ((const TransKeyOps *) ctx->model_scratch)->seed_cb(st->key, cc->period);
}
void tkey_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st, bool *fp) {
    (void) fp;
    ((const TransKeyOps *) ctx->model_scratch)->move_cb(st->key, cc->period);
}
void tkey_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < cc->period; i++) dst->key[i] = src->key[i];
}

// Neighbour move shared by the permutation-key types (swap dominant, with short
// reverses and block moves), preserving the permutation property.
void perm_move(int *key, int K) {
    if (K < 2) return;
    double r = frand();
    if (r < 0.70) {
        int a = rand_int(0, K), b = rand_int(0, K);
        int t = key[a]; key[a] = key[b]; key[b] = t;
    } else if (r < 0.85) {
        int max_blk = min(K, 8);
        int blk = rand_int(2, max_blk + 1);
        int s = rand_int(0, K - blk + 1);
        for (int a = s, b = s + blk - 1; a < b; a++, b--) {
            int t = key[a]; key[a] = key[b]; key[b] = t;
        }
    } else {
        int max_blk = min(K, 8);
        int blk = rand_int(1, max_blk + 1);
        int s = rand_int(0, K - blk + 1);
        int d = rand_int(0, K - blk + 1);
        if (d == s) return;
        int tmp[8];
        for (int a = 0; a < blk; a++) tmp[a] = key[s + a];
        if (d < s) { for (int a = s - 1; a >= d; a--) key[a + blk] = key[a]; }
        else       { for (int a = s + blk; a < d + blk; a++) key[a - blk] = key[a]; }
        for (int a = 0; a < blk; a++) key[d + a] = tmp[a];
    }
}

// Restart seed shared by the permutation-key types: a random permutation.
void perm_seed(int *key, int K) {
    for (int i = 0; i < K; i++) key[i] = i;
    shuffle(key, K);
}


int sweep_keylen(const SolverCtx *ctx, const SolverConfig *cc) { (void)ctx; (void)cc; return 0; }
void sweep_noop_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) { (void)ctx; (void)cc; (void)st; }
void sweep_noop_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) { (void)cc; (void)src; (void)dst; }
// Integer square root with exact-square test: returns N where N*N == x, else -1.
int exact_isqrt(int x) {
    if (x < 0) return -1;
    int n = (int)(sqrt((double)x) + 0.5);
    for (int d = -1; d <= 1; d++)
        if ((n + d) >= 0 && (n + d) * (n + d) == x) return n + d;
    return -1;
}


// =====================================================================
//  Held-Karp exact maximum-weight Hamiltonian path (see trans_common.h)
// =====================================================================
//
// Standard subset-DP: dp[mask][i] is the best score of a path visiting exactly the
// nodes in `mask` and ending at node i; the answer is max_i dp[full][i], with the
// order recovered by parent pointers. Single-threaded program => the DP tables are
// file-static scratch (not stack), sized to the 2^R cap. O(R^2 * 2^R).
#define HK_NEG (-1e300)
static double hk_dp[1 << HELD_KARP_MAX_NODES][HELD_KARP_MAX_NODES];
static int    hk_par[1 << HELD_KARP_MAX_NODES][HELD_KARP_MAX_NODES];

double held_karp_best_path(int R, const double *indiv, const double *delta, int *order_out) {
    if (R < 1) return 0.0;
    if (R == 1) { order_out[0] = 0; return indiv[0]; }
    if (R > HELD_KARP_MAX_NODES) {       // too large: identity order, identity-path score
        double s = indiv[0];
        order_out[0] = 0;
        for (int i = 1; i < R; i++) { order_out[i] = i; s += indiv[i] + delta[(i - 1) * R + i]; }
        return s;
    }

    int size = 1 << R;
    for (int mask = 0; mask < size; mask++)
        for (int i = 0; i < R; i++) { hk_dp[mask][i] = HK_NEG; hk_par[mask][i] = -1; }
    for (int i = 0; i < R; i++) hk_dp[1 << i][i] = indiv[i];

    for (int mask = 0; mask < size; mask++) {
        for (int i = 0; i < R; i++) {
            double base = hk_dp[mask][i];
            if (base == HK_NEG || !((mask >> i) & 1)) continue;
            const double *di = &delta[i * R];
            for (int j = 0; j < R; j++) {
                if ((mask >> j) & 1) continue;
                int nm = mask | (1 << j);
                double v = base + di[j] + indiv[j];
                if (v > hk_dp[nm][j]) { hk_dp[nm][j] = v; hk_par[nm][j] = i; }
            }
        }
    }

    int full = size - 1, bi = 0;
    double best = HK_NEG;
    for (int i = 0; i < R; i++)
        if (hk_dp[full][i] > best) { best = hk_dp[full][i]; bi = i; }

    // Walk parents back, then reverse into order_out.
    int tmp[HELD_KARP_MAX_NODES], n = 0, mask = full, i = bi;
    while (i != -1) { tmp[n++] = i; int pi = hk_par[mask][i]; mask ^= (1 << i); i = pi; }
    for (int k = 0; k < n; k++) order_out[k] = tmp[n - 1 - k];
    return best;
}

// Cached word-set keyed on the loaded dictionary pointer (single-threaded program).
static WordSet *g_trans_ws = NULL;
static char   **g_trans_ws_dict = NULL;
WordSet *trans_word_set(SharedData *shared) {
    if (!shared || shared->dict == NULL || shared->n_dict_words <= 0) return NULL;
    if (g_trans_ws && g_trans_ws_dict == shared->dict) return g_trans_ws;
    if (g_trans_ws) word_set_free(g_trans_ws);
    g_trans_ws = word_set_build(shared->dict, shared->n_dict_words);
    g_trans_ws_dict = shared->dict;
    return g_trans_ws;
}

// Per-row additive objective: raw within-word n-gram sum + optional dictionary
// word-coverage reward (both additive across a join, so the seam stays exact).
static double row_objective(const int *t, int len, const float *ngram_data, int ngram_size,
                            const WordSet *ws, double wword) {
    double s = ngram_sum_raw(t, len, ngram_data, ngram_size);
    if (ws && wword != 0.0) s += wword * word_coverage(t, len, ws);
    return s;
}

double seam_best_row_order(int R, int *const rows[], const int rowlen[],
    const float *ngram_data, int ngram_size, const WordSet *ws, double wword,
    double *indiv_buf, double *delta_buf, int *order_out) {

    if (R < 1) return 0.0;
    if (R > HELD_KARP_MAX_NODES) {        // best-L not attempted above the cap
        for (int i = 0; i < R; i++) order_out[i] = i;
        double s = 0.0;
        for (int i = 0; i < R; i++) s += row_objective(rows[i], rowlen[i], ngram_data, ngram_size, ws, wword);
        return s;
    }

    for (int a = 0; a < R; a++)
        indiv_buf[a] = row_objective(rows[a], rowlen[a], ngram_data, ngram_size, ws, wword);

    // seam delta[a][b] = obj(row_a ++ row_b) - indiv[a] - indiv[b]; the join
    // creates exactly the windows / boundary tokens that straddle the a->b boundary.
    static int joinbuf[2 * MAX_CIPHER_LENGTH];
    for (int a = 0; a < R; a++) {
        for (int b = 0; b < R; b++) {
            if (a == b) { delta_buf[a * R + b] = 0.0; continue; }
            int n = 0;
            for (int k = 0; k < rowlen[a]; k++) joinbuf[n++] = rows[a][k];
            for (int k = 0; k < rowlen[b]; k++) joinbuf[n++] = rows[b][k];
            double j = row_objective(joinbuf, n, ngram_data, ngram_size, ws, wword);
            delta_buf[a * R + b] = j - indiv_buf[a] - indiv_buf[b];
        }
    }
    return held_karp_best_path(R, indiv_buf, delta_buf, order_out);
}



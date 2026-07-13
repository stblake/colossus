#include "columnar_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Dedicated columnar transposition solver (TRANSCOL / TRANSCOL2)
//
//  Where solve_general_transposition hill-climbs the full N-length
//  permutation key and leans on a structure-score guard, this solver
//  optimizes only the small per-stage column-order permutation (length
//  K, the column count). The search space is K! (or K1!*K2! for double),
//  every candidate is a genuine columnar layout by construction, and a
//  column swap maps one-to-one onto the key -- the AZDecrypt approach.
// =====================================================================

// Decrypt one or two stacked columnar stages. cipher[] -> out[].
// Encryption applies stage 0 then stage 1, so we invert stage 1 then stage 0.
static void decrypt_columnar_stages(int cipher[], int len, int nstages,
    int K[2], int order[2][MAX_COLS], int dir[2], int out[]) {
    if (nstages == 1) {
        decrypt_columnar(cipher, len, K[0], order[0], dir[0], out);
    } else {
        int tmp[MAX_CIPHER_LENGTH];
        decrypt_columnar(cipher, len, K[1], order[1], dir[1], tmp);  // undo outer (2nd applied)
        decrypt_columnar(tmp,    len, K[0], order[0], dir[0], out);  // undo inner (1st applied)
    }
}

// Perturb a single stage's column order with one permutation-preserving move
// (swap dominant, plus short reverse and short block-move). The direction is
// only flipped when the user asked to search both (search_dir == COL_READ_BOTH);
// otherwise it stays pinned to the requested read direction.
static void perturbate_column_order(int order[], int K, int *dir, int search_dir) {
    if (search_dir == COL_READ_BOTH && frand() < 0.05) {
        *dir = 1 - *dir;
        return;
    }
    if (K < 2) return;

    double r = frand();
    if (r < 0.70) {
        // Swap two columns (the dominant move).
        int a = rand_int(0, K), b = rand_int(0, K);
        int t = order[a]; order[a] = order[b]; order[b] = t;
    } else if (r < 0.85) {
        // Reverse a short segment of the order.
        int max_blk = min(K, 8);
        int blk = rand_int(2, max_blk + 1);
        int s = rand_int(0, K - blk + 1);
        for (int a = s, b = s + blk - 1; a < b; a++, b--) {
            int t = order[a]; order[a] = order[b]; order[b] = t;
        }
    } else {
        // Cut a short block and re-insert it elsewhere (a range rotation).
        int max_blk = min(K, 8);
        int blk = rand_int(1, max_blk + 1);
        int s = rand_int(0, K - blk + 1);
        int d = rand_int(0, K - blk + 1);
        if (d == s) return;
        int tmp[8];
        for (int a = 0; a < blk; a++) tmp[a] = order[s + a];
        if (d < s) {
            for (int a = s - 1; a >= d; a--) order[a + blk] = order[a];
        } else {
            for (int a = s + blk; a < d + blk; a++) order[a - blk] = order[a];
        }
        for (int a = 0; a < blk; a++) order[d + a] = tmp[a];
    }
}

// Seed a fresh restart: random column order(s), direction(s) per read_direction.
static void columnar_seed(ColossusConfig *cfg, int nstages,
    int K[2], int order[2][MAX_COLS], int dir[2]) {
    for (int s = 0; s < nstages; s++) {
        for (int c = 0; c < K[s]; c++) order[s][c] = c;
        shuffle(order[s], K[s]);
        dir[s] = (cfg->read_direction == COL_READ_BOTH) ? rand_int(0, 2) : cfg->read_direction;
    }
}

// ---- columnar model (TRANSCOL / TRANSCOL2; cipher-agnostic engine) ---------
// SHAPE_ANNEAL over the per-stage column order. The state lives in SolverState:
//   aux[0] = nstages; aux[1..2] = K per stage; aux[3..4] = read direction per
//   stage; key[] reinterpreted as order[2][MAX_COLS] (stride MAX_COLS).
// Single columnar enumerates one config per column count K (the K-sweep); double
// columnar is a single config that randomises (K1,K2) per restart in seed().

// Column-count search range, clamped to [2, len/2] and the array bound.
static void columnar_krange(const SolverCtx *ctx, int *lo, int *hi) {
    int l = ctx->cfg->min_cols, h = ctx->cfg->max_cols;
    int cap = ctx->cipher_len / 2;
    if (cap < 2) cap = 2;
    if (cap > MAX_COLS) cap = MAX_COLS;
    if (l < 2) l = 2;
    if (h > cap) h = cap;
    if (l > h) l = h;
    *lo = l; *hi = h;
}

static int columnar_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    if (ctx->cfg->cipher_type == TRANSCOL2) {
        out[0].period = 0; out[0].j = 0; out[0].k = 0;
        return 1;                                   // (K1,K2) randomised per restart
    }
    int lo, hi;
    columnar_krange(ctx, &lo, &hi);
    int n = 0;
    for (int K = lo; K <= hi && n < cap; K++) {      // one config per column count
        out[n].period = K; out[n].j = 0; out[n].k = 0;
        n++;
    }
    return n;
}

static void columnar_model_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    int nstages = (ctx->cfg->cipher_type == TRANSCOL2) ? 2 : 1;
    int *K   = &st->aux[1];
    int *dir = &st->aux[3];
    int (*order)[MAX_COLS] = (int (*)[MAX_COLS]) st->key;
    st->aux[0] = nstages;
    if (nstages == 1) {
        K[0] = cc->period;
    } else {
        int lo, hi;
        columnar_krange(ctx, &lo, &hi);
        K[0] = rand_int(lo, hi + 1);
        K[1] = rand_int(lo, hi + 1);
    }
    columnar_seed(ctx->cfg, nstages, K, order, dir);
}

static void columnar_model_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                                   bool *force_primary) {
    (void) cc; (void) force_primary;
    int nstages = st->aux[0];
    int (*order)[MAX_COLS] = (int (*)[MAX_COLS]) st->key;
    int s = (nstages == 1) ? 0 : rand_int(0, 2);
    int dir = st->aux[3 + s];
    perturbate_column_order(order[s], st->aux[1 + s], &dir, ctx->cfg->read_direction);
    st->aux[3 + s] = dir;
}

static void columnar_model_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int nstages = src->aux[0];
    const int (*sord)[MAX_COLS] = (const int (*)[MAX_COLS]) src->key;
    int (*dord)[MAX_COLS] = (int (*)[MAX_COLS]) dst->key;
    dst->aux[0] = nstages;
    for (int s = 0; s < nstages; s++) {
        dst->aux[1 + s] = src->aux[1 + s];          // K
        dst->aux[3 + s] = src->aux[3 + s];          // dir
        int K = src->aux[1 + s];
        for (int c = 0; c < K; c++) dord[s][c] = sord[s][c];
    }
}

static void columnar_model_decrypt(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                                   int *out, double *score_adjust) {
    (void) cc; (void) score_adjust;
    int nstages = st->aux[0];
    int K[2]   = { st->aux[1], st->aux[2] };
    int dir[2] = { st->aux[3], st->aux[4] };
    int (*order)[MAX_COLS] = (int (*)[MAX_COLS]) st->key;
    decrypt_columnar_stages(ctx->cipher, ctx->cipher_len, nstages, K, order, dir, out);
}

static void columnar_model_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    int nstages = st->aux[0];
    int K[2]   = { st->aux[1], st->aux[2] };
    int dir[2] = { st->aux[3], st->aux[4] };
    int (*order)[MAX_COLS] = (int (*)[MAX_COLS]) st->key;
    int buf[MAX_CIPHER_LENGTH];
    decrypt_columnar_stages(ctx->cipher, ctx->cipher_len, nstages, K, order, dir, buf);

    double elapsed = engine_elapsed_sec(stats);
    double n_iter_per_sec = (elapsed > 0.) ? ((double) stats->n_iterations)/elapsed : 0.;
    printf("\n%.2f\t[sec]\n", elapsed);
    printf("%.0fK\t[it/sec/thread]\n", 1.e-3*n_iter_per_sec);
    printf("%d\t[restarts]\n", stats->n_restarts);
    printf("%d\t[backtracks]\n", stats->n_backtracks);
    printf("%d\t[slips]\n", stats->n_slips);
    printf("%.4f\t[entropy]\n", entropy(buf, ctx->cipher_len));
    printf("%.2f\t[score]\n", score);
    for (int s = 0; s < nstages; s++)
        printf("stage %d: K=%d dir=%s\t[params]\n", s + 1, K[s],
            dir[s] == COL_READ_BT ? "bt" : "tb");
    printf("\n");
    print_text(buf, ctx->cipher_len); printf("\n");
    fflush(stdout);
}

static void columnar_model_report(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    SharedData *shared = ctx->shared;
    int cipher_len = ctx->cipher_len;
    int *cipher_indices = ctx->cipher;
    char *cribtext_str = ctx->cribtext;
    int n_cribs = ctx->n_cribs;
    int nstages = st->aux[0];
    int best_K[2]   = { st->aux[1], st->aux[2] };
    int best_dir[2] = { st->aux[3], st->aux[4] };
    int (*best_order)[MAX_COLS] = (int (*)[MAX_COLS]) st->key;
    int n_words_found = 0;

    if (nstages == 1)
        printf("\ntranscol: single columnar, %d columns, read %s\n",
            best_K[0], best_dir[0] == COL_READ_BT ? "bottom-to-top" : "top-to-bottom");
    else
        printf("\ntranscol2: double columnar, %d x %d columns\n", best_K[0], best_K[1]);

    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < cipher_len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[cipher_len] = '\0';

    if (cfg->dictionary_present && shared->dict != NULL) {
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);
    }

    printf("\nResult Score: %.2f | Words: %d\n", score, n_words_found);

    print_text(cipher_indices, cipher_len);
    printf("\n");
    print_text(decrypted, cipher_len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, cipher_len);
    printf("%s\n", cribtext_str);

    if (PARTIAL_CRIB_MATCH && n_cribs > 0) {
        for (int i = 0; i < cipher_len; i++) {
            if (cribtext_str[i] == '_') {
                printf("_");
            } else {
                int diff = abs(decrypted[i] - (g_char_to_idx[toupper((unsigned char)cribtext_str[i]) & 127]));
                if (diff < 10) printf("%d", diff); else printf("*");
            }
        }
    }
    printf("\n");

    // Recovered column order(s), for reproduction.
    for (int s = 0; s < nstages; s++) {
        printf("stage %d (K=%d, dir=%s) order:", s + 1, best_K[s],
            best_dir[s] == COL_READ_BT ? "bt" : "tb");
        for (int c = 0; c < best_K[s]; c++) printf(" %d", best_order[s][c]);
        printf("\n");
    }
    printf("\n");

    // One-liner summary.
    if (cfg->dictionary_present) {
        printf(">>> %.2f, %d, %d, ", score, n_words_found, cfg->cipher_type);
    } else {
        printf(">>> %.2f, %d, ", score, cfg->cipher_type);
    }
    if (nstages == 1) printf("%d, ", best_K[0]);
    else printf("%d, %d, ", best_K[0], best_K[1]);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(cipher_indices, cipher_len);
    printf(", ");
    print_text(decrypted, cipher_len);
    printf("\n");
}

static const CipherModel COLUMNAR_MODEL = {
    .name = "columnar",
    .shape = SHAPE_ANNEAL,
    .needs_hist = false,
    .enumerate_configs = columnar_enumerate,
    .key_len = NULL,
    .seed = columnar_model_seed,
    .perturb = columnar_model_perturb,
    .copy_state = columnar_model_copy,
    .decrypt = columnar_model_decrypt,
    .report = columnar_model_report,
    .report_verbose = columnar_model_report_verbose,
};

void solve_columnar(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str; // ciphertext is carried as cipher_indices.

    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a columnar solve.\n\n");
        return ;
    }

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    run_solver(&COLUMNAR_MODEL, &ctx);
}




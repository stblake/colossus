#include "permutation_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"


// General transposition solver (AZDecrypt-style)
//
// Solves an arbitrary columnar / route transposition by hill-climbing the full
// permutation key directly, rather than the fixed parameters of a specific
// transform family. The key is an array key[i] = source position, so the
// candidate plaintext is decrypted[i] = ciphertext[key[i]]. Neighbour moves
// (swap, segment reverse, block move) all preserve the permutation; scoring,
// slip, backtracking and restarts reuse the program's standard machinery.

// Apply a permutation key: decrypted[i] = cipher[key[i]].
static void apply_permutation(int cipher_indices[], int key[], int len, int decrypted[]) {
    for (int i = 0; i < len; i++) decrypted[i] = cipher_indices[key[i]];
}

// Structural regularity of a permutation key, in [0,1]. Columnar / route
// transpositions are periodic: for the period p (the column count), plaintext
// positions i and i+p sit in the same column one row apart, so their ciphertext
// positions differ by a constant stride. We scan candidate periods and, for
// each, histogram the strides key[i+p]-key[i]; the best (most concentrated)
// period gives the score = fraction of pairs in the single modal stride.
// English-looking but structurally random permutations score low; a true
// (keyed) columnar key scores high at its period. This is the n-gram-gaming
// guard that AZDecrypt gets from its periodic-redundancy rule.
static double key_structure_score(int key[], int len, int *out_period) {
    static int hist[2 * MAX_CIPHER_LENGTH];
    if (out_period) *out_period = 0;
    if (len < 2) return 0.0;

    int max_period = len / 3;
    if (max_period > 24) max_period = 24;   // bound cost; covers realistic column counts
    if (max_period < 1) max_period = 1;

    double best_frac = 0.0;
    for (int p = 1; p <= max_period; p++) {
        int npairs = len - p;
        if (npairs < 1) break;
        int modal = 0;
        for (int i = 0; i + p < len; i++) {
            int b = key[i + p] - key[i] + len;   // offset into [0, 2*len)
            int c = ++hist[b];
            if (c > modal) modal = c;
        }
        // Re-zero only the bins we touched (cheaper than a full memset per period).
        for (int i = 0; i + p < len; i++) hist[key[i + p] - key[i] + len] = 0;

        double frac = (double) modal / (double) npairs;
        if (frac > best_frac) { best_frac = frac; if (out_period) *out_period = p; }
    }
    return best_frac;
}

// Perturb a permutation key with one of several permutation-preserving moves.
// target_period is the column count currently detected in the key (0 if none);
// the column-swap move uses it so it reorders the actual columns rather than
// guessing a period at random.
static void perturbate_permutation(int key[], int len, int target_period) {
    int max_p = min(len / 2, 40);

    // Column swap: swap the two sets of key entries {.., r*p+c1, ..} and
    // {.., r*p+c2, ..}. With p = the detected column count this reorders whole
    // columns in one step without disturbing the within-column structure — the
    // move that lets a columnar-seeded key reach the exact column order. Once a
    // period is detected it is the dominant move; small position moves would
    // need ~one-column-height steps to do the same and erode the structure.
    bool do_colswap = (target_period >= 2 && target_period <= max_p && max_p >= 2 && frand() < 0.6);
    if (do_colswap) {
        int p = target_period;
        int c1 = rand_int(0, p), c2 = rand_int(0, p);
        if (c1 == c2) return;
        for (int r = 0; ; r++) {
            int i1 = r * p + c1, i2 = r * p + c2;
            if (i1 >= len || i2 >= len) break;   // only swap rows where both columns exist
            int t = key[i1]; key[i1] = key[i2]; key[i2] = t;
        }
        return;
    }

    int move = rand_int(0, 3);
    if (move == 0) {
        // Swap two positions.
        int i = rand_int(0, len), j = rand_int(0, len);
        int t = key[i]; key[i] = key[j]; key[j] = t;
    } else if (move == 1) {
        // Reverse a short segment.
        int max_blk = min(len, 14);
        int blk = rand_int(2, max_blk + 1);
        int s = rand_int(0, len - blk + 1);
        for (int a = s, b = s + blk - 1; a < b; a++, b--) {
            int t = key[a]; key[a] = key[b]; key[b] = t;
        }
    } else {
        // Cut a short block and re-insert it elsewhere (a range rotation).
        int max_blk = min(len, 14);
        int blk = rand_int(1, max_blk + 1);
        int s = rand_int(0, len - blk + 1);          // block start
        int d = rand_int(0, len - blk + 1);          // destination start
        if (d == s) return;
        int tmp[14];
        for (int a = 0; a < blk; a++) tmp[a] = key[s + a];
        if (d < s) {
            for (int a = s - 1; a >= d; a--) key[a + blk] = key[a];
        } else {
            for (int a = s + blk; a < d + blk; a++) key[a - blk] = key[a];
        }
        for (int a = 0; a < blk; a++) key[d + a] = tmp[a];
    }
}

// Seed a key with a columnar-transposition layout: the ciphertext is read as
// `p` columns (row-major fill, so the final row may be short) taken in column
// order `ord`. key[plaintext_pos] = ciphertext_pos. Seeding restarts from such
// structured keys turns the intractable free-permutation search into "find the
// right period and column order", which the climber refines easily.
static void build_columnar_seed(int key[], int len, int p, int ord[]) {
    int pos = 0;
    for (int k = 0; k < p; k++) {
        int c = ord[k];
        for (int r = 0; r * p + c < len; r++) key[r * p + c] = pos++;
    }
}

// ---- general transposition model (TYPE transposition; cipher-agnostic) -----
// SHAPE_ANNEAL over the full length-N permutation key (st->key, key_len = N). The
// AZDecrypt periodic-redundancy guard (key_structure_score) is folded in two
// ways: its value becomes the decrypt score_adjust, and the period it detects is
// carried in st->aux[0] so the next perturb's column-swap move can target it.
// The RNG draw order matches the original climber, so the search is bit-identical.

static int permutation_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx; (void) cap;
    out[0].period = 0; out[0].j = 0; out[0].k = 0;
    return 1;
}

static void permutation_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    int len = ctx->cipher_len;
    st->key_len = len;
    // Mostly a structured columnar seed (random period and column order),
    // occasionally a fully random permutation so route transpositions stay reachable.
    if (frand() < 0.85) {
        int max_p = min(len / 2, 40);
        if (max_p < 2) max_p = 2;
        int p = rand_int(2, max_p + 1);
        int ord[64];
        if (p > 64) p = 64;
        for (int c = 0; c < p; c++) ord[c] = c;
        shuffle(ord, p);
        build_columnar_seed(st->key, len, p, ord);
    } else {
        for (int i = 0; i < len; i++) st->key[i] = i;
        shuffle(st->key, len);
    }
    st->aux[0] = 0;     // detected period; decrypt() fills it from key_structure_score
}

static void permutation_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                                bool *force_primary) {
    (void) cc; (void) force_primary;
    perturbate_permutation(st->key, ctx->cipher_len, st->aux[0]);
}

static void permutation_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int len = src->key_len;
    dst->key_len = len;
    dst->aux[0] = src->aux[0];
    for (int i = 0; i < len; i++) dst->key[i] = src->key[i];
}

static void permutation_decrypt(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                                int *out, double *score_adjust) {
    (void) cc;
    apply_permutation(ctx->cipher, st->key, ctx->cipher_len, out);
    int period = 0;
    double s = key_structure_score(st->key, ctx->cipher_len, &period);
    st->aux[0] = period;                                     // carry period for the next perturb
    *score_adjust = ctx->cfg->weight_structure * s;
}

static void permutation_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    int buf[MAX_CIPHER_LENGTH];
    apply_permutation(ctx->cipher, (int *) st->key, ctx->cipher_len, buf);
    double elapsed = engine_elapsed_sec(stats);
    double n_iter_per_sec = (elapsed > 0.) ? ((double) stats->n_iterations)/elapsed : 0.;
    printf("\n%.2f\t[sec]\n", elapsed);
    printf("%.0fK\t[it/sec/thread]\n", 1.e-3*n_iter_per_sec);
    printf("%d\t[restarts]\n", stats->n_restarts);
    printf("%d\t[backtracks]\n", stats->n_backtracks);
    printf("%d\t[slips]\n", stats->n_slips);
    printf("%.4f\t[entropy]\n", entropy(buf, ctx->cipher_len));
    printf("%.2f\t[score]\n", score);
    printf("%d\t[period]\n", st->aux[0]);
    printf("\n");
    print_text(buf, ctx->cipher_len); printf("\n");
    print_solution_check(buf, ctx->cipher_len);
    fflush(stdout);
}

static void permutation_report(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    SharedData *shared = ctx->shared;
    int cipher_len = ctx->cipher_len;
    int *cipher_indices = ctx->cipher;
    char *cribtext_str = ctx->cribtext;
    int n_cribs = ctx->n_cribs;
    int n_words_found = 0;

    printf("\ntransposition: permutation key of length %d\n", cipher_len);

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
    print_solution_check(decrypted, cipher_len);
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

    // Recovered key (source position for each output position), for reproduction.
    printf("key: ");
    for (int i = 0; i < cipher_len; i++) printf("%d%s", st->key[i], (i + 1 < cipher_len) ? " " : "");
    printf("\n\n");

    // One-liner summary.
    if (cfg->dictionary_present) {
        printf(">>> %.2f, %d, %d, %s, ", score, n_words_found, cfg->cipher_type,
            cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    } else {
        printf(">>> %.2f, %d, %s, ", score, cfg->cipher_type,
            cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    }
    print_text(cipher_indices, cipher_len);
    printf(", ");
    print_text(decrypted, cipher_len);
    printf("\n");
}

static const CipherModel PERMUTATION_MODEL = {
    .name = "transposition",
    .shape = SHAPE_ANNEAL,
    .needs_hist = false,
    .enumerate_configs = permutation_enumerate,
    .key_len = NULL,
    .seed = permutation_seed,
    .perturb = permutation_perturb,
    .copy_state = permutation_copy,
    .decrypt = permutation_decrypt,
    .report = permutation_report,
    .report_verbose = permutation_report_verbose,
};

void solve_general_transposition(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs) {

    (void) ciphertext_str; // ciphertext is carried as cipher_indices.

    if (cipher_len < 3) {
        printf("\n\nERROR: ciphertext too short for a transposition solve.\n\n");
        return ;
    }

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    run_solver(&PERMUTATION_MODEL, &ctx);
}



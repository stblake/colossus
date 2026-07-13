#include "hill_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Hill solver (TYPE hill)
// =====================================================================
//
// The Hill cipher multiplies each block of k plaintext letters by a k x k key matrix
// mod 26. Breaking it ciphertext-only is two coupled problems: recover the block size k
// and recover the matrix. There is no IoC-style estimator for k, so the solver simply
// sweeps candidate sizes (one engine config per k in [2 .. max]) and lets the n-gram
// score pick the winner -- a wrong k decrypts to gibberish and loses. -period pins k.
//
// For each k the engine hill-climbs / anneals the matrix. The crucial trick: the state
// carried in st->key IS the DECRYPTION matrix, applied straight to the ciphertext
// (decrypted = D * cipher mod 26), so the hot path never inverts a matrix. The true
// plaintext came from an invertible encryption key K, so the climb converges on D = K^-1.
// The matrix is inverted only once, at report time, to display the recovered encryption
// key (cosmetic -- the graded output is the recovered plaintext).
//
// A matrix is a bijection iff it is invertible mod 26. A SINGULAR D is NOT a bijection:
// it folds the ciphertext onto a sub-lattice and decrypts to a low-entropy repetitive
// string (the zero matrix -> all A's) that OUT-scores real plaintext on n-grams, so the
// climb is attracted into that collapse rather than out-competing it. The decrypt hook
// therefore penalises singular candidates via score_adjust (HILL_SINGULAR_PENALTY) -- the
// Hill analogue of the homophonic anti-collapse term. Like Playfair/Bifid/Trifid, Hill
// effectively needs the discriminating -logprob (mean log-probability) fitness.

#define HILL_DEFAULT_MAX_K 5

// A singular decryption matrix (det not coprime to 26) is NOT a bijection: it folds the
// ciphertext onto a sub-lattice and decrypts to a low-entropy, repetitive plaintext (the
// zero matrix -> all A's; a rank-deficient matrix -> NANANA...) that out-scores real
// plaintext on raw n-grams. The climb is therefore attracted into the collapse basin
// rather than repelled from it. Real Hill decryption matrices are invertible, so any
// singular candidate is penalised decisively below the worst plausible plaintext score
// (mean log-prob lives in ~[-8,-2]); the row-add perturbation preserves the determinant,
// so the search still moves freely within the invertible subset.
#define HILL_SINGULAR_PENALTY (-1000.0)

typedef struct {
    int k_list[HILL_MAX_K];   // candidate block sizes (config order)
    int n_k;                  // number of candidate block sizes
} HillScratch;

// One config per candidate block size; period carries k, key carries the k*k matrix.
static int hill_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const HillScratch *h = (const HillScratch *) ctx->model_scratch;
    int n = h->n_k;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = h->k_list[i];
        out[i].j = 0; out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

// Seed: a uniformly random decryption matrix (every entry in 0..25). Invertibility is not
// required for a seed -- the climb finds a good (and ultimately invertible) matrix.
static void hill_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int k = cc->period, n = k * k;
    for (int i = 0; i < n; i++) st->key[i] = rand_int(0, ALPHABET_SIZE);
    st->key_len = n;
}

// Neighbour move on the matrix: 85% change one element to a *different* random value
// (the dominant fine move); 10% randomize a whole row; 5% add a random nonzero multiple
// of one row to another mod 26 (a coarse, row-coupling jump that escapes basins the
// single-element move cannot). For k == 1 only the single-element move applies.
static void hill_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                         SolverState *st, bool *force_primary) {
    (void) ctx; (void) force_primary;
    int k = cc->period, n = k * k;
    double r = frand();
    if (r < 0.85 || k == 1) {                      // change one element (never a no-op)
        int idx = rand_int(0, n);
        int cur = st->key[idx];
        int nv = rand_int(0, ALPHABET_SIZE - 1);   // 0..24, then skip `cur` -> 0..25 \ {cur}
        if (nv >= cur) nv++;
        st->key[idx] = nv;
    } else if (r < 0.95) {                          // randomize a whole row
        int row = rand_int(0, k);
        for (int c = 0; c < k; c++) st->key[row * k + c] = rand_int(0, ALPHABET_SIZE);
    } else {                                        // row_dst += t * row_src  (mod 26)
        int rs = rand_int(0, k), rd = rand_int(0, k);
        if (rs == rd) rd = (rd + 1) % k;
        int t = 1 + rand_int(0, ALPHABET_SIZE - 1);
        for (int c = 0; c < k; c++) {
            int v = st->key[rd * k + c] + t * st->key[rs * k + c];
            st->key[rd * k + c] = ((v % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
        }
    }
}

static void hill_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    // Copy a fixed upper bound (k*k <= HILL_MAX_KEY) so the engine need not know k.
    for (int i = 0; i < HILL_MAX_KEY; i++) dst->key[i] = src->key[i];
}

static void hill_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, int *out, double *score_adjust) {
    hill_mat_mul_blocks(st->key, cc->period, ctx->cipher, ctx->cipher_len, out);
    int det = hill_det_mod(st->key, cc->period);
    *score_adjust = (hill_mod_inverse(det, ALPHABET_SIZE) == 0) ? HILL_SINGULAR_PENALTY : 0.0;
}

// Render a k x k matrix as an indented box of mod-26 integers.
static void hill_print_matrix(const int mat[], int k, const char *indent) {
    for (int r = 0; r < k; r++) {
        printf("%s", indent);
        for (int c = 0; c < k; c++) printf("%2d ", mat[r * k + c]);
        printf("\n");
    }
}

static void hill_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    printf("\n  k=%d, decryption matrix:\n", cc->period);
    hill_print_matrix(st->key, cc->period, "    ");
    report_transposition_verbose(ctx, score, decrypted, stats, "hill");
}

static void hill_report(const SolverCtx *ctx, const SolverConfig *cc,
                        const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    int len = ctx->cipher_len, k = cc->period, n = k * k;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[len] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // The state IS the decryption matrix; invert it to recover the encryption key.
    int enc[HILL_MAX_KEY];
    int invertible = hill_mat_inverse(st->key, k, enc);

    // Compact letter renderings (A..Z, row-major) for the one-line summary.
    char dmat[HILL_MAX_KEY + 1], emat[HILL_MAX_KEY + 1];
    for (int i = 0; i < n; i++) dmat[i] = index_to_char(st->key[i]);
    dmat[n] = '\0';
    if (invertible) {
        for (int i = 0; i < n; i++) emat[i] = index_to_char(enc[i]);
        emat[n] = '\0';
    } else {
        strcpy(emat, "SINGULAR");
    }

    printf("\nResult Score: %.2f | Words: %d | k=%d | decrypt-matrix=%s | encrypt-key=%s\n",
        score, n_words_found, k, dmat, emat);

    print_cipher(ctx->cipher, len, NULL);
    printf("\n");
    print_text(decrypted, len);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, len);
    printf("%s\n", ctx->cribtext);

    printf("\nrecovered %dx%d decryption matrix (row major, mod 26):\n", k, k);
    hill_print_matrix(st->key, k, "    ");
    if (invertible) {
        printf("encryption key = (decryption matrix)^-1 mod 26:\n");
        hill_print_matrix(enc, k, "    ");
    } else {
        printf("(decryption matrix is singular mod 26 -- no encryption key)\n");
    }

    // Publish the recovered solution for callers that pass a SolveResult.
    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = k;       // report the recovered block size here
        vec_copy(decrypted, ctx->result->decrypted, len);
        ctx->result->decrypted_len = len;
    }

    // One-liner summary: >>> score, [words,] type, k=, dmat=, emat=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, k=%d, dmat=%s, emat=%s, ",
            score, n_words_found, cfg->cipher_type, k, dmat, emat);
    else
        printf(">>> %.2f, %d, k=%d, dmat=%s, emat=%s, ",
            score, cfg->cipher_type, k, dmat, emat);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, len, NULL);
    printf(", ");
    print_text(decrypted, len);
    printf("\n");
}

static const CipherModel HILL_MODEL = {
    .name = "hill", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = hill_enumerate, .key_len = NULL,
    .seed = hill_seed, .perturb = hill_perturb, .copy_state = hill_copy,
    .decrypt = hill_decrypt_hook, .report = hill_report,
    .report_verbose = hill_report_verbose,
};

void solve_hill(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    // Hill is mod-26 by construction; it needs the full 26-letter alphabet.
    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Hill needs the 26-letter alphabet (got %d). Do not exclude letters.\n\n",
               g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Hill solve.\n\n");
        return;
    }
    // Hill letters must be solid: reject any sentinel (space/punctuation).
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Hill ciphertext must be solid letters (try -skipspaces).\n\n", i);
            return;
        }

    HillScratch scratch;
    if (cfg->period_present) {
        int k = cfg->period;
        if (k < 1 || k > HILL_MAX_K) {
            printf("\n\nERROR: Hill block size k=%d out of range [1..%d].\n\n", k, HILL_MAX_K);
            return;
        }
        scratch.k_list[0] = k;
        scratch.n_k = 1;
        if (cfg->verbose) printf("\nhill: block size pinned to k=%d\n", k);
    } else {
        int max_k = (cfg->max_period > 0) ? cfg->max_period : HILL_DEFAULT_MAX_K;
        if (max_k > HILL_MAX_K) max_k = HILL_MAX_K;
        if (max_k < 2) max_k = 2;
        int nk = 0;
        for (int k = 2; k <= max_k; k++)
            if (cipher_len >= 2 * k) scratch.k_list[nk++] = k;   // need a block to repeat
        if (nk == 0) scratch.k_list[nk++] = 2;
        scratch.n_k = nk;
    }

    if (cfg->verbose)
        printf("\nhill: %d positions, 26-letter alphabet, %d candidate block size(s)\n",
            cipher_len, scratch.n_k);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;          // hill_report fills it (may be NULL for CLI use)

    run_solver(&HILL_MODEL, &ctx);
}

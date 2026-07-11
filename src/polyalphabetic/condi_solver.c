#include "condi_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Condi solver (TYPE condi)
// =====================================================================
//
// Condi (condi.c) enciphers over a keyed alphabet sigma (a 26-permutation) with a PLAINTEXT-
// FEEDBACK running key: the shift for each letter is the position of the preceding plaintext
// letter in sigma; the first letter uses a starter offset. idx(ct_i) = (idx(pt_i) + idx(pt_{i-1})
// + 1) mod 26, so a given plaintext letter enciphers under a position-dependent shift -> the
// cipher is effectively polyalphabetic (flatter monograms than a simple substitution), and the
// only key is sigma (26!) plus the starter (0..25).
//
// There is NO period (a feedback stream cipher) and NO per-column / monogram decoupling: the
// running key is plaintext-derived, hence entangled with sigma (unlike Gromark's primer, which
// fixes the running key independent of sigma and admits a Hungarian pre-pass). The attack shipped
// here is a FREE-permutation sigma anneal (random shuffle seed, cell-swap moves, no anti-collapse
// penalty since every sigma is a bijection) with the starter searched: the 26 starter values are
// ENUMERATED as engine configs (a wrong starter cascades the causal decrypt into gibberish -> low
// n-gram score) and the n-gram score across configs picks the winner. -startkey pins a single
// starter. sigma is a free permutation (not a keyed keyword+tail) because the ACA alphabet may be
// cyclically SHIFTED, and a rotation genuinely changes the plaintext (idx(ct) shifts by r but the
// feedback sum by 2r, so no cancellation) -- the free permutation absorbs the rotation for free.
//
// !!! DOCUMENTED STRUCTURAL LIMITATION (see tests/test_condi_solver.c) !!!
// The plaintext feedback makes local search INEFFECTIVE for Condi. Because idx(pt_i) feeds the next
// offset, one sigma cell swap re-derives the ENTIRE downstream plaintext (an alternating-sign shift
// cascades from the first changed position), so the true sigma is an ISOLATED NEEDLE: measured on a
// real solve, the true key scores ~3.40 while its BEST single-swap neighbour scores only ~2.48 and
// the mean neighbour ~1.37 (the random floor is ~1.01) -- there is NO basin, hence no gradient for
// anneal / shotgun / PSO, blind OR crib-assisted. So this solver does NOT recover Condi blind at any
// budget; it is kept as an honest bounded attempt (and the substrate for a future crib-anchored
// constraint solver, the tractable attack: a crib gives linear congruences pos[ct_i] = pos[pt_i] +
// pos[pt_{i-1}] + 1 mod 26 that pin sigma directly). Cribs are still passed through (positional:
// decrypted[i] == plaintext[i]) and -logprob is honoured.

#define CONDI_ALPHA  ALPHABET_SIZE      // 26-letter keyed alphabet

typedef struct {
    int  n;                              // cipher length (== plaintext length)
    int  n_configs;                      // number of starter configs
    const int *starters;                 // [n_configs]  starter offset 0..25 per config
} CondiScratch;

// Final config store (single-threaded), filled by solve_condi before run_solver.
static int g_condi_starters[CONDI_ALPHA];

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int condi_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const CondiScratch *a = (const CondiScratch *) ctx->model_scratch;
    int n = a->n_configs;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = 0;
        out[i].j = i; out[i].k = 0;
        out[i].aux[0] = a->starters[i];   // starter offset for this config
        out[i].aux[1] = 0;
    }
    return n;
}

static void condi_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    for (int i = 0; i < CONDI_ALPHA; i++) st->key[i] = i;
    shuffle(st->key, CONDI_ALPHA);        // random keyed alphabet
    st->aux[0] = cc->aux[0];              // starter (pinned per config, never perturbed)
    st->key_len = CONDI_ALPHA;
}

// Simple-substitution swap on sigma (a 26-permutation with no grid geometry -> plain cell swaps).
static void condi_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                          SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int a = rand_int(0, CONDI_ALPHA), b = rand_int(0, CONDI_ALPHA);
    int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
}

static void condi_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < CONDI_ALPHA; i++) dst->key[i] = src->key[i];
    dst->aux[0] = src->aux[0];
    dst->key_len = src->key_len;
}

static void condi_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                               SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    int n = ctx->cipher_len;
    int sinv[CONDI_ALPHA];
    for (int i = 0; i < CONDI_ALPHA; i++) sinv[st->key[i]] = i;
    condi_decrypt(ctx->cipher, n, st->key, sinv, st->aux[0], out);
    *score_adjust = 0.0;
}

// ===================================================================
//  Reporting
// ===================================================================

static void condi_alpha_string(const int sigma[], char out[]) {
    for (int i = 0; i < CONDI_ALPHA; i++) out[i] = index_to_char(sigma[i]);
    out[CONDI_ALPHA] = '\0';
}

static void condi_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    char alpha[CONDI_ALPHA + 1]; condi_alpha_string(st->key, alpha);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  starter %d, score=%.4f  [%.1fs, %d restarts]\n    alphabet=%s\n",
        st->aux[0], score, elapsed, stats->n_restarts, alpha);
    fflush(stdout);
}

static void condi_report(const SolverCtx *ctx, const SolverConfig *cc,
                         const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len, starter = st->aux[0];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char alpha[CONDI_ALPHA + 1]; condi_alpha_string(st->key, alpha);

    printf("\nResult Score: %.2f | Words: %d | starter=%d | alphabet=%s\n",
        score, n_words_found, starter, alpha);

    print_cipher(ctx->cipher, n, NULL);
    printf("\n");
    print_text(decrypted, n);
    printf("\n%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->condi_start = starter;
        for (int i = 0; i < CONDI_ALPHA; i++) ctx->result->plaintext_keyword[i] = st->key[i];
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, starter=, alphabet=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, starter=%d, alphabet=%s, ",
            score, n_words_found, cfg->cipher_type, starter, alpha);
    else
        printf(">>> %.2f, %d, starter=%d, alphabet=%s, ",
            score, cfg->cipher_type, starter, alpha);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel CONDI_MODEL = {
    .name = "condi", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = condi_enumerate, .key_len = NULL,
    .seed = condi_seed, .perturb = condi_perturb, .copy_state = condi_copy,
    .decrypt = condi_decrypt_hook, .report = condi_report,
    .report_verbose = condi_report_verbose,
};

// ===================================================================
//  Entry point
// ===================================================================

void solve_condi(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Condi needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 8) {
        printf("\n\nERROR: ciphertext too short for a Condi solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Condi ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    // Starter configs: -startkey pins one, else enumerate 0..25 (the n-gram score selects it).
    int nc = 0;
    if (cfg->startkey_present) {
        g_condi_starters[nc++] = ((cfg->startkey % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    } else {
        for (int s = 0; s < ALPHABET_SIZE; s++) g_condi_starters[nc++] = s;
    }

    if (cfg->verbose)
        printf("\ncondi: %d letters, free-permutation sigma anneal, %d starter config(s)\n",
            cipher_len, nc);

    CondiScratch a;
    a.n = cipher_len;
    a.n_configs = nc;
    a.starters = g_condi_starters;

    // Cribs are supported: the cipher is positional, so crib positions map straight onto the
    // decrypted plaintext (decrypted[i] == plaintext[i]).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &a;
    ctx.result = result;

    run_solver(&CONDI_MODEL, &ctx);
}

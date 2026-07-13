#include "bazeries_solver.h"
#include "bazeries.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Bazeries solver (TYPE bazeries)
// =====================================================================
//
// Bazeries ("simple substitution plus transposition", ACA) is keyed by ONE number N < 10^6.
// That single number drives BOTH a digit-grouped reversal transposition (group sizes cycle
// through N's decimal digits) and a fixed monoalphabetic map between a column-major plaintext
// square and a keyed (N spelled out) row-major ciphertext square (see bazeries.c).
//
// Colossus is optimisation-only, so rather than add an exhaustive-search engine driver, the
// number is CLIMBED as an ordinary CipherModel state: the climbed lane is N's decimal digits
// (key[0..D-1], leading digit 1..9), with one engine config per digit count D in [minD..maxD]
// (the union over D = 1..6 covers all of 1..999999). Because N drives both stages, every
// candidate is a complete deterministic decryption -- nothing is left coupled the way ADFGVX's
// square and column order are -- but the digit landscape is rugged: flipping one digit changes
// the spelled square AND a transposition group size, so the n-gram score alone gives sparse
// gradient over a < 10^6 space.
//
// The key enabler is a SQUARE-QUALITY decoupling reward folded into score_adjust, the analog of
// ADFGVX's structural IoC term. The inverse substitution is monoalphabetic, and a transposition
// only permutes letters (it leaves the monogram multiset unchanged), so the decrypt's mean
// English-monogram fit depends only on the square (i.e. on which N is spelled), NOT on the
// transposition digits. Rewarding it pulls the digit climb toward numbers whose spelled square
// has the correct frequency assignment -- a smooth-ish prior over N -- after which the full
// n-gram score discriminates the exact number (right square AND right digit pattern). Like the
// other square types it effectively needs -logprob. Cribs are not used (the transposition
// scrambles plaintext positions). The recovered plaintext is what the regression suite checks.

// Weight of the square-quality monogram reward. The decrypt's mean monogram fit spans ~0.04
// (a wrong, frequency-scrambling square) to ~0.066 (English), so this scales that ~0.026 span
// into a nudge comparable to the n-gram contrast of a wrong square, while staying a near-constant
// offset once the square is right (so the n-gram score then decides the digit pattern). Tuned
// against tests/test_bazeries_solver.c.
#define BAZ_MONO_WEIGHT 20.0

typedef struct {
    int n;            // cipher length (== plaintext length)
    int minD, maxD;   // digit-count (D) sweep range, in [1 .. BAZERIES_MAX_DIGITS]
} BazScratch;

// --- digit lane <-> number helpers -------------------------------------------------

static long baz_state_number(const SolverState *st) {
    int D = st->aux[0];
    long N = 0;
    for (int i = 0; i < D; i++) N = N * 10 + st->key[i];
    return N;
}

// --- model hooks ------------------------------------------------------------------

static int baz_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const BazScratch *a = (const BazScratch *) ctx->model_scratch;
    int c = 0;
    for (int D = a->minD; D <= a->maxD && c < cap; D++) {     // one config per digit count
        out[c].period = D; out[c].j = 0; out[c].k = 0;
        out[c].aux[0] = 0; out[c].aux[1] = 0;
        c++;
    }
    return c;
}

static void baz_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int D = cc->period;
    st->key[0] = rand_int(1, 10);                            // leading digit 1..9
    for (int i = 1; i < D; i++) st->key[i] = rand_int(0, 10);
    st->aux[0] = D;
    st->key_len = D;
}

static void baz_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                        SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int D = st->aux[0];
    if (D == 1 || frand() < 0.85) {                          // change one digit to a new value
        int pos = rand_int(0, D), cur = st->key[pos], nv;
        do { nv = (pos == 0) ? rand_int(1, 10) : rand_int(0, 10); } while (nv == cur);
        st->key[pos] = nv;
    } else {                                                 // change two digits (a bigger jump)
        for (int t = 0; t < 2; t++) {
            int pos = rand_int(0, D);
            st->key[pos] = (pos == 0) ? rand_int(1, 10) : rand_int(0, 10);
        }
    }
}

static void baz_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int D = src->aux[0];
    for (int i = 0; i < D; i++) dst->key[i] = src->key[i];
    dst->aux[0] = D;
    dst->key_len = D;
}

static void baz_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    int n = ctx->cipher_len;
    long N = baz_state_number(st);
    bazeries_decrypt(ctx->cipher, n, N, out);
    // Square-quality reward: the decrypt's mean monogram fit is transposition-independent, so it
    // gives the digit climb a gradient toward the correct (spelled) square.
    double mono = 0.0;
    for (int i = 0; i < n; i++) mono += g_monograms[out[i]];
    *score_adjust = BAZ_MONO_WEIGHT * mono / (double) n;
}

// --- reporting --------------------------------------------------------------------

static void baz_square_string(long N, char out[]) {
    int square[BAZERIES_GRID];
    bazeries_build_square(N, square);
    for (int i = 0; i < BAZERIES_GRID; i++) out[i] = index_to_char(square[i]);
    out[BAZERIES_GRID] = '\0';
}

static void baz_keyword_string(long N, char out[]) {
    int kw[BAZERIES_MAX_SPELL];
    int m = bazeries_spell(N, kw);
    for (int i = 0; i < m; i++) out[i] = index_to_char(kw[i]);
    out[m] = '\0';
}

static void baz_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    long N = baz_state_number(st);
    char sq[BAZERIES_GRID + 1]; baz_square_string(N, sq);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  number=%ld (D=%d), score=%.4f  [%.1fs, %d restarts]\n    square=%s\n",
        N, st->aux[0], score, elapsed, stats->n_restarts, sq);
    fflush(stdout);
}

static void baz_report(const SolverCtx *ctx, const SolverConfig *cc,
                       const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len;
    long N = baz_state_number(st);

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char sq[BAZERIES_GRID + 1]; baz_square_string(N, sq);
    char kw[BAZERIES_MAX_SPELL + 1]; baz_keyword_string(N, kw);

    printf("\nResult Score: %.2f | Words: %d | number=%ld | keyword=%s | square=%s\n",
        score, n_words_found, N, kw, sq);
    print_cipher(ctx->cipher, n, NULL);
    printf("\n");
    print_text(decrypted, n);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, n);
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = st->aux[0];             // digit count
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, number=, square=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, number=%ld, square=%s, ",
            score, n_words_found, cfg->cipher_type, N, sq);
    else
        printf(">>> %.2f, %d, number=%ld, square=%s, ",
            score, cfg->cipher_type, N, sq);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel BAZERIES_MODEL = {
    .name = "bazeries", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = baz_enumerate, .key_len = NULL,
    .seed = baz_seed, .perturb = baz_perturb, .copy_state = baz_copy,
    .decrypt = baz_decrypt_hook, .report = baz_report,
    .report_verbose = baz_report_verbose,
};

// =====================================================================
//  Entry point
// =====================================================================

void solve_bazeries(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != BAZERIES_GRID) {
        printf("\n\nERROR: Bazeries needs the 25-letter J->I alphabet (got %d). "
               "Run -type bazeries so the alphabet is forced.\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Bazeries solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Bazeries ciphertext must be solid letters.\n\n", i);
            return;
        }

    // Cribs are not used (the transposition scrambles plaintext positions).
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    // Digit-count sweep: -period pins a single D; otherwise sweep 1..BAZERIES_MAX_DIGITS.
    int minD = 1, maxD = BAZERIES_MAX_DIGITS;
    if (cfg->period_present) {
        minD = maxD = cfg->period;
        if (minD < 1) minD = maxD = 1;
        if (maxD > BAZERIES_MAX_DIGITS) minD = maxD = BAZERIES_MAX_DIGITS;
    }

    BazScratch scratch = { .n = cipher_len, .minD = minD, .maxD = maxD };

    if (cfg->verbose)
        printf("\nbazeries: %d letters, key number digit-count sweep D in [%d..%d]\n",
            cipher_len, minD, maxD);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&BAZERIES_MODEL, &ctx);
}

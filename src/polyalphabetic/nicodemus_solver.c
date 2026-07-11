//
// Nicodemus solver
// ================
//
// Nicodemus is a substitution+transposition composite keyed by one keyword of length P:
// the keyword gives both the per-column Vigenere/Variant/Beaufort shifts AND the per-block
// columnar read order (its alphabetical rank order). The naive state -- a free keyword of P
// letters -- is rugged: one letter change moves a shift AND the order, and the n-gram score
// only rewards a key when order and all P shifts are simultaneously near-correct.
//
// Instead we decouple, exactly as the default -optimalcycle path does for the polyalphabetic
// ciphers: the annealed state is the COLUMN ORDER alone (a permutation of 0..P-1), and the P
// per-column shifts are DERIVED deterministically for each candidate order. After undoing the
// columnar (nicodemus_detranspose), the detransposed grid-column g = i % P is a Caesar sample
// under one shift; we pick the shift maximising that column's monogram fit against g_monograms.
// Because every column is fit to English MONOGRAMS regardless of order, the order search is
// driven by the n-gram (quadgram) score -- cross-column digraphs only form at the true order.
// Solving this general (order, shifts) form also cracks the ACA cipher (whose true key is the
// special case order == argsort(shifts)); the recovered plaintext is what matters.
//
// One engine config per (period P, block height H) pair (both swept; -period / -blockheight
// pin them). Cribs are not used (the per-block columnar scrambles plaintext positions).

#include "nicodemus_solver.h"
#include "nicodemus.h"
#include "engine.h"
#include "trans_common.h"
#include "scoring.h"

#define NICO_DEFAULT_MAXP   12   // default top of the period sweep when -maxcols is left at its global default
#define NICO_MIN_BLOCK_H     2   // bottom of the block-height sweep
#define NICO_DEFAULT_MAXH    8   // default top of the block-height sweep

typedef struct { int n; int minP, maxP, minH, maxH; int variant; } NicoScratch;

static int nico_variant_of(int cipher_type) {
    if (cipher_type == NICODEMUS_VARIANT)  return NICO_VARIANT;
    if (cipher_type == NICODEMUS_BEAUFORT) return NICO_BEAU;
    return NICO_VIG;
}

// Detranspose under `order`, derive the P per-column shifts by monogram fit, then inverse-
// substitute into plain[]. shifts_out[0..P-1] receives the recovered shifts.
static void nico_decrypt_core(const SolverCtx *ctx, int P, int H, const int *order,
                              int variant, int *shifts_out, int *plain) {
    int n = ctx->cipher_len;
    static _Thread_local int desub[MAX_CIPHER_LENGTH];
    nicodemus_detranspose(ctx->cipher, n, P, H, order, desub);

    for (int col = 0; col < P; col++) {
        int hist[ALPHABET_SIZE] = {0};
        for (int i = col; i < n; i += P) hist[desub[i]]++;
        double best = -1e300; int best_s = 0;
        for (int s = 0; s < ALPHABET_SIZE; s++) {
            double sc = 0.0;
            for (int c = 0; c < ALPHABET_SIZE; c++)
                if (hist[c]) sc += hist[c] * g_monograms[nicodemus_inv_sub(c, s, variant)];
            if (sc > best) { best = sc; best_s = s; }
        }
        shifts_out[col] = best_s;
    }
    nicodemus_inv_substitute(desub, n, P, shifts_out, variant, plain);
}

// ---- engine hooks ----------------------------------------------------------

static int nico_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const NicoScratch *a = (const NicoScratch *) ctx->model_scratch;
    int c = 0;
    for (int P = a->minP; P <= a->maxP; P++) {
        if (P < 2) continue;
        for (int H = a->minH; H <= a->maxH && c < cap; H++) {
            out[c].period = P; out[c].j = 0; out[c].k = 0;
            out[c].aux[0] = H; out[c].aux[1] = 0;
            c++;
        }
    }
    return c;
}

static void nico_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int P = cc->period;
    perm_seed(st->key, P);          // random column order
    st->aux[0] = P;
    st->aux[1] = cc->aux[0];        // block height H
    st->key_len = P;
}

static void nico_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                         SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    perm_move(st->key, st->aux[0]);
}

static void nico_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int P = src->aux[0];
    for (int i = 0; i < P; i++) dst->key[i] = src->key[i];
    dst->aux[0] = P;
    dst->aux[1] = src->aux[1];
    dst->key_len = P;
}

static void nico_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, int *out, double *score_adjust) {
    const NicoScratch *a = (const NicoScratch *) ctx->model_scratch;
    int P = st->aux[0], H = st->aux[1];
    int shifts[MAX_COLS];
    (void) cc;
    nico_decrypt_core(ctx, P, H, st->key, a->variant, shifts, out);
    *score_adjust = 0.0;
}

// Recovered substitution key as letters (one per column), plus the read order as digits.
static void nico_key_strings(const int shifts[], const int order[], int P,
                             char *keyword, char *orderstr) {
    for (int c = 0; c < P; c++) keyword[c] = index_to_char(shifts[c]);
    keyword[P] = '\0';
    int o = 0;
    for (int j = 0; j < P; j++) o += snprintf(orderstr + o, 8, "%d%s", order[j], j + 1 < P ? " " : "");
}

static void nico_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const NicoScratch *a = (const NicoScratch *) ctx->model_scratch;
    (void) cc; (void) decrypted;
    int P = st->aux[0], H = st->aux[1];
    int shifts[MAX_COLS], plain[MAX_CIPHER_LENGTH];
    nico_decrypt_core(ctx, P, H, st->key, a->variant, shifts, plain);
    char keyword[MAX_COLS + 1], orderstr[4 * MAX_COLS];
    nico_key_strings(shifts, st->key, P, keyword, orderstr);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  P=%d H=%d keyword=%s order=[%s] score=%.4f  [%.1fs, %d restarts]\n",
        P, H, keyword, orderstr, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void nico_report(const SolverCtx *ctx, const SolverConfig *cc,
                        const SolverState *st, double score, int *decrypted) {
    const NicoScratch *a = (const NicoScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len, P = st->aux[0], H = st->aux[1];
    int shifts[MAX_COLS], plain[MAX_CIPHER_LENGTH];
    (void) cc;
    nico_decrypt_core(ctx, P, H, st->key, a->variant, shifts, plain);

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char keyword[MAX_COLS + 1], orderstr[4 * MAX_COLS];
    nico_key_strings(shifts, st->key, P, keyword, orderstr);
    const char *vname = (a->variant == NICO_VARIANT) ? "variant"
                      : (a->variant == NICO_BEAU)    ? "beaufort" : "vigenere";

    printf("\nResult Score: %.2f | Words: %d | P=%d | H=%d | sub=%s | keyword=%s | order=[%s]\n",
        score, n_words_found, P, H, vname, keyword, orderstr);
    print_cipher(ctx->cipher, n, NULL);
    printf("\n");
    print_text(decrypted, n);
    printf("\n%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = P;
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, P=%d, H=%d, keyword=%s, order=[%s], ",
            score, n_words_found, cfg->cipher_type, P, H, keyword, orderstr);
    else
        printf(">>> %.2f, %d, P=%d, H=%d, keyword=%s, order=[%s], ",
            score, cfg->cipher_type, P, H, keyword, orderstr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel NICODEMUS_MODEL = {
    .name = "nicodemus", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = nico_enumerate, .key_len = NULL,
    .seed = nico_seed, .perturb = nico_perturb, .copy_state = nico_copy,
    .decrypt = nico_decrypt_hook, .report = nico_report,
    .report_verbose = nico_report_verbose,
};

// ---- entry point -----------------------------------------------------------

void solve_nicodemus(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;
    (void) crib_indices; (void) crib_positions; (void) n_cribs;   // cribs unused

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Nicodemus needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 8) {
        printf("\n\nERROR: ciphertext too short for a Nicodemus solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Nicodemus ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    NicoScratch a;
    a.n = cipher_len;
    a.variant = nico_variant_of(cfg->cipher_type);

    // Period sweep: -period pins; else -mincols..-maxcols, defaulting the top to
    // NICO_DEFAULT_MAXP when -maxcols is left at the global default (Nicodemus keywords
    // are short). Clamp to the array bound and the ciphertext length.
    if (cfg->period_present) { a.minP = a.maxP = cfg->period; }
    else {
        a.minP = cfg->min_cols < 2 ? 2 : cfg->min_cols;
        a.maxP = (cfg->max_cols == 30) ? NICO_DEFAULT_MAXP : cfg->max_cols;
    }
    if (a.maxP > MAX_COLS) a.maxP = MAX_COLS;
    if (a.maxP > cipher_len) a.maxP = cipher_len;
    if (a.maxP < a.minP) a.maxP = a.minP;

    // Block-height sweep: -blockheight pins; else NICO_MIN_BLOCK_H..-maxblockheight (default).
    if (cfg->block_height > 0) { a.minH = a.maxH = cfg->block_height; }
    else {
        a.minH = NICO_MIN_BLOCK_H;
        a.maxH = (cfg->max_block_height > 0) ? cfg->max_block_height : NICO_DEFAULT_MAXH;
    }
    if (a.maxH < a.minH) a.maxH = a.minH;

    if (cfg->verbose)
        printf("\nnicodemus (%s): %d letters, column-order anneal + derived shifts, "
               "P %d..%d, H %d..%d\n", (a.variant == NICO_VARIANT) ? "variant"
               : (a.variant == NICO_BEAU) ? "beaufort" : "vigenere",
               cipher_len, a.minP, a.maxP, a.minH, a.maxH);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, 0);
    ctx.model_scratch = &a;
    ctx.result = result;
    run_solver(&NICODEMUS_MODEL, &ctx);
}

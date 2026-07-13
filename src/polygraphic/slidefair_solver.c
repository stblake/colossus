#include "slidefair_solver.h"
#include "slidefair.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Slidefair solver (TYPES slidefair / slidefair-var / slidefair-beau)
// =====================================================================
//
// Slidefair (ACA "periodic digraphic Vigenere/Variant/Beaufort") enciphers the plaintext in
// consecutive DIGRAPHS over a two-row slide (slidefair.c). Digraph i (cipher letters 2i, 2i+1) is
// keyed by keyword letter key[i mod P], so writing P digraphs per row puts one key letter per
// column. The entire key is P key letters in 0..25 (a 26^P space), and the cipher is positional
// (decrypted[i] is plaintext[i], so cribs apply 1:1).
//
// The state is the cycleword lane: the P key letters (cycleword[0..P-1], each 0..25); one engine
// config per period P (IoC period estimation is useless through the digraphic pairing, so P is
// swept and the n-gram score picks it -- the rigid pairing makes a wrong P decrypt to gibberish).
// Because a digraph is enciphered ENTIRELY by its column key, every digraph in column c decrypts
// from key[c] alone -- the columns are independent given their keys. That lets the per-column
// monogram-fit key (decrypt the column's digraphs for each candidate key, pick the one whose
// decrypted letters best match English monograms -- the analog of derive_optimal_cycleword) WARM-
// START the seed; the n-gram (quadgram) score then drives the anneal, since cross-column digraphs
// only form at the true keys, and corrects any column the monogram fit mis-set. No score_adjust is
// needed (every cycleword yields a valid bijective decrypt; n-grams discriminate). Like the rest of
// the Vigenere family it rides the reward-only quadgram table (no -logprob needed).
//
// One primitive + one solver serve all three variants, branched on cfg->cipher_type. Vigenere and
// Variant are NOT separately identifiable (Variant key 26-k == Vigenere key k, and each per-column
// key is derived freely), so either solver cracks a shift-Slidefair; only Beaufort is distinct.

#define SLIDEFAIR_DEFAULT_MAXP 12   // default top of the period sweep when -maxcols is left at default
#define SLIDEFAIR_SEED_WARM    0.6  // P(a column is seeded at its monogram-best key vs random)

typedef struct { int n; int minP, maxP; int type; } SlidefairScratch;

// Per-column monogram-fit key: the key letter whose decrypted column letters best match English.
// Independent of the other columns (each digraph is enciphered solely by its column key).
static int slidefair_best_key(const int *cipher, int n, int P, int col, int type) {
    int ndg = n / 2;
    double best = -1e300;
    int best_k = 0;
    for (int kk = 0; kk < ALPHABET_SIZE; kk++) {
        double sc = 0.0;
        for (int i = col; i < ndg; i += P) {                  // every digraph in this column
            int p1, p2;
            slidefair_pair_dec(cipher[2 * i], cipher[2 * i + 1], kk, type, &p1, &p2);
            sc += g_monograms[p1] + g_monograms[p2];
        }
        if (sc > best) { best = sc; best_k = kk; }
    }
    return best_k;
}

// --- model hooks ------------------------------------------------------------------

static int slidefair_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const SlidefairScratch *a = (const SlidefairScratch *) ctx->model_scratch;
    int c = 0;
    for (int P = a->minP; P <= a->maxP && c < cap; P++) {     // one config per period
        out[c].period = P; out[c].j = 0; out[c].k = 0;
        out[c].aux[0] = 0; out[c].aux[1] = 0;
        c++;
    }
    return c;
}

static void slidefair_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const SlidefairScratch *a = (const SlidefairScratch *) ctx->model_scratch;
    int P = cc->period;
    for (int c = 0; c < P; c++) {
        // Warm-start most columns at the monogram-best key; randomise the rest so restarts
        // diverge and the anneal/n-gram pass can override a mis-set column.
        if (frand() < SLIDEFAIR_SEED_WARM)
            st->cycleword[c] = slidefair_best_key(ctx->cipher, ctx->cipher_len, P, c, a->type);
        else
            st->cycleword[c] = rand_int(0, ALPHABET_SIZE);
    }
    st->aux[0] = P;
    st->key_len = 0;                                 // cycleword lane only; key lane unused
}

static void slidefair_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int P = st->aux[0];
    int col = rand_int(0, P);
    int cur = st->cycleword[col], nv;
    do { nv = rand_int(0, ALPHABET_SIZE); } while (nv == cur);   // a different key for one column
    st->cycleword[col] = nv;
}

static void slidefair_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int P = src->aux[0];
    for (int i = 0; i < P; i++) dst->cycleword[i] = src->cycleword[i];
    dst->aux[0] = P;
    dst->key_len = 0;
}

static void slidefair_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                   SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const SlidefairScratch *a = (const SlidefairScratch *) ctx->model_scratch;
    slidefair_decrypt(out, ctx->cipher, ctx->cipher_len, st->cycleword, st->aux[0], a->type);
    *score_adjust = 0.0;
}

// --- reporting --------------------------------------------------------------------

// The recovered key as a keyword: each per-column key letter is the keyword letter directly.
// (For Vigenere/Variant the sign is not separately identifiable, so a Vigenere solve of a Variant
// cipher reports the complementary keyword -- both decrypt to the same plaintext.)
static void slidefair_key_string(const int keys[], int P, char out[]) {
    for (int c = 0; c < P; c++) out[c] = index_to_char(keys[c]);
    out[P] = '\0';
}

static const char *slidefair_variant_name(int type) {
    if (type == SLIDEFAIR_VAR)  return "Variant";
    if (type == SLIDEFAIR_BEAU) return "Beaufort";
    return "Vigenere";
}

static void slidefair_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    const SlidefairScratch *a = (const SlidefairScratch *) ctx->model_scratch;
    int P = st->aux[0];
    char kw[MAX_COLS + 1]; slidefair_key_string(st->cycleword, P, kw);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  %s P=%d keyword=%s score=%.4f  [%.1fs, %d restarts]\n",
        slidefair_variant_name(a->type), P, kw, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void slidefair_report(const SolverCtx *ctx, const SolverConfig *cc,
                             const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const SlidefairScratch *a = (const SlidefairScratch *) ctx->model_scratch;
    int n = ctx->cipher_len, P = st->aux[0];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char kw[MAX_COLS + 1]; slidefair_key_string(st->cycleword, P, kw);

    printf("\nResult Score: %.2f | Words: %d | %s | P=%d | keyword=%s\n",
        score, n_words_found, slidefair_variant_name(a->type), P, kw);
    print_cipher(ctx->cipher, n, NULL);
    printf("\n");
    print_text(decrypted, n);
    printf("\n");
    print_solution_check(decrypted, n);
    print_spaces_line(g_spaces_table, decrypted, n);
    printf("%s\n", ctx->cribtext);

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
        printf(">>> %.2f, %d, %d, P=%d, keyword=%s, ",
            score, n_words_found, cfg->cipher_type, P, kw);
    else
        printf(">>> %.2f, %d, P=%d, keyword=%s, ",
            score, cfg->cipher_type, P, kw);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel SLIDEFAIR_MODEL = {
    .name = "slidefair", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = slidefair_enumerate, .key_len = NULL,
    .seed = slidefair_seed, .perturb = slidefair_perturb, .copy_state = slidefair_copy,
    .decrypt = slidefair_decrypt_hook, .report = slidefair_report,
    .report_verbose = slidefair_report_verbose,
};

// =====================================================================
//  Entry point
// =====================================================================

void solve_slidefair(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Slidefair needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Slidefair solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Slidefair ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    // Period sweep: -period pins; else -mincols..-maxcols, defaulting the top to
    // SLIDEFAIR_DEFAULT_MAXP when -maxcols is left at its global default (keywords are short).
    // A period needs at least one full digraph per column, so clamp maxP to cipher_len/2.
    SlidefairScratch a;
    a.n = cipher_len;
    a.type = cfg->cipher_type;
    if (cfg->period_present) { a.minP = a.maxP = cfg->period; }
    else {
        a.minP = cfg->min_cols < 1 ? 1 : cfg->min_cols;
        if (a.minP < 1) a.minP = 1;
        a.maxP = (cfg->max_cols == 30) ? SLIDEFAIR_DEFAULT_MAXP : cfg->max_cols;
    }
    if (a.maxP > MAX_COLS) a.maxP = MAX_COLS;
    if (a.maxP > cipher_len / 2) a.maxP = cipher_len / 2;
    if (a.maxP < a.minP) a.maxP = a.minP;

    if (cfg->verbose)
        printf("\nslidefair: %d letters, %s, per-column key anneal, P %d..%d\n",
            cipher_len, slidefair_variant_name(a.type), a.minP, a.maxP);

    // Cribs are supported: the cipher is positional, so crib positions map straight onto the
    // decrypted plaintext (decrypted[i] == plaintext[i]).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &a;
    ctx.result = result;

    run_solver(&SLIDEFAIR_MODEL, &ctx);
}

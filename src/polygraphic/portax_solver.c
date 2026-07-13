#include "portax_solver.h"
#include "portax.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Portax solver (TYPE portax)
// =====================================================================
//
// Portax (ACA "periodic digraphic Porta") enciphers the plaintext in VERTICAL PAIRS over a Porta
// slide (portax.c). The plaintext is laid row-major at width P (= keyword length), rows are taken
// in pairs (2g, 2g+1), and the pair in column c is enciphered as a unit by the column key letter
// -- of which only the Porta SHIFT key/2 (0..12) matters. So the entire key is P shifts in 0..12
// (a 13^P space, the Porta cycleword domain), and the cipher is positional (decrypted[i] is
// plaintext[i], so cribs apply 1:1).
//
// The state is the cycleword lane: the P shifts (cycleword[0..P-1], each 0..12); one engine config
// per period P (IoC period estimation is useless through the digraphic pairing, so P is swept and
// the n-gram score picks it -- the rigid pairing makes a wrong P decrypt to gibberish). Because a
// vertical pair is enciphered ENTIRELY by its column key, every pair in column c decrypts from
// shift[c] alone -- the columns are independent given their shifts. That lets the per-column
// monogram-fit shift (decrypt the column's pairs for each candidate shift, pick the one whose
// decrypted letters best match English monograms -- the analog of derive_optimal_cycleword) WARM-
// START the seed; the n-gram (quadgram) score then drives the anneal, since cross-column digraphs
// only form at the true shifts, and corrects any column the monogram fit mis-set. No score_adjust
// is needed (every cycleword yields a valid bijective decrypt; n-grams discriminate). Like the
// other Porta-family ciphers it rides the reward-only quadgram table, but -logprob helps on short
// or hard ciphers.

#define PORTAX_DEFAULT_MAXP 12   // default top of the period sweep when -maxcols is left at default
#define PORTAX_SEED_WARM    0.6  // P(a column is seeded at its monogram-best shift vs random)

typedef struct { int n; int minP, maxP; } PortaxScratch;

// Per-column monogram-fit shift: the shift whose decrypted column letters best match English.
// Independent of the other columns (each pair is enciphered solely by its column key).
static int portax_best_shift(const int *cipher, int n, int P, int col) {
    int block = 2 * P;
    double best = -1e300;
    int best_s = 0;
    for (int s = 0; s < PORTAX_HALF; s++) {
        double sc = 0.0;
        for (int b = 0; b < n; b += block) {
            int it = b + col, ib = b + P + col;
            if (ib >= n) break;                       // no bottom partner in this (ragged) block
            int x, y;
            portax_pair(cipher[it], cipher[ib], s, &x, &y);
            sc += g_monograms[x] + g_monograms[y];
        }
        if (sc > best) { best = sc; best_s = s; }
    }
    return best_s;
}

// --- model hooks ------------------------------------------------------------------

static int portax_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const PortaxScratch *a = (const PortaxScratch *) ctx->model_scratch;
    int c = 0;
    for (int P = a->minP; P <= a->maxP && c < cap; P++) {     // one config per period
        out[c].period = P; out[c].j = 0; out[c].k = 0;
        out[c].aux[0] = 0; out[c].aux[1] = 0;
        c++;
    }
    return c;
}

static void portax_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    int P = cc->period;
    for (int c = 0; c < P; c++) {
        // Warm-start most columns at the monogram-best shift; randomise the rest so restarts
        // diverge and the anneal/n-gram pass can override a mis-set column.
        if (frand() < PORTAX_SEED_WARM)
            st->cycleword[c] = portax_best_shift(ctx->cipher, ctx->cipher_len, P, c);
        else
            st->cycleword[c] = rand_int(0, PORTAX_HALF);
    }
    st->aux[0] = P;
    st->key_len = 0;                                 // cycleword lane only; key lane unused
}

static void portax_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                           SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int P = st->aux[0];
    int col = rand_int(0, P);
    int cur = st->cycleword[col], nv;
    do { nv = rand_int(0, PORTAX_HALF); } while (nv == cur);   // a different shift for one column
    st->cycleword[col] = nv;
}

static void portax_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int P = src->aux[0];
    for (int i = 0; i < P; i++) dst->cycleword[i] = src->cycleword[i];
    dst->aux[0] = P;
    dst->key_len = 0;
}

static void portax_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    portax_apply(ctx->cipher, ctx->cipher_len, st->cycleword, st->aux[0], out);
    *score_adjust = 0.0;
}

// --- reporting --------------------------------------------------------------------

// The recovered key as a representative keyword: shift s -> the even letter of its Porta pair
// (A,C,E,...,Y). U<->V and the other pair members are not separately identifiable.
static void portax_key_string(const int shifts[], int P, char out[]) {
    for (int c = 0; c < P; c++) out[c] = index_to_char(2 * shifts[c]);
    out[P] = '\0';
}

static void portax_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    int P = st->aux[0];
    char kw[MAX_COLS + 1]; portax_key_string(st->cycleword, P, kw);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  P=%d keyword=%s score=%.4f  [%.1fs, %d restarts]\n",
        P, kw, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void portax_report(const SolverCtx *ctx, const SolverConfig *cc,
                          const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len, P = st->aux[0];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char kw[MAX_COLS + 1]; portax_key_string(st->cycleword, P, kw);

    printf("\nResult Score: %.2f | Words: %d | P=%d | keyword=%s\n",
        score, n_words_found, P, kw);
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

static const CipherModel PORTAX_MODEL = {
    .name = "portax", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = portax_enumerate, .key_len = NULL,
    .seed = portax_seed, .perturb = portax_perturb, .copy_state = portax_copy,
    .decrypt = portax_decrypt_hook, .report = portax_report,
    .report_verbose = portax_report_verbose,
};

// =====================================================================
//  Entry point
// =====================================================================

void solve_portax(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Portax needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Portax solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Portax ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    // Period sweep: -period pins; else -mincols..-maxcols, defaulting the top to
    // PORTAX_DEFAULT_MAXP when -maxcols is left at its global default (Portax keywords are
    // short). A period needs at least one full row-pair, so clamp maxP to cipher_len/2.
    PortaxScratch a;
    a.n = cipher_len;
    if (cfg->period_present) { a.minP = a.maxP = cfg->period; }
    else {
        a.minP = cfg->min_cols < 1 ? 1 : cfg->min_cols;
        if (a.minP < 1) a.minP = 1;
        a.maxP = (cfg->max_cols == 30) ? PORTAX_DEFAULT_MAXP : cfg->max_cols;
    }
    if (a.maxP > MAX_COLS) a.maxP = MAX_COLS;
    if (a.maxP > cipher_len / 2) a.maxP = cipher_len / 2;
    if (a.maxP < a.minP) a.maxP = a.minP;

    if (cfg->verbose)
        printf("\nportax: %d letters, per-column Porta-shift anneal, P %d..%d\n",
            cipher_len, a.minP, a.maxP);

    // Cribs are supported: the cipher is positional, so crib positions map straight onto the
    // decrypted plaintext (decrypted[i] == plaintext[i]).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &a;
    ctx.result = result;

    run_solver(&PORTAX_MODEL, &ctx);
}

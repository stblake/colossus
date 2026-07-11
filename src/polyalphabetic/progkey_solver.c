#include "progkey_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Progressive Key solver (TYPE progkey / progkey-var / progkey-beau)
// =====================================================================
//
// The Progressive Key cipher (progkey.c) is a periodic base cipher (Vig/Var/Beau) under a
// P-letter keyword whose key DRIFTS by a constant (g*prog) every group g = i/P. The whole key
// is therefore P per-column shifts (0..25, the keyword) plus a single progression index `prog`
// (0..25). IoC period estimation fails (within a column each group carries a different drifted
// shift, so columns are not monoalphabetic -- the same situation as autokey), so the PERIOD is
// brute-forced and the PROGRESSION enumerated: one engine config per (P, prog) pair, with prog
// carried in cc->aux[0].
//
// The key to efficiency is the same de-coupling the polyalphabetic -optimalcycle path uses:
// for a FIXED prog, DE-PROGRESSING the ciphertext (undoing only the drift pass, leaving the
// primary base cipher C1) makes every column an independent Caesar/Beaufort sample under its
// own keyword shift. So the per-column monogram-fit shift WARM-STARTS the seed (decrypt the
// de-progressed column for each candidate shift, keep the one whose decrypted letters best
// match English monograms -- the analog of derive_optimal_cycleword), and the n-gram (quadgram)
// score then drives the anneal AND, across all (P, prog) configs, selects the true period and
// progression (a wrong P or prog leaves columns drifted -> poor monogram fit -> gibberish ->
// low n-gram score). No score_adjust is needed (every cycleword yields a valid bijective
// decrypt; n-grams discriminate). Like the rest of the Vigenere family it rides the reward-only
// quadgram table (no -logprob needed); -logprob helps on short or hard ciphers.

#define PROGKEY_DEFAULT_MAXP 15   // default top of the period sweep when -maxcols is left at default
#define PROGKEY_SEED_WARM    0.7  // P(a column is seeded at its monogram-best shift vs random)

typedef struct { int n; int minP, maxP; int progLo, progHi; int base; } ProgkeyScratch;

// Per-column monogram-fit base shift for a candidate progression `prog`: de-progress this
// column's cells (undo the drift pass) and pick the keyword shift whose decrypted letters best
// match English monograms. Independent of the other columns once the cipher is de-progressed.
static int progkey_best_shift(const int *cipher, int n, int P, int prog, int col, int base) {
    double best = -1e300;
    int best_s = 0;
    for (int s = 0; s < ALPHABET_SIZE; s++) {
        double sc = 0.0;
        for (int i = col; i < n; i += P) {
            int g = i / P;
            int c1 = progkey_base_decrypt(cipher[i], (g * prog) % ALPHABET_SIZE, base); // undo drift
            int pt = progkey_base_decrypt(c1, s, base);                                 // undo keyword
            sc += g_monograms[pt];
        }
        if (sc > best) { best = sc; best_s = s; }
    }
    return best_s;
}

// --- model hooks ------------------------------------------------------------------

static int progkey_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const ProgkeyScratch *a = (const ProgkeyScratch *) ctx->model_scratch;
    int c = 0;
    for (int P = a->minP; P <= a->maxP; P++) {
        for (int prog = a->progLo; prog <= a->progHi && c < cap; prog++) {
            out[c].period = P; out[c].j = 0; out[c].k = 0;
            out[c].aux[0] = prog; out[c].aux[1] = 0;
            c++;
        }
    }
    return c;
}

static void progkey_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const ProgkeyScratch *a = (const ProgkeyScratch *) ctx->model_scratch;
    int P = cc->period, prog = cc->aux[0], base = a->base;
    for (int c = 0; c < P; c++) {
        // Warm-start most columns at the monogram-best shift; randomise the rest so restarts
        // diverge and the anneal/n-gram pass can override a mis-set column.
        if (ctx->cfg->optimal_cycleword && frand() < PROGKEY_SEED_WARM)
            st->cycleword[c] = progkey_best_shift(ctx->cipher, ctx->cipher_len, P, prog, c, base);
        else
            st->cycleword[c] = rand_int(0, ALPHABET_SIZE);
    }
    st->aux[0] = P;
    st->aux[1] = prog;
    st->key_len = 0;                                 // cycleword lane only; key lane unused
}

static void progkey_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int P = st->aux[0];
    int col = rand_int(0, P);
    int cur = st->cycleword[col], nv;
    do { nv = rand_int(0, ALPHABET_SIZE); } while (nv == cur);   // a different shift for one column
    st->cycleword[col] = nv;
}

static void progkey_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int P = src->aux[0];
    for (int i = 0; i < P; i++) dst->cycleword[i] = src->cycleword[i];
    dst->aux[0] = P;
    dst->aux[1] = src->aux[1];
    dst->key_len = 0;
}

static void progkey_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                 SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    int base = progkey_base(ctx->cfg->cipher_type);
    progkey_decrypt(out, ctx->cipher, ctx->cipher_len,
                    st->cycleword, st->aux[0], st->aux[1], base);
    *score_adjust = 0.0;
}

// --- reporting --------------------------------------------------------------------

// The recovered keyword: each per-column base shift 0..25 is a key letter directly (A..Z).
static void progkey_key_string(const int shifts[], int P, char out[]) {
    for (int c = 0; c < P; c++) out[c] = index_to_char(shifts[c]);
    out[P] = '\0';
}

static const char *progkey_base_name(int base) {
    return base == PROGKEY_BASE_VAR ? "variant"
         : base == PROGKEY_BASE_BEAU ? "beaufort" : "vigenere";
}

static void progkey_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    int P = st->aux[0];
    char kw[MAX_CYCLEWORD_LEN + 1]; progkey_key_string(st->cycleword, P, kw);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  P=%d prog=%d keyword=%s score=%.4f  [%.1fs, %d restarts]\n",
        P, st->aux[1], kw, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void progkey_report(const SolverCtx *ctx, const SolverConfig *cc,
                           const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len, P = st->aux[0], prog = st->aux[1];
    int base = progkey_base(cfg->cipher_type);

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char kw[MAX_CYCLEWORD_LEN + 1]; progkey_key_string(st->cycleword, P, kw);

    printf("\nResult Score: %.2f | Words: %d | base=%s | P=%d | prog=%d | keyword=%s\n",
        score, n_words_found, progkey_base_name(base), P, prog, kw);
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
        ctx->result->progression = prog;
        for (int i = 0; i < P; i++) ctx->result->cycleword[i] = st->cycleword[i];
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, base=%s, P=%d, prog=%d, keyword=%s, ",
            score, n_words_found, cfg->cipher_type, progkey_base_name(base), P, prog, kw);
    else
        printf(">>> %.2f, %d, base=%s, P=%d, prog=%d, keyword=%s, ",
            score, cfg->cipher_type, progkey_base_name(base), P, prog, kw);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel PROGKEY_MODEL = {
    .name = "progkey", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = progkey_enumerate, .key_len = NULL,
    .seed = progkey_seed, .perturb = progkey_perturb, .copy_state = progkey_copy,
    .decrypt = progkey_decrypt_hook, .report = progkey_report,
    .report_verbose = progkey_report_verbose,
};

// =====================================================================
//  Entry point
// =====================================================================

void solve_progkey(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Progressive Key needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Progressive Key solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Progressive Key ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    ProgkeyScratch a;
    a.n = cipher_len;
    a.base = progkey_base(cfg->cipher_type);

    // Period sweep: -period pins; else 1..(-maxcols), defaulting the top to PROGKEY_DEFAULT_MAXP
    // when -maxcols is left at its global default (progressive-key keywords are short). IoC
    // estimation is useless through the drift, so the period is brute-forced like autokey.
    if (cfg->cycleword_len_present) { a.minP = a.maxP = cfg->cycleword_len; }
    else if (cfg->period_present)   { a.minP = a.maxP = cfg->period; }
    else {
        a.minP = (cfg->min_cols >= 1) ? cfg->min_cols : 1;
        a.maxP = (cfg->max_cols == 30) ? PROGKEY_DEFAULT_MAXP : cfg->max_cols;
    }
    if (a.minP < 1) a.minP = 1;
    if (a.maxP > MAX_CYCLEWORD_LEN) a.maxP = MAX_CYCLEWORD_LEN;
    if (a.maxP > cipher_len) a.maxP = cipher_len;
    if (a.maxP < a.minP) a.maxP = a.minP;

    // Progression sweep: -progression pins; else the full 0..PROGKEY_MAX_PROG range (0 = plain
    // periodic cipher). A wrong prog leaves the columns drifted, so the n-gram score discards it.
    if (cfg->progression_present) {
        a.progLo = a.progHi = ((cfg->progression % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    } else {
        a.progLo = 0; a.progHi = PROGKEY_MAX_PROG;
    }

    if (cfg->verbose)
        printf("\nprogkey: %d letters, base=%s, per-column monogram-warm anneal, P %d..%d, prog %d..%d\n",
            cipher_len, progkey_base_name(a.base), a.minP, a.maxP, a.progLo, a.progHi);

    // Cribs are supported: the cipher is positional, so crib positions map straight onto the
    // decrypted plaintext (decrypted[i] == plaintext[i]).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &a;
    ctx.result = result;

    run_solver(&PROGKEY_MODEL, &ctx);
}

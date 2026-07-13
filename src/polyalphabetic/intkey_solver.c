#include "intkey_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Interrupted Key solver (TYPE intkey / intkey-var / intkey-beau)
// =====================================================================
//
// The Interrupted Key cipher (intkey.c) is a periodic base cipher (Vig/Var/Beau) under a P-letter
// keyword whose key index RESETS to the first key letter at break points. The whole key is P
// per-column shifts (the keyword) plus WHERE the breaks fall. IoC period estimation fails (the
// reset makes columns non-stationary, like autokey / progkey), so the PERIOD is swept, and the
// interruption STRATEGY (how breaks arise) is enumerated -- the n-gram score picks the winner.
//
// The efficiency lever is the same de-coupling the -optimalcycle path uses, plus a Gromark-style
// PRE-PASS. For the letter-interruptor strategies:
//   - CIPHERTEXT interruptor: the break positions are a function of the ciphertext ALONE (reset
//     after the interruptor CT letter), so the key index k[i] at every position is KNOWN once the
//     interruptor letter is fixed -- independent of the keyword. Grouping positions by k[i] makes
//     each column an independent Caesar/Beaufort sample, recovered by a monogram fit (the analog
//     of derive_optimal_cycleword). So for each (period, interruptor letter) we DERIVE the
//     monogram-best keyword deterministically and n-gram-score the decrypt -- a cheap pre-score.
//   - PLAINTEXT interruptor: the breaks depend on the produced plaintext (hence the keyword), so
//     there is no keyword-independent column map. An EM warm-start (fit keyword under the plain-
//     periodic assumption, causal-decrypt to get the actual columns, refit) gives a good seed,
//     then the same pre-score ranks it.
// The top-K (P, strategy, interruptor) by pre-score become engine configs, each annealing the
// keyword from a warm restart. A wrong period/interruptor/scheme leaves columns mis-fit ->
// gibberish -> low n-gram score, so the engine selects the true one. score_adjust stays 0 (every
// keyword yields a valid bijective decrypt; n-grams discriminate). Like the rest of the Vigenere
// family it rides the reward-only quadgram table (no -logprob needed; -logprob helps when short).
//
// SUPPLIED-BREAKS (random / word-division, -breaks): the mask is known, so it decouples exactly
// like the CT strategy (one config, monogram-fit warm start). JOINT (-intscheme joint): the mask
// is unknown, so the keyword AND an N-bit break-mask (in the key lane) are annealed together --
// genuinely blind but reliable only on longer text / with cribs (characterized, not asserted).

#define IK_DEFAULT_MAXP  15    // default top of the period sweep when -maxcols is left at default
#define IK_SEED_WARM     0.7   // P(a column is seeded at its monogram-best shift vs random)
#define IK_DEFAULT_TOPK  12    // pre-pass configs kept to anneal when -nprimers is not given
#define IK_MAX_CONFIGS   256   // hard cap on enumerated configs
#define IK_EM_ITERS      3     // plaintext-interruptor EM warm-start refinement passes
#define IK_JOINT_BREAK_PROB 0.25  // in JOINT, P(perturb the break-mask vs the keyword)

typedef struct {
    int n;
    int base;
    int n_configs;
    const int *strat;          // [n_configs]  IK_STRAT_*
    const int *interruptor;    // [n_configs]  0..25 (or -1)
    const int *period;         // [n_configs]
    const int *warm;           // [n_configs * MAX_CYCLEWORD_LEN]  warm keyword per config
    const int *user_break;     // [n]  supplied break mask (IK_STRAT_BREAKS), or NULL
} IntkeyScratch;

// Final config store (single-threaded), filled by solve_intkey before run_solver.
static int g_ik_strat[IK_MAX_CONFIGS];
static int g_ik_int[IK_MAX_CONFIGS];
static int g_ik_period[IK_MAX_CONFIGS];
static int g_ik_warm[IK_MAX_CONFIGS * MAX_CYCLEWORD_LEN];
static int g_ik_user_break[MAX_CIPHER_LENGTH];

// Scratch for the pre-pass / warm-start (single-threaded).
static int    g_ik_karr[MAX_CIPHER_LENGTH];
static int    g_ik_pt[MAX_CIPHER_LENGTH];
static int    g_ik_mask_scratch[MAX_CIPHER_LENGTH];
static double g_ik_colscore[MAX_CYCLEWORD_LEN * ALPHABET_SIZE];

// CT-strategy break-mask cache for the decrypt hot path (rebuilt when the config index changes).
static _Thread_local int g_ik_maskcache[MAX_CIPHER_LENGTH];
static _Thread_local int g_ik_maskcache_cfg = -1;

// ===================================================================
//  Column fit + strategy-specific key-index derivation
// ===================================================================

// Fill k_arr[i] = the key index used at position i for a given break mask (is_break[i]==1 forces
// k=0 at i; then increment mod P). Mirrors intkey_decrypt_mask's index walk.
static void ik_kindices_from_mask(const int *is_break, int n, int P, int *k_arr) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (is_break[i]) k = 0;
        k_arr[i] = k;
        if (++k >= P) k = 0;
    }
}

// Causal plaintext-interruptor walk: record the key index used at each position (resetting after
// the produced plaintext letter equals `interruptor`), and optionally the decrypted plaintext.
static void ik_ptint_kindices(const int *cipher, int n, const int *keyword, int P, int base,
                              int interruptor, int *k_arr, int *pt_out) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        int p = progkey_base_decrypt(cipher[i], keyword[k], base);
        if (pt_out) pt_out[i] = p;
        k_arr[i] = k;
        if (p == interruptor) k = 0;
        else if (++k >= P) k = 0;
    }
}

// Given the per-position key index k_arr, pick, for each column c in 0..P-1, the base shift whose
// decrypted letters best match English monograms (independent per column). The analog of
// derive_optimal_cycleword, but keyed on an arbitrary column map rather than a fixed stride.
static void ik_fit_keyword(const int *cipher, int n, int base, const int *k_arr, int P,
                           int *keyword_out) {
    for (int x = 0; x < P * ALPHABET_SIZE; x++) g_ik_colscore[x] = 0.0;
    for (int i = 0; i < n; i++) {
        double *row = &g_ik_colscore[k_arr[i] * ALPHABET_SIZE];
        for (int s = 0; s < ALPHABET_SIZE; s++)
            row[s] += g_monograms[progkey_base_decrypt(cipher[i], s, base)];
    }
    for (int c = 0; c < P; c++) {
        double best = -1e300; int bs = 0;
        double *row = &g_ik_colscore[c * ALPHABET_SIZE];
        for (int s = 0; s < ALPHABET_SIZE; s++) if (row[s] > best) { best = row[s]; bs = s; }
        keyword_out[c] = bs;
    }
}

// Pre-score one candidate (period, strategy, interruptor): derive the monogram-best keyword and
// return the n-gram score of the resulting decrypt, writing the warm keyword. Strategy CT/BREAKS
// use a known mask; PT uses the EM warm-start.
static double ik_prescore(const int *cipher, int n, int base, int P, int strat, int interruptor,
                          const int *user_break, const float *ngram, int ngram_size,
                          int *warm_out) {
    if (strat == IK_STRAT_PT) {
        for (int i = 0; i < n; i++) g_ik_karr[i] = i % P;          // plain-periodic assumption
        ik_fit_keyword(cipher, n, base, g_ik_karr, P, warm_out);
        for (int it = 0; it < IK_EM_ITERS; it++) {
            ik_ptint_kindices(cipher, n, warm_out, P, base, interruptor, g_ik_karr, NULL);
            ik_fit_keyword(cipher, n, base, g_ik_karr, P, warm_out);
        }
        intkey_decrypt_ptint(g_ik_pt, (int *) cipher, n, warm_out, P, base, interruptor);
        return ngram_score(g_ik_pt, n, (float *) ngram, ngram_size);
    }

    const int *mask;
    if (strat == IK_STRAT_BREAKS) {
        mask = user_break;
    } else { // IK_STRAT_CT
        intkey_build_mask_ct(g_ik_mask_scratch, cipher, n, interruptor);
        mask = g_ik_mask_scratch;
    }
    ik_kindices_from_mask(mask, n, P, g_ik_karr);
    ik_fit_keyword(cipher, n, base, g_ik_karr, P, warm_out);
    intkey_decrypt_mask(g_ik_pt, (int *) cipher, n, warm_out, P, base, mask);
    return ngram_score(g_ik_pt, n, (float *) ngram, ngram_size);
}

// ===================================================================
//  Top-K accumulator for the CT/PT pre-pass
// ===================================================================

#define IK_TOPK_MAX 64
typedef struct {
    int    K, count, worst;
    double score[IK_TOPK_MAX];
    int    strat[IK_TOPK_MAX];
    int    interruptor[IK_TOPK_MAX];
    int    period[IK_TOPK_MAX];
    int    warm[IK_TOPK_MAX][MAX_CYCLEWORD_LEN];
} IkTopK;

static void ik_topk_init(IkTopK *t, int K) {
    t->K = (K > IK_TOPK_MAX) ? IK_TOPK_MAX : (K < 1 ? 1 : K);
    t->count = 0; t->worst = 0;
}

static void ik_topk_refresh_worst(IkTopK *t) {
    int w = 0;
    for (int i = 1; i < t->count; i++) if (t->score[i] < t->score[w]) w = i;
    t->worst = w;
}

static void ik_topk_consider(IkTopK *t, double score, int strat, int interruptor, int P,
                             const int warm[]) {
    int slot;
    if (t->count < t->K) slot = t->count++;
    else if (score > t->score[t->worst]) slot = t->worst;
    else return;
    t->score[slot] = score;
    t->strat[slot] = strat;
    t->interruptor[slot] = interruptor;
    t->period[slot] = P;
    for (int i = 0; i < P; i++) t->warm[slot][i] = warm[i];
    if (t->count == t->K) ik_topk_refresh_worst(t);
}

// ===================================================================
//  Model hooks
// ===================================================================

static int ik_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const IntkeyScratch *a = (const IntkeyScratch *) ctx->model_scratch;
    int n = a->n_configs;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = a->period[i];
        out[i].j = i;                       // config index (warm-keyword / mask-cache lookup)
        out[i].k = 0;
        out[i].aux[0] = a->strat[i];
        out[i].aux[1] = a->interruptor[i];
    }
    return n;
}

static void ik_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const IntkeyScratch *a = (const IntkeyScratch *) ctx->model_scratch;
    int P = cc->period, strat = cc->aux[0], interruptor = cc->aux[1], idx = cc->j, n = ctx->cipher_len;
    const int *warm = &a->warm[idx * MAX_CYCLEWORD_LEN];
    for (int c = 0; c < P; c++) {
        if (ctx->cfg->optimal_cycleword && frand() < IK_SEED_WARM)
            st->cycleword[c] = warm[c];
        else
            st->cycleword[c] = rand_int(0, ALPHABET_SIZE);
    }
    st->aux[0] = P; st->aux[1] = strat; st->aux[2] = interruptor; st->aux[3] = idx; st->aux[4] = n;
    if (strat == IK_STRAT_JOINT) {
        for (int i = 0; i < n; i++) st->key[i] = 0;   // start with no interruptions (plain periodic)
        st->key_len = n;
    } else {
        st->key_len = 0;
    }
}

static void ik_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, bool *force_primary) {
    (void) cc; (void) force_primary;
    int P = st->aux[0], strat = st->aux[1], n = ctx->cipher_len;
    if (strat == IK_STRAT_JOINT && n > 1 && frand() < IK_JOINT_BREAK_PROB) {
        int pos = rand_int(1, n);                     // position 0 is always a group start
        st->key[pos] ^= 1;                            // toggle a break
    } else {
        int col = rand_int(0, P);
        int cur = st->cycleword[col], nv;
        do { nv = rand_int(0, ALPHABET_SIZE); } while (nv == cur);
        st->cycleword[col] = nv;
    }
}

static void ik_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int P = src->aux[0], strat = src->aux[1], n = src->aux[4];
    for (int i = 0; i < P; i++) dst->cycleword[i] = src->cycleword[i];
    for (int i = 0; i < 5; i++) dst->aux[i] = src->aux[i];
    if (strat == IK_STRAT_JOINT) {
        for (int i = 0; i < n; i++) dst->key[i] = src->key[i];
        dst->key_len = n;
    } else {
        dst->key_len = 0;
    }
}

static void ik_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const IntkeyScratch *a = (const IntkeyScratch *) ctx->model_scratch;
    int P = st->aux[0], strat = st->aux[1], interruptor = st->aux[2], idx = st->aux[3];
    int n = ctx->cipher_len, base = a->base;

    if (strat == IK_STRAT_PT) {
        intkey_decrypt_ptint(out, ctx->cipher, n, st->cycleword, P, base, interruptor);
    } else if (strat == IK_STRAT_JOINT) {
        intkey_decrypt_mask(out, ctx->cipher, n, st->cycleword, P, base, st->key);
    } else if (strat == IK_STRAT_BREAKS) {
        intkey_decrypt_mask(out, ctx->cipher, n, st->cycleword, P, base, a->user_break);
    } else { // IK_STRAT_CT: cache the keyword-independent mask per config
        if (g_ik_maskcache_cfg != idx) {
            intkey_build_mask_ct(g_ik_maskcache, ctx->cipher, n, interruptor);
            g_ik_maskcache_cfg = idx;
        }
        intkey_decrypt_mask(out, ctx->cipher, n, st->cycleword, P, base, g_ik_maskcache);
    }
    *score_adjust = 0.0;
}

// ===================================================================
//  Reporting
// ===================================================================

static void ik_key_string(const int shifts[], int P, char out[]) {
    for (int c = 0; c < P; c++) out[c] = index_to_char(shifts[c]);
    out[P] = '\0';
}

static const char *ik_base_name(int base) {
    return base == IK_BASE_VAR ? "variant" : base == IK_BASE_BEAU ? "beaufort" : "vigenere";
}

static const char *ik_scheme_name(int strat) {
    return strat == IK_STRAT_PT ? "pt-interruptor"
         : strat == IK_STRAT_BREAKS ? "supplied-breaks"
         : strat == IK_STRAT_JOINT ? "joint" : "ct-interruptor";
}

// Interruptor descriptor: a letter for CT/PT, "n/a" otherwise.
static void ik_int_string(int strat, int interruptor, char out[]) {
    if ((strat == IK_STRAT_CT || strat == IK_STRAT_PT) && interruptor >= 0)
        { out[0] = index_to_char(interruptor); out[1] = '\0'; }
    else { out[0] = 'n'; out[1] = '/'; out[2] = 'a'; out[3] = '\0'; }
}

static void ik_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    int P = st->aux[0], strat = st->aux[1], interruptor = st->aux[2];
    char kw[MAX_CYCLEWORD_LEN + 1]; ik_key_string(st->cycleword, P, kw);
    char is[8]; ik_int_string(strat, interruptor, is);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  P=%d scheme=%s interruptor=%s keyword=%s score=%.4f  [%.1fs, %d restarts]\n",
        P, ik_scheme_name(strat), is, kw, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void ik_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    (void) cc;
    const IntkeyScratch *a = (const IntkeyScratch *) ctx->model_scratch;
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len, P = st->aux[0], strat = st->aux[1], interruptor = st->aux[2];
    int base = a->base;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char kw[MAX_CYCLEWORD_LEN + 1]; ik_key_string(st->cycleword, P, kw);
    char is[8]; ik_int_string(strat, interruptor, is);

    printf("\nResult Score: %.2f | Words: %d | base=%s | P=%d | scheme=%s | interruptor=%s | keyword=%s\n",
        score, n_words_found, ik_base_name(base), P, ik_scheme_name(strat), is, kw);
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
        ctx->result->interruptor = ((strat == IK_STRAT_CT || strat == IK_STRAT_PT) ? interruptor : -1);
        ctx->result->intscheme = strat;
        for (int i = 0; i < P; i++) ctx->result->cycleword[i] = st->cycleword[i];
        vec_copy(decrypted, ctx->result->decrypted, n);
        ctx->result->decrypted_len = n;
    }

    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, base=%s, P=%d, scheme=%s, interruptor=%s, keyword=%s, ",
            score, n_words_found, cfg->cipher_type, ik_base_name(base), P, ik_scheme_name(strat), is, kw);
    else
        printf(">>> %.2f, %d, base=%s, P=%d, scheme=%s, interruptor=%s, keyword=%s, ",
            score, cfg->cipher_type, ik_base_name(base), P, ik_scheme_name(strat), is, kw);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel INTKEY_MODEL = {
    .name = "intkey", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = ik_enumerate, .key_len = NULL,
    .seed = ik_seed, .perturb = ik_perturb, .copy_state = ik_copy,
    .decrypt = ik_decrypt_hook, .report = ik_report,
    .report_verbose = ik_report_verbose,
};

// ===================================================================
//  Break-file parsing (supplied-breaks / random / word-division scheme)
// ===================================================================

// Read a whitespace/comma-separated list of 0-based group-start positions from `path` into a break
// mask (mask[pos]=1). Returns the largest group length seen (>=1) or -1 on error. The largest
// group length is the natural keyword length P (the ACA "entire keyword used at least once" rule).
static int ik_load_breaks(const char *path, int n, int *mask_out) {
    for (int i = 0; i < n; i++) mask_out[i] = 0;
    FILE *fp = fopen(path, "r");
    if (!fp) { printf("\n\nERROR: cannot open -breaks file '%s'.\n\n", path); return -1; }
    int pos, any = 0;
    while (fscanf(fp, " %d", &pos) == 1) {
        if (pos > 0 && pos < n) { mask_out[pos] = 1; any = 1; }
        else if (pos == 0) any = 1;   // position 0 is an implicit group start; accepted, no-op
    }
    fclose(fp);
    if (!any) { printf("\n\nERROR: -breaks file '%s' had no valid positions.\n\n", path); return -1; }

    int maxlen = 1, run = 0;          // largest run between consecutive breaks
    for (int i = 0; i < n; i++) {
        if (mask_out[i]) run = 0;
        run++;
        if (run > maxlen) maxlen = run;
    }
    return maxlen;
}

// ===================================================================
//  Entry point
// ===================================================================

void solve_intkey(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Interrupted Key needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for an Interrupted Key solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Interrupted Key ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    int base = intkey_base(cfg->cipher_type);
    int n = cipher_len;

    // Period sweep: -cyclewordlen / -period pin; else 1..(-maxcols), defaulting the top to
    // IK_DEFAULT_MAXP when -maxcols is at its global default (interrupted keywords are short).
    int minP, maxP;
    if (cfg->cycleword_len_present) { minP = maxP = cfg->cycleword_len; }
    else if (cfg->period_present)   { minP = maxP = cfg->period; }
    else {
        minP = (cfg->min_cols >= 1) ? cfg->min_cols : 1;
        maxP = (cfg->max_cols == 30) ? IK_DEFAULT_MAXP : cfg->max_cols;
    }
    if (minP < 2) minP = 2;                          // P==1 is a trivial Caesar; skip in the blind sweep
    if (maxP > MAX_CYCLEWORD_LEN) maxP = MAX_CYCLEWORD_LEN;
    if (maxP > n) maxP = n;
    if (maxP < minP) maxP = minP;

    // Interruptor letter range: -interruptor pins one, else 0..25.
    int intLo = 0, intHi = ALPHABET_SIZE - 1;
    if (cfg->interruptor_present) {
        int L = ((cfg->interruptor % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
        intLo = intHi = L;
    }

    // Which strategies to attack.
    bool do_ct = false, do_pt = false, do_breaks = false, do_joint = false;
    if (cfg->intscheme_present) {
        switch (cfg->intscheme) {
            case IK_STRAT_CT:     do_ct = true; break;
            case IK_STRAT_PT:     do_pt = true; break;
            case IK_STRAT_BREAKS: do_breaks = true; break;
            case IK_STRAT_JOINT:  do_joint = true; break;
            default:              do_ct = do_pt = true; break;
        }
    } else if (cfg->breaks_present) {
        do_breaks = true;                            // -breaks implies the supplied-breaks scheme
    } else {
        do_ct = do_pt = true;                        // blind default: enumerate CT + PT
    }

    // Supplied breaks (random / word-division): load the mask; P is the max group length (or -period).
    int breaksP = 0;
    if (do_breaks) {
        if (!cfg->breaks_present) {
            printf("\n\nERROR: -intscheme breaks requires -breaks <file>.\n\n");
            return;
        }
        breaksP = ik_load_breaks(cfg->breaks_file, n, g_ik_user_break);
        if (breaksP < 1) return;
        if (cfg->period_present || cfg->cycleword_len_present) breaksP = maxP;   // honour an explicit pin
    }

    int K = (cfg->n_primers > 0) ? cfg->n_primers : IK_DEFAULT_TOPK;

    // --- Build the config list ---
    int nc = 0;

    if (do_ct || do_pt) {
        // Gromark-style pre-pass: score every (period, strategy, interruptor), keep the top-K.
        static IkTopK top;
        ik_topk_init(&top, K);
        static int warm[MAX_CYCLEWORD_LEN];
        for (int P = minP; P <= maxP; P++) {
            for (int L = intLo; L <= intHi; L++) {
                if (do_ct) {
                    double sc = ik_prescore(cipher_indices, n, base, P, IK_STRAT_CT, L,
                        NULL, shared->ngram_data, cfg->ngram_size, warm);
                    ik_topk_consider(&top, sc, IK_STRAT_CT, L, P, warm);
                }
                if (do_pt) {
                    double sc = ik_prescore(cipher_indices, n, base, P, IK_STRAT_PT, L,
                        NULL, shared->ngram_data, cfg->ngram_size, warm);
                    ik_topk_consider(&top, sc, IK_STRAT_PT, L, P, warm);
                }
            }
        }
        // Sort kept configs by pre-score descending (K small; simple selection).
        int order[IK_TOPK_MAX];
        for (int i = 0; i < top.count; i++) order[i] = i;
        for (int aidx = 0; aidx < top.count; aidx++) {
            int b = aidx;
            for (int c = aidx + 1; c < top.count; c++)
                if (top.score[order[c]] > top.score[order[b]]) b = c;
            int t = order[aidx]; order[aidx] = order[b]; order[b] = t;
        }
        for (int i = 0; i < top.count && nc < IK_MAX_CONFIGS; i++) {
            int s = order[i];
            g_ik_strat[nc] = top.strat[s];
            g_ik_int[nc]   = top.interruptor[s];
            g_ik_period[nc] = top.period[s];
            for (int c = 0; c < top.period[s]; c++)
                g_ik_warm[nc * MAX_CYCLEWORD_LEN + c] = top.warm[s][c];
            nc++;
        }
        if (cfg->verbose && top.count)
            printf("\nintkey: pre-pass kept top %d (best n-gram fit %.4f)\n",
                top.count, top.score[order[0]]);
    }

    if (do_breaks && nc < IK_MAX_CONFIGS) {
        int *warm = &g_ik_warm[nc * MAX_CYCLEWORD_LEN];
        ik_kindices_from_mask(g_ik_user_break, n, breaksP, g_ik_karr);
        ik_fit_keyword(cipher_indices, n, base, g_ik_karr, breaksP, warm);
        g_ik_strat[nc] = IK_STRAT_BREAKS; g_ik_int[nc] = -1; g_ik_period[nc] = breaksP;
        nc++;
    }

    if (do_joint) {
        for (int P = minP; P <= maxP && nc < IK_MAX_CONFIGS; P++) {
            int *warm = &g_ik_warm[nc * MAX_CYCLEWORD_LEN];
            for (int i = 0; i < n; i++) g_ik_karr[i] = i % P;      // plain-periodic warm keyword
            ik_fit_keyword(cipher_indices, n, base, g_ik_karr, P, warm);
            g_ik_strat[nc] = IK_STRAT_JOINT; g_ik_int[nc] = -1; g_ik_period[nc] = P;
            nc++;
        }
    }

    if (nc < 1) {
        printf("\n\nERROR: Interrupted Key produced no candidate configurations.\n\n");
        return;
    }

    if (cfg->verbose)
        printf("intkey: %d letters, base=%s, P %d..%d, %d configs (ct=%d pt=%d breaks=%d joint=%d)\n",
            n, ik_base_name(base), minP, maxP, nc, do_ct, do_pt, do_breaks, do_joint);

    IntkeyScratch a;
    a.n = n; a.base = base; a.n_configs = nc;
    a.strat = g_ik_strat; a.interruptor = g_ik_int; a.period = g_ik_period; a.warm = g_ik_warm;
    a.user_break = do_breaks ? g_ik_user_break : NULL;
    g_ik_maskcache_cfg = -1;                         // invalidate the CT mask cache

    // Cribs are supported: the cipher is positional, so crib positions map straight onto the
    // decrypted plaintext (decrypted[i] == plaintext[i]).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &a;
    ctx.result = result;

    run_solver(&INTKEY_MODEL, &ctx);
}

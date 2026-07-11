#include "gromark_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Gromark / Periodic Gromark solver (TYPE gromark / gromark-periodic)
// =====================================================================
//
// Gromark composes a keyed 26-letter substitution sigma with a chain-addition running key d[]
// from a P-digit primer (C[i] = sigma[(p[i] + d[i]) mod 26]; Periodic Gromark adds a per-group
// offset). Breaking it is a COUPLED problem -- recover sigma AND the primer -- but unlike ADFGVX /
// Nihilist-Sub there is NO square-independent reward that isolates the primer: the additive shift
// sits INSIDE the permutation sigma, so it cannot be stripped in cipher space, and chain addition
// amplifies a single wrong primer digit into a fully-wrong running key, so naive joint annealing
// of (sigma, primer) never accepts a primer move. The two types are attacked by DIFFERENT models:
//
// BASIC Gromark (GROMARK_MODEL) -- a PRIMER PRE-PASS (the analog of bifid_estimate_periods). The
//   primer space is small (10^P); for a FIXED primer the chain d[] is known, so recovering sigma
//   collapses to a monoalphabetic-substitution-with-known-per-position-shift: group cipher
//   positions by symbol c (all share sigma_inv[c]=u, plaintext (u - d[i]) mod 26) and pick the
//   assignment c->u maximizing an English monogram fit (a 26x26 Hungarian) -> a provisional sigma.
//   The primer is RANKED by the n-gram score of that provisional decrypt; the top-K primers become
//   engine configs, each annealing sigma from a random restart (warm-started from the provisional
//   sigma). score_adjust stays 0 (the primer is pinned per config). State: key[0..25] = sigma;
//   aux[0]=P, aux[2]=config index (-> primer in model_scratch).
//
// PERIODIC Gromark (GROMARK_KW_MODEL) -- a KEYWORD anneal. Its whole key is ONE keyword of P
//   distinct letters (~28 bits) that derives sigma, the primer (the keyword letters' ranks) AND
//   the per-group offsets (their positions in sigma) together, so they cannot be searched as
//   independent unknowns. Annealing the keyword directly (rebuilding everything each decrypt) is a
//   tiny, navigable keyspace. See the GROMARK_KW_MODEL section below.

#define GM_SIGMA      ALPHABET_SIZE      // 26-letter keyed alphabet
#define GM_TOPK_MAX   192                // cap on enumerated primer configs (each annealed)
#define GM_DIGITS     10                 // running-key digits 0..9

typedef struct {
    int  n;                              // cipher length (== plaintext length)
    int  n_primers;                      // number of configs (top-K primers)
    const int *periods;                  // [n_primers]  P per config
    const int *primers;                  // [n_primers * GROMARK_MAX_PRIMER]
    const int *warm;                     // [n_primers * GM_SIGMA] provisional sigma per config
} GromarkScratch;

// Per-config primer/sigma store (single-threaded), filled by the pre-pass.
static int g_gm_periods[GM_TOPK_MAX];
static int g_gm_primers[GM_TOPK_MAX * GROMARK_MAX_PRIMER];
static int g_gm_warm[GM_TOPK_MAX * GM_SIGMA];

// Running-key cache for the current config (rebuilt when the config index changes).
static _Thread_local int g_gm_chain[MAX_CIPHER_LENGTH];
static _Thread_local int g_gm_chain_cfg = -1;

// ===================================================================
//  Pre-pass: rank primers by an assignment-based frequency attack
// ===================================================================

// English monogram log-weights for the assignment cost (natural log of g_monograms).
static double g_logmono[GM_SIGMA];
static void gm_init_logmono(void) {
    for (int i = 0; i < GM_SIGMA; i++) {
        double m = g_monograms[i];
        g_logmono[i] = log(m > 1e-9 ? m : 1e-9);
    }
}

// Max-weight perfect assignment on an N x N matrix (N <= GM_SIGMA), gain[r*N + c].
// Fills col_of_row[r] = matched column. O(N^3) Hungarian (min-cost on negated gain,
// the classic potentials + augmenting-path formulation, 1-indexed internally).
static void gm_hungarian_max(const double *gain, int N, int col_of_row[]) {
    const double INF = 1e18;
    double u[GM_SIGMA + 1], v[GM_SIGMA + 1], minv[GM_SIGMA + 1];
    int p[GM_SIGMA + 1], way[GM_SIGMA + 1];
    bool used[GM_SIGMA + 1];
    for (int i = 0; i <= N; i++) { u[i] = v[i] = 0.0; p[i] = 0; way[i] = 0; }

    for (int i = 1; i <= N; i++) {
        p[0] = i;
        int j0 = 0;
        for (int j = 0; j <= N; j++) { minv[j] = INF; used[j] = false; }
        do {
            used[j0] = true;
            int i0 = p[j0], j1 = -1;
            double delta = INF;
            for (int j = 1; j <= N; j++) if (!used[j]) {
                double cur = (-gain[(i0 - 1) * N + (j - 1)]) - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j = 0; j <= N; j++) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
    }
    for (int j = 1; j <= N; j++) col_of_row[p[j] - 1] = j - 1;
}

// Top-K accumulator (best n-gram scores seen, with their primer/period/provisional sigma).
typedef struct {
    int    K, count, worst;       // worst = index of the current minimum-score entry
    double score[GM_TOPK_MAX];
    int    period[GM_TOPK_MAX];
    int    primer[GM_TOPK_MAX][GROMARK_MAX_PRIMER];
    int    sigma[GM_TOPK_MAX][GM_SIGMA];
} GmTopK;

static void gm_topk_init(GmTopK *t, int K) {
    t->K = (K > GM_TOPK_MAX) ? GM_TOPK_MAX : K;
    t->count = 0; t->worst = 0;
}

static void gm_topk_refresh_worst(GmTopK *t) {
    int w = 0;
    for (int i = 1; i < t->count; i++) if (t->score[i] < t->score[w]) w = i;
    t->worst = w;
}

static void gm_topk_consider(GmTopK *t, double score, int P,
                             const int primer[], const int sigma[]) {
    int slot = -1;
    if (t->count < t->K) {
        slot = t->count++;
    } else if (score > t->score[t->worst]) {
        slot = t->worst;
    } else {
        return;
    }
    t->score[slot] = score;
    t->period[slot] = P;
    for (int i = 0; i < P; i++) t->primer[slot][i] = primer[i];
    for (int i = 0; i < GM_SIGMA; i++) t->sigma[slot][i] = sigma[i];
    if (t->count == t->K) gm_topk_refresh_worst(t);
}

// Evaluate one candidate basic-Gromark primer: with the running key d[] known, recovering the
// keyed alphabet sigma is a monoalphabetic-substitution-with-known-per-position-shift problem --
// group positions by cipher symbol c (all share sigma_inv[c]=u, plaintext (u - d[i]) mod 26) and
// pick the assignment c->u maximizing an English monogram fit (a 26x26 Hungarian). The primer is
// then ranked by the n-gram score of the resulting provisional decrypt. Returns the score and
// writes the provisional sigma. (Periodic Gromark does NOT use this -- it anneals the keyword.)
static double gm_eval_primer(const int cipher[], int n, int P, const int primer[],
                             const float *ngram, int ngram_size, int sigma_out[]) {
    static _Thread_local int d[MAX_CIPHER_LENGTH];
    gromark_chain_key(primer, P, n, d);

    // gain[c][u] = sum over positions with cipher==c of logmono[(u - d[i]) mod 26].
    static double gain[GM_SIGMA * GM_SIGMA];
    for (int x = 0; x < GM_SIGMA * GM_SIGMA; x++) gain[x] = 0.0;
    for (int i = 0; i < n; i++) {
        const double *lm = g_logmono;
        double *row = &gain[cipher[i] * GM_SIGMA];
        for (int u = 0; u < GM_SIGMA; u++) row[u] += lm[(u - d[i] + GM_SIGMA) % GM_SIGMA];
    }
    int col_of_row[GM_SIGMA], sinv[GM_SIGMA];
    gm_hungarian_max(gain, GM_SIGMA, col_of_row);
    for (int c = 0; c < GM_SIGMA; c++) { sinv[c] = col_of_row[c]; sigma_out[col_of_row[c]] = c; }

    static int pt[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) pt[i] = (sinv[cipher[i]] - d[i] + GM_SIGMA) % GM_SIGMA;
    return ngram_score(pt, n, (float *) ngram, ngram_size);
}

int gromark_rank_primers(const int cipher[], int n, int variant,
                         int fixed_period, int minP, int maxP,
                         const float *ngram, int ngram_size, int K,
                         int out_periods[], int out_primers[], int out_warm[],
                         bool verbose) {
    (void) variant; (void) minP; (void) maxP;          // basic Gromark only (periodic anneals the keyword)
    gm_init_logmono();
    static GmTopK top;
    gm_topk_init(&top, K);

    int P = fixed_period;
    long space = 1; for (int i = 0; i < P; i++) space *= 10;   // 10^P
    int primer[GROMARK_MAX_PRIMER], sigma[GM_SIGMA];
    for (long val = 0; val < space; val++) {
        long x = val;
        for (int i = P - 1; i >= 0; i--) { primer[i] = (int) (x % 10); x /= 10; }
        double sc = gm_eval_primer(cipher, n, P, primer, ngram, ngram_size, sigma);
        gm_topk_consider(&top, sc, P, primer, sigma);
    }

    // Sort the kept entries by score descending (simple selection; K is small).
    int order[GM_TOPK_MAX];
    for (int i = 0; i < top.count; i++) order[i] = i;
    for (int a = 0; a < top.count; a++) {
        int b = a;
        for (int c = a + 1; c < top.count; c++)
            if (top.score[order[c]] > top.score[order[b]]) b = c;
        int t = order[a]; order[a] = order[b]; order[b] = t;
    }

    for (int k = 0; k < top.count; k++) {
        int s = order[k];
        out_periods[k] = top.period[s];
        for (int i = 0; i < GROMARK_MAX_PRIMER; i++)
            out_primers[k * GROMARK_MAX_PRIMER + i] = (i < top.period[s]) ? top.primer[s][i] : 0;
        for (int i = 0; i < GM_SIGMA; i++) out_warm[k * GM_SIGMA + i] = top.sigma[s][i];
    }
    if (verbose) {
        printf("\nGromark primer pre-pass: kept top %d of the primer space "
               "(best n-gram fit %.4f, period(s)", top.count, top.count ? top.score[order[0]] : 0.0);
        int last = -1;
        for (int k = 0; k < top.count; k++) if (out_periods[k] != last) { printf(" %d", out_periods[k]); last = out_periods[k]; }
        printf(")\n");
    }
    return top.count;
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int gm_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const GromarkScratch *a = (const GromarkScratch *) ctx->model_scratch;
    int n = a->n_primers;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = a->periods[i];
        out[i].j = i; out[i].k = 0;
        out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

static void gm_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const GromarkScratch *a = (const GromarkScratch *) ctx->model_scratch;
    int P = cc->period;
    const int *w = &a->warm[cc->j * GM_SIGMA];
    // Warm-start ~40% of restarts from the pre-pass's provisional sigma; else a random alphabet.
    bool warm = frand() < 0.40;
    if (warm) {
        for (int i = 0; i < GM_SIGMA; i++) st->key[i] = w[i];
    } else {
        for (int i = 0; i < GM_SIGMA; i++) st->key[i] = i;
        shuffle(st->key, GM_SIGMA);
    }
    st->aux[0] = P;
    st->aux[2] = cc->j;
    st->key_len = GM_SIGMA;
}

// Simple-substitution swap on sigma (sigma has no Polybius geometry, so plain cell swaps --
// NOT bifid's row/column grid moves, which would be a silent correctness bug).
static void gm_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                       SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int a = rand_int(0, GM_SIGMA), b = rand_int(0, GM_SIGMA);
    int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
}

static void gm_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < GM_SIGMA; i++) dst->key[i] = src->key[i];
    for (int i = 0; i < 3; i++) dst->aux[i] = src->aux[i];
    dst->key_len = src->key_len;
}

static void gm_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, int *out, double *score_adjust) {
    const GromarkScratch *a = (const GromarkScratch *) ctx->model_scratch;
    int P = st->aux[0], cfgidx = st->aux[2], n = ctx->cipher_len;

    // The running key depends only on the (pinned) primer -> cache it per config.
    if (g_gm_chain_cfg != cfgidx) {
        gromark_chain_key(&a->primers[cfgidx * GROMARK_MAX_PRIMER], P, n, g_gm_chain);
        g_gm_chain_cfg = cfgidx;
    }
    int sinv[GM_SIGMA];
    for (int i = 0; i < GM_SIGMA; i++) sinv[st->key[i]] = i;
    gromark_decrypt_core(ctx->cipher, n, sinv, g_gm_chain, NULL, P, out);   // basic: no offsets
    *score_adjust = 0.0;
    (void) cc;
}

// ===================================================================
//  Reporting
// ===================================================================

static void gm_alpha_string(const int sigma[], char out[]) {
    for (int i = 0; i < GM_SIGMA; i++) out[i] = index_to_char(sigma[i]);
    out[GM_SIGMA] = '\0';
}

static void gm_primer_string(const int primer[], int P, char out[]) {
    int w = 0;
    for (int i = 0; i < P; i++) w += sprintf(out + w, "%d", primer[i] % 10);
}

static void gm_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    const GromarkScratch *a = (const GromarkScratch *) ctx->model_scratch;
    int P = st->aux[0], cfgidx = st->aux[2];
    char alpha[GM_SIGMA + 1], primer[GROMARK_MAX_PRIMER + 1];
    gm_alpha_string(st->key, alpha);
    gm_primer_string(&a->primers[cfgidx * GROMARK_MAX_PRIMER], P, primer);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  period %d, primer %s, score=%.4f  [%.1fs, %d restarts]\n    alphabet=%s\n",
        P, primer, score, elapsed, stats->n_restarts, alpha);
    (void) cc; (void) decrypted;
    fflush(stdout);
}

static void gm_report(const SolverCtx *ctx, const SolverConfig *cc,
                      const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    const GromarkScratch *a = (const GromarkScratch *) ctx->model_scratch;
    int n = ctx->cipher_len, P = st->aux[0], cfgidx = st->aux[2];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char alpha[GM_SIGMA + 1], primer[GROMARK_MAX_PRIMER + 1];
    gm_alpha_string(st->key, alpha);
    gm_primer_string(&a->primers[cfgidx * GROMARK_MAX_PRIMER], P, primer);

    printf("\nResult Score: %.2f | Words: %d | period=%d | primer=%s | alphabet=%s\n",
        score, n_words_found, P, primer, alpha);

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

    // One-liner summary: >>> score, [words,] type, period=, primer=, alphabet=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, period=%d, primer=%s, alphabet=%s, ",
            score, n_words_found, cfg->cipher_type, P, primer, alpha);
    else
        printf(">>> %.2f, %d, period=%d, primer=%s, alphabet=%s, ",
            score, cfg->cipher_type, P, primer, alpha);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
    (void) cc;
}

static const CipherModel GROMARK_MODEL = {
    .name = "gromark", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = gm_enumerate, .key_len = NULL,
    .seed = gm_seed, .perturb = gm_perturb, .copy_state = gm_copy,
    .decrypt = gm_decrypt_hook, .report = gm_report,
    .report_verbose = gm_report_verbose,
};

// ===================================================================
//  Periodic Gromark: a KEYWORD anneal (the whole key is one keyword)
// ===================================================================
//
// In Periodic Gromark the ENTIRE key is a single keyword of P distinct letters: it derives the
// K2M mixed alphabet sigma, the primer (the keyword letters' alphabetical ranks), AND the per-
// group offsets (their positions in sigma) all together. Treating sigma / primer / offsets as
// independent unknowns (as the basic-Gromark pre-pass does) blows a ~28-bit key up into a 26! x
// P! x 26^P coupled space the primer cannot be ranked in. So Periodic Gromark is solved by
// annealing the KEYWORD directly -- a tiny, navigable keyspace -- rebuilding sigma/primer/offsets
// from it each decrypt. One engine config per candidate period P (swept like a fractionation
// period); no primer pre-pass. State: st->key[0..P-1] = the keyword letter indices (distinct).

typedef struct { int n; int minP, maxP; } GromarkKwScratch;

static int gmkw_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const GromarkKwScratch *a = (const GromarkKwScratch *) ctx->model_scratch;
    int c = 0;
    for (int P = a->minP; P <= a->maxP && c < cap; P++) {
        if (P < 2 || P > 9) continue;
        out[c].period = P; out[c].j = 0; out[c].k = 0; out[c].aux[0] = 0; out[c].aux[1] = 0;
        c++;
    }
    return c;
}

static void gmkw_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    int P = cc->period;
    int perm[GM_SIGMA];
    for (int i = 0; i < GM_SIGMA; i++) perm[i] = i;
    shuffle(perm, GM_SIGMA);
    for (int i = 0; i < P; i++) st->key[i] = perm[i];      // P distinct random letters
    st->aux[0] = P;
    st->key_len = P;
}

static void gmkw_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                         SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int P = st->aux[0];
    if (P >= 2 && frand() < 0.40) {                        // swap two keyword positions
        int a = rand_int(0, P), b = rand_int(0, P);
        int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
    } else {                                               // replace one letter with an unused one
        bool used[GM_SIGMA] = {false};
        for (int i = 0; i < P; i++) used[st->key[i]] = true;
        int pool[GM_SIGMA], m = 0;
        for (int i = 0; i < GM_SIGMA; i++) if (!used[i]) pool[m++] = i;
        if (m > 0) st->key[rand_int(0, P)] = pool[rand_int(0, m)];
    }
}

static void gmkw_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    int P = src->aux[0];
    for (int i = 0; i < P; i++) dst->key[i] = src->key[i];
    dst->aux[0] = P;
    dst->key_len = P;
}

static void gmkw_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                              SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    int P = st->aux[0], n = ctx->cipher_len;
    int sigma[GM_SIGMA], sinv[GM_SIGMA], primer[GROMARK_MAX_PRIMER], offsets[GROMARK_MAX_PRIMER];
    static _Thread_local int d[MAX_CIPHER_LENGTH];
    gromark_build_from_keyword_idx(st->key, P, sigma, primer, offsets);
    for (int i = 0; i < GM_SIGMA; i++) sinv[sigma[i]] = i;
    gromark_chain_key(primer, P, n, d);
    gromark_decrypt_core(ctx->cipher, n, sinv, d, offsets, P, out);
    *score_adjust = 0.0;
}

static void gmkw_kw_string(const int kw[], int P, char out[]) {
    for (int i = 0; i < P; i++) out[i] = index_to_char(kw[i]);
    out[P] = '\0';
}

static void gmkw_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) ctx; (void) cc; (void) decrypted;
    int P = st->aux[0];
    char kw[GROMARK_MAX_PRIMER + 1]; gmkw_kw_string(st->key, P, kw);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  period %d, keyword %s, score=%.4f  [%.1fs, %d restarts]\n",
        P, kw, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void gmkw_report(const SolverCtx *ctx, const SolverConfig *cc,
                        const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    int n = ctx->cipher_len, P = st->aux[0];
    int sigma[GM_SIGMA], primer[GROMARK_MAX_PRIMER], offsets[GROMARK_MAX_PRIMER];
    gromark_build_from_keyword_idx(st->key, P, sigma, primer, offsets);

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char kw[GROMARK_MAX_PRIMER + 1]; gmkw_kw_string(st->key, P, kw);
    char alpha[GM_SIGMA + 1]; for (int i = 0; i < GM_SIGMA; i++) alpha[i] = index_to_char(sigma[i]);
    alpha[GM_SIGMA] = '\0';
    char primstr[GROMARK_MAX_PRIMER + 1]; for (int i = 0; i < P; i++) primstr[i] = '0' + primer[i] % 10;
    primstr[P] = '\0';

    printf("\nResult Score: %.2f | Words: %d | period=%d | keyword=%s | primer=%s | alphabet=%s\n",
        score, n_words_found, P, kw, primstr, alpha);
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
        printf(">>> %.2f, %d, %d, period=%d, keyword=%s, primer=%s, alphabet=%s, ",
            score, n_words_found, cfg->cipher_type, P, kw, primstr, alpha);
    else
        printf(">>> %.2f, %d, period=%d, keyword=%s, primer=%s, alphabet=%s, ",
            score, cfg->cipher_type, P, kw, primstr, alpha);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, n, NULL);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
    (void) cc;
}

static const CipherModel GROMARK_KW_MODEL = {
    .name = "gromark-periodic", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = gmkw_enumerate, .key_len = NULL,
    .seed = gmkw_seed, .perturb = gmkw_perturb, .copy_state = gmkw_copy,
    .decrypt = gmkw_decrypt_hook, .report = gmkw_report,
    .report_verbose = gmkw_report_verbose,
};

// ===================================================================
//  Entry point
// ===================================================================

void solve_gromark(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Gromark needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 8) {
        printf("\n\nERROR: ciphertext too short for a Gromark solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: ciphertext has a non-alphabet symbol at position %d; "
                   "Gromark ciphertext must be solid letters.\n\n", i);
            return;
        }

    int variant = cfg->cipher_type;            // GROMARK or GROMARK_PERIODIC

    // Cribs are not used (the offset/primer coupling makes a soft crib give no gradient).
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    // --- Periodic Gromark: anneal the KEYWORD directly (one config per swept period) ---
    if (variant == GROMARK_PERIODIC) {
        GromarkKwScratch kw;
        kw.n = cipher_len;
        if (cfg->period_present) { kw.minP = kw.maxP = cfg->period; }
        else { kw.minP = 4; kw.maxP = (cfg->max_period > 0) ? cfg->max_period : 8; if (kw.maxP > 9) kw.maxP = 9; }
        if (cfg->verbose)
            printf("\ngromark (periodic): %d letters, keyword anneal, period(s) %d..%d\n",
                cipher_len, kw.minP, kw.maxP);
        SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
            cipher_indices, cipher_len, crib_indices, crib_positions, 0);
        ctx.model_scratch = &kw;
        ctx.result = result;
        run_solver(&GROMARK_KW_MODEL, &ctx);
        return;
    }

    // --- Basic Gromark: primer pre-pass over the 10^P space, then anneal sigma per top-K ---
    int fixed_period = cfg->period_present ? cfg->period : GROMARK_PRIMER_LEN;
    int minP = fixed_period, maxP = fixed_period;

    // Length-adaptive K: more candidate primers for shorter (harder) text. Each top-K primer
    // is a full alphabet-anneal, so cap it. -nprimers overrides.
    int K;
    if (cfg->n_primers > 0) {
        K = cfg->n_primers;
    } else {
        if (cipher_len >= 200) K = 24;
        else if (cipher_len >= 150) K = 48;
        else if (cipher_len >= 120) K = 96;
        else K = GM_TOPK_MAX;
    }
    if (K > GM_TOPK_MAX) K = GM_TOPK_MAX;

    if (cfg->verbose)
        printf("\ngromark (basic): %d letters, primer space 10^%d, top-K=%d\n",
            cipher_len, fixed_period, K);

    int n_primers = gromark_rank_primers(cipher_indices, cipher_len, variant,
        fixed_period, minP, maxP, shared->ngram_data, cfg->ngram_size, K,
        g_gm_periods, g_gm_primers, g_gm_warm, cfg->verbose);
    if (n_primers < 1) {
        printf("\n\nERROR: Gromark primer pre-pass produced no candidates.\n\n");
        return;
    }

    GromarkScratch scratch;
    scratch.n = cipher_len;
    scratch.n_primers = n_primers;
    scratch.periods = g_gm_periods;
    scratch.primers = g_gm_primers;
    scratch.warm = g_gm_warm;
    g_gm_chain_cfg = -1;                        // invalidate the running-key cache

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&GROMARK_MODEL, &ctx);
}

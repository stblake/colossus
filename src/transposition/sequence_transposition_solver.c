#include "sequence_transposition_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Sequence Transposition solver (TYPE sequence-transposition / seqtrans / st)
// =====================================================================
//
// Sequence Transposition is a COLUMNAR transposition whose column assignment comes from a
// Gromark-style chain-addition digit sequence: a P-digit primer generates one digit 0..9 per
// plaintext letter (SS[i]), the letter goes into column SS[i], and the 10 columns are drawn off
// in a keyword-defined order. See sequence_transposition.c for the primitive.
//
// The cryptanalytic unknowns are:
//   * the PRIMER (default 5 digits, 10^5 space). In the ACA convention the primer is transmitted
//     openly at the head of the message, so it is normally KNOWN. When it is supplied (-primer)
//     the solver skips straight to the transposition; when it is not, a Gromark-style pre-pass
//     ranks the whole primer space and keeps the top-K as engine configs.
//   * the READ-ORDER PERMUTATION of the 10 digit columns (10! = 3.6M). This is a genuine
//     transposition with gradient -- swapping the read order of two columns moves a contiguous
//     block of the decrypt -- so it is a small permutation climb (state = pi[0..9], a permutation
//     of 0..9; pi[k] = the digit column emitted at read step k). The keyword itself is NOT
//     uniquely recoverable (many keywords share a ranking), so the solver reports the read order.
//
// For a fixed primer the digit sequence SS[] and per-column counts are constant, so they are
// cached per config (thread-local); a decrypt is then O(N): lay the columns out in pi order to
// get each column's start offset, then redistribute the ciphertext to plaintext positions.
//
// SHAPE_ANNEAL over pi (like the columnar solver), so -method shotgun/anneal/pso all apply.

#define SEQ_BUCKETS   SEQ_TRANS_BUCKETS      // 10 digit columns (0..9)
#define SEQ_TOPK_MAX  256                    // cap on kept primer configs (blind case)

// Blind pre-pass mini-solve budget (per candidate primer): a few restarts of a short swap
// hill-climb over the 10-element read order. Deliberately cheap -- it only RANKS primers; the
// kept top-K are searched thoroughly by the engine afterwards.
#define SEQ_PREPASS_RESTARTS 3
#define SEQ_PREPASS_ITERS    80

typedef struct {
    int  n;                                  // cipher length (== plaintext length)
    int  primer_len;                         // P (chain-addition lag / primer digit count)
    int  n_primers;                          // number of configs
    const int *primers;                      // [n_primers * SEQ_TRANS_MAX_PRIMER] digits per config
} SeqScratch;

// Per-config cache of the digit sequence + column counts (rebuilt when the config index changes).
static _Thread_local int g_seq_ss[MAX_CIPHER_LENGTH];
static _Thread_local int g_seq_count[SEQ_BUCKETS];
static _Thread_local int g_seq_cfg = -1;

// Pre-pass / single-primer store (single-threaded; filled before the search starts).
static int g_seq_primers[SEQ_TOPK_MAX * SEQ_TRANS_MAX_PRIMER];

// Build SS[] (one digit per position) and the per-column counts for a primer.
static void seq_build_ss(const int primer[], int P, int n, int ss[], int count[]) {
    gromark_chain_key(primer, P, n, ss);
    for (int d = 0; d < SEQ_BUCKETS; d++) count[d] = 0;
    for (int i = 0; i < n; i++) count[ss[i]]++;
}

// Decrypt from a cached (ss, count): each column's start offset is its cumulative position in
// the pi read order; redistribute the ciphertext back to plaintext positions. O(n).
static void seq_decrypt_from_cache(const int cipher[], int n, const int ss[],
                                   const int count[], const int pi[], int out[]) {
    int next[SEQ_BUCKETS], run = 0;
    for (int k = 0; k < SEQ_BUCKETS; k++) { int d = pi[k]; next[d] = run; run += count[d]; }
    for (int i = 0; i < n; i++) out[i] = cipher[next[ss[i]]++];
}

// ===================================================================
//  Blind primer pre-pass
// ===================================================================

// Quick best-read-order score for a fixed primer: a few restarts of a swap hill-climb over the
// 10-element read order, returning the best n-gram score (and, optionally, the best pi). Used
// only to RANK primers -- the kept top-K are annealed properly by the engine.
static double seq_best_pi(const int cipher[], int n, const int ss[], const int count[],
                          const float *ng, int ngsz, int restarts, int iters, int best_pi_out[]) {
    static _Thread_local int dec[MAX_CIPHER_LENGTH];
    int pi[SEQ_BUCKETS], best_pi[SEQ_BUCKETS];
    double best = -1e30;
    for (int r = 0; r < restarts; r++) {
        perm_seed(pi, SEQ_BUCKETS);
        seq_decrypt_from_cache(cipher, n, ss, count, pi, dec);
        double cur = ngram_score(dec, n, (float *) ng, ngsz);
        for (int it = 0; it < iters; it++) {
            int a = rand_int(0, SEQ_BUCKETS), b = rand_int(0, SEQ_BUCKETS);
            if (a == b) continue;
            int t = pi[a]; pi[a] = pi[b]; pi[b] = t;                    // swap
            seq_decrypt_from_cache(cipher, n, ss, count, pi, dec);
            double s = ngram_score(dec, n, (float *) ng, ngsz);
            if (s > cur) cur = s;                                       // greedy accept
            else { t = pi[a]; pi[a] = pi[b]; pi[b] = t; }               // revert
        }
        if (cur > best) { best = cur; for (int k = 0; k < SEQ_BUCKETS; k++) best_pi[k] = pi[k]; }
    }
    if (best_pi_out) for (int k = 0; k < SEQ_BUCKETS; k++) best_pi_out[k] = best_pi[k];
    return best;
}

int seq_trans_rank_primers(const int cipher[], int n, int primer_len,
                           const float *ngram, int ngram_size, int K, int out_primers[]) {
    int P = primer_len;
    if (K > SEQ_TOPK_MAX) K = SEQ_TOPK_MAX;
    if (K < 1) K = 1;

    static double top_score[SEQ_TOPK_MAX];
    static int    top_primer[SEQ_TOPK_MAX * SEQ_TRANS_MAX_PRIMER];
    int count_kept = 0, worst = 0;

    long space = 1;
    for (int i = 0; i < P; i++) space *= 10;                            // 10^P
    int primer[SEQ_TRANS_MAX_PRIMER];
    static _Thread_local int ss[MAX_CIPHER_LENGTH];
    int cnt[SEQ_BUCKETS];

    for (long val = 0; val < space; val++) {
        long x = val;
        for (int i = P - 1; i >= 0; i--) { primer[i] = (int)(x % 10); x /= 10; }
        seq_build_ss(primer, P, n, ss, cnt);
        double sc = seq_best_pi(cipher, n, ss, cnt, ngram, ngram_size,
                                SEQ_PREPASS_RESTARTS, SEQ_PREPASS_ITERS, NULL);

        int slot = -1;
        if (count_kept < K) {
            slot = count_kept++;
        } else if (sc > top_score[worst]) {
            slot = worst;
        }
        if (slot >= 0) {
            top_score[slot] = sc;
            for (int i = 0; i < P; i++) top_primer[slot * SEQ_TRANS_MAX_PRIMER + i] = primer[i];
            if (count_kept == K) {                                      // refresh the worst slot
                worst = 0;
                for (int i = 1; i < K; i++) if (top_score[i] < top_score[worst]) worst = i;
            }
        }
    }

    // Sort kept entries by score descending (selection; K small).
    int order[SEQ_TOPK_MAX];
    for (int i = 0; i < count_kept; i++) order[i] = i;
    for (int a = 0; a < count_kept; a++) {
        int b = a;
        for (int c = a + 1; c < count_kept; c++) if (top_score[order[c]] > top_score[order[b]]) b = c;
        int t = order[a]; order[a] = order[b]; order[b] = t;
    }
    for (int k = 0; k < count_kept; k++) {
        int s = order[k];
        for (int i = 0; i < SEQ_TRANS_MAX_PRIMER; i++)
            out_primers[k * SEQ_TRANS_MAX_PRIMER + i] =
                (i < P) ? top_primer[s * SEQ_TRANS_MAX_PRIMER + i] : 0;
    }
    return count_kept;
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int seq_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const SeqScratch *a = (const SeqScratch *) ctx->model_scratch;
    int n = a->n_primers;
    if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = a->primer_len;      // P
        out[i].j = i;                       // config index -> primer
        out[i].k = 0; out[i].aux[0] = 0; out[i].aux[1] = 0;
    }
    return n;
}

static void seq_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    perm_seed(st->key, SEQ_BUCKETS);        // random read order
    st->aux[0] = cc->j;                     // config index (-> primer)
    st->aux[1] = cc->period;                // P
    st->key_len = SEQ_BUCKETS;
}

static void seq_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                        bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    perm_move(st->key, SEQ_BUCKETS);        // swap-dominant permutation move (+ reverse/block-move)
}

static void seq_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int k = 0; k < SEQ_BUCKETS; k++) dst->key[k] = src->key[k];
    dst->aux[0] = src->aux[0];
    dst->aux[1] = src->aux[1];
    dst->key_len = src->key_len;
}

static void seq_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                             int *out, double *score_adjust) {
    (void) cc;
    const SeqScratch *a = (const SeqScratch *) ctx->model_scratch;
    int cfgidx = st->aux[0], P = st->aux[1], n = ctx->cipher_len;
    if (g_seq_cfg != cfgidx) {              // rebuild SS[]/counts for this (pinned) primer
        seq_build_ss(&a->primers[cfgidx * SEQ_TRANS_MAX_PRIMER], P, n, g_seq_ss, g_seq_count);
        g_seq_cfg = cfgidx;
    }
    seq_decrypt_from_cache(ctx->cipher, n, g_seq_ss, g_seq_count, st->key, out);
    *score_adjust = 0.0;
}

// ===================================================================
//  Reporting
// ===================================================================

static void seq_pi_string(const int pi[], char out[]) {
    int w = 0;
    for (int k = 0; k < SEQ_BUCKETS; k++) w += sprintf(out + w, "%s%d", k ? " " : "", pi[k]);
}

static void seq_primer_string(const int primer[], int P, char out[]) {
    for (int i = 0; i < P; i++) out[i] = (char)('0' + primer[i] % 10);
    out[P] = '\0';
}

static void seq_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) decrypted;
    const SeqScratch *a = (const SeqScratch *) ctx->model_scratch;
    int cfgidx = st->aux[0], P = st->aux[1];
    char primer[SEQ_TRANS_MAX_PRIMER + 1], pistr[4 * SEQ_BUCKETS];
    seq_primer_string(&a->primers[cfgidx * SEQ_TRANS_MAX_PRIMER], P, primer);
    seq_pi_string(st->key, pistr);
    double elapsed = engine_elapsed_sec(stats);
    printf("\n  primer %s, read order [%s], score=%.4f  [%.1fs, %d restarts]\n",
        primer, pistr, score, elapsed, stats->n_restarts);
    fflush(stdout);
}

static void seq_report(const SolverCtx *ctx, const SolverConfig *cc,
                       const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    const SeqScratch *a = (const SeqScratch *) ctx->model_scratch;
    int n = ctx->cipher_len, cfgidx = st->aux[0], P = st->aux[1];

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char primer[SEQ_TRANS_MAX_PRIMER + 1], pistr[4 * SEQ_BUCKETS];
    seq_primer_string(&a->primers[cfgidx * SEQ_TRANS_MAX_PRIMER], P, primer);
    seq_pi_string(st->key, pistr);

    printf("\nResult Score: %.2f | Words: %d | primer=%s | read order=[%s]\n",
        score, n_words_found, primer, pistr);

    print_text(ctx->cipher, n);
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

    // One-liner summary: >>> score, [words,] type, primer=, read=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, primer=%s, read=[%s], ",
            score, n_words_found, cfg->cipher_type, primer, pistr);
    else
        printf(">>> %.2f, %d, primer=%s, read=[%s], ",
            score, cfg->cipher_type, primer, pistr);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(ctx->cipher, n);
    printf(", ");
    print_text(decrypted, n);
    printf("\n");
}

static const CipherModel SEQUENCE_TRANSPOSITION_MODEL = {
    .name = "sequence-transposition", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = seq_enumerate, .key_len = NULL,
    .seed = seq_seed, .perturb = seq_perturb, .copy_state = seq_copy,
    .decrypt = seq_decrypt_hook, .report = seq_report,
    .report_verbose = seq_report_verbose,
};

// ===================================================================
//  Entry point
// ===================================================================

void solve_sequence_transposition(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;

    if (cipher_len < SEQ_BUCKETS) {
        printf("\n\nERROR: ciphertext too short for a Sequence Transposition solve.\n\n");
        return;
    }
    // Non-letter characters (spaces / punctuation) are carried as reversible negative
    // sentinels and ride through the transposition as ordinary cells: the chain-addition
    // sequence buckets each position regardless of its token, seq_decrypt only PERMUTES
    // tokens (their values are never used as array indices), and ngram_score projects to
    // the letters-only stream (sentinels transparent). So spaces/punctuation are preserved
    // -- they fall into place when the read order is right -- exactly like the other
    // space-significant transposition types. Pass -skipspaces to strip them upstream and
    // solve the solid-letter stream instead. A live in-alphabet value >= g_alpha would be
    // a decode bug, so keep guarding that.
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Sequence Transposition: symbol out of alphabet at "
                   "position %d.\n\n", i);
            return;
        }

    int P = (cfg->seq_primer_len > 0) ? cfg->seq_primer_len
          : (cfg->period_present ? cfg->period : SEQ_TRANS_PRIMER_LEN);
    if (P < 2) P = 2;
    if (P > SEQ_TRANS_MAX_PRIMER) P = SEQ_TRANS_MAX_PRIMER;

    SeqScratch scratch;
    scratch.n = cipher_len;
    scratch.primer_len = P;
    scratch.primers = g_seq_primers;

    if (cfg->seq_primer_len > 0) {
        // Primer supplied (the ACA convention): a single config, pure read-order search.
        for (int i = 0; i < P; i++) g_seq_primers[i] = cfg->seq_primer[i] % 10;
        scratch.n_primers = 1;
        if (cfg->verbose)
            printf("\nsequence-transposition: %d letters, primer given, read-order search (10!)\n",
                cipher_len);
    } else {
        // Blind: rank the 10^P primer space, keep top-K as engine configs. Length-adaptive K
        // (more candidates for shorter, harder text); -nprimers overrides.
        int K;
        if (cfg->n_primers > 0)       K = cfg->n_primers;
        else if (cipher_len >= 200)   K = 24;
        else if (cipher_len >= 150)   K = 64;
        else if (cipher_len >= 120)   K = 128;
        else                          K = SEQ_TOPK_MAX;
        if (K > SEQ_TOPK_MAX) K = SEQ_TOPK_MAX;

        if (cfg->verbose)
            printf("\nsequence-transposition: %d letters, blind primer pre-pass "
                   "(space 10^%d, top-K=%d)\n", cipher_len, P, K);

        int nk = seq_trans_rank_primers(cipher_indices, cipher_len, P,
            shared->ngram_data, cfg->ngram_size, K, g_seq_primers);
        if (nk < 1) {
            printf("\n\nERROR: Sequence Transposition primer pre-pass produced no candidates.\n\n");
            return;
        }
        scratch.n_primers = nk;
    }

    g_seq_cfg = -1;                         // invalidate the per-config cache

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;

    run_solver(&SEQUENCE_TRANSPOSITION_MODEL, &ctx);
}

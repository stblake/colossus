#include "keyphrase_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Key Phrase solver -- partition anneal + inner beam-Viterbi. See keyphrase_solver.h.
// =====================================================================

typedef struct {
    int   *cipher;             // letters-only ciphertext stream (0..25), length nl
    int    nl;
    int    ngram_size;
    float *ngram;
    int    beam;

    int    K;                  // distinct ciphertext letters that occur
    int    obs[ALPHABET_SIZE]; // group index -> the observed ciphertext letter (0..25)
    int    ct_group[ALPHABET_SIZE]; // ciphertext letter -> group index 0..K-1 (-1 if never occurs)

    int    eng_rank[ALPHABET_SIZE]; // plaintext letters by descending English monogram (warm start)
    int    ctf_rank[ALPHABET_SIZE]; // groups 0..K-1 by descending ciphertext frequency (warm start)

    int   *orig;               // original layout (letters + sentinels) for the spaced report
    int    orig_len;
} KeyphraseScratch;

// Decode scratch (per-worker), off the stack. g_kp_dec doubles as the report's decode buffer.
static _Thread_local unsigned char g_kp_sym[MAX_CIPHER_LENGTH][KP_BEAM_MAX];   // chosen letter per (pos, beam)
static _Thread_local unsigned char g_kp_back[MAX_CIPHER_LENGTH][KP_BEAM_MAX];  // parent beam index
static _Thread_local int g_kp_dec[MAX_CIPHER_LENGTH];

// ---------------------------------------------------------------------
//  Key moves (non-empty-preserving: every observed ct letter keeps >= 1 plaintext letter)
// ---------------------------------------------------------------------

static void kp_group_counts(const int grp[KEYPHRASE_LEN], int cnt[ALPHABET_SIZE], int K) {
    for (int g = 0; g < K; g++) cnt[g] = 0;
    for (int p = 0; p < KEYPHRASE_LEN; p++) cnt[grp[p]]++;
}

void keyphrase_move(int grp[KEYPHRASE_LEN], int K) {
    if (K <= 1) return;
    if (frand() < 0.65) {
        // Reassign one plaintext letter to a different group, keeping its source group non-empty.
        int cnt[ALPHABET_SIZE]; kp_group_counts(grp, cnt, K);
        int p, tries = 0;
        do { p = rand_int(0, KEYPHRASE_LEN); tries++; } while (cnt[grp[p]] <= 1 && tries < 40);
        if (cnt[grp[p]] <= 1) {                                   // all singletons: fall back to a swap
            int i = rand_int(0, KEYPHRASE_LEN), j = rand_int(0, KEYPHRASE_LEN);
            int t = grp[i]; grp[i] = grp[j]; grp[j] = t;
            return;
        }
        int g; do { g = rand_int(0, K); } while (g == grp[p]);
        grp[p] = g;
    } else {
        // Swap two plaintext letters' groups (fine-grained exploration; preserves group sizes).
        int i = rand_int(0, KEYPHRASE_LEN), j = rand_int(0, KEYPHRASE_LEN);
        int t = grp[i]; grp[i] = grp[j]; grp[j] = t;
    }
}

// ---------------------------------------------------------------------
//  Inner decode: beam Viterbi over the ambiguous letters-only stream
// ---------------------------------------------------------------------
//
// Each ciphertext letter maps (via its group) to a set of candidate plaintext letters; we choose
// one per position to maximize the n-gram score. A window is scored once, when the letter
// completing it is placed (from the n-th letter on), matching ngram_score's window walk. A small
// beam keeps the top-`beam` partial decodes per position. Fills out[0..nl-1], returns nl.
static int kp_decode(const KeyphraseScratch *a, const int grp[KEYPHRASE_LEN], int beam, int out[]) {
    const int n = a->ngram_size, nl = a->nl;
    const float *ngram = a->ngram;
    if (beam > KP_BEAM_MAX) beam = KP_BEAM_MAX;

    long POW = 1; for (int i = 0; i < n - 1; i++) POW *= ALPHABET_SIZE;   // 26^(n-1)

    // group -> candidate plaintext letters (rebuilt from grp[] each decode).
    int cand[ALPHABET_SIZE][ALPHABET_SIZE], ncand[ALPHABET_SIZE];
    for (int g = 0; g < a->K; g++) ncand[g] = 0;
    for (int p = 0; p < KEYPHRASE_LEN; p++) {
        int g = grp[p];
        if (g < 0 || g >= a->K) continue;                        // defensive
        cand[g][ncand[g]++] = p;
    }

    long   cur_key[KP_BEAM_MAX];
    double cur_score[KP_BEAM_MAX];
    int    cur_cnt = 1;
    cur_key[0] = 0; cur_score[0] = 0.0;

    for (int t = 0; t < nl; t++) {
        int c = a->cipher[t];
        int g = a->ct_group[c];
        int nc = (g >= 0) ? ncand[g] : 0;
        int addwin = (t >= n - 1);
        if (nc == 0) {                                           // empty group (shouldn't occur): filler A
            for (int b = 0; b < cur_cnt; b++) {
                g_kp_sym[t][b] = 0; g_kp_back[t][b] = (unsigned char) b;
                cur_key[b] = (cur_key[b] * ALPHABET_SIZE + 0) % POW;
            }
            continue;
        }
        long   nk_top[KP_BEAM_MAX];
        double ns_top[KP_BEAM_MAX];
        int    from_top[KP_BEAM_MAX], sym_top[KP_BEAM_MAX], top = 0;
        for (int b = 0; b < cur_cnt; b++) {
            for (int ci = 0; ci < nc; ci++) {
                int x = cand[g][ci];
                long widx = cur_key[b] * ALPHABET_SIZE + x;
                double ns = cur_score[b] + (addwin ? ngram_weight_at(ngram, (int) widx) : 0.0);
                long nk = widx % POW;
                if (top < beam) {
                    int p = top++;
                    while (p > 0 && ns_top[p - 1] < ns) {
                        ns_top[p] = ns_top[p - 1]; nk_top[p] = nk_top[p - 1];
                        from_top[p] = from_top[p - 1]; sym_top[p] = sym_top[p - 1]; p--;
                    }
                    ns_top[p] = ns; nk_top[p] = nk; from_top[p] = b; sym_top[p] = x;
                } else if (ns > ns_top[beam - 1]) {
                    int p = beam - 1;
                    while (p > 0 && ns_top[p - 1] < ns) {
                        ns_top[p] = ns_top[p - 1]; nk_top[p] = nk_top[p - 1];
                        from_top[p] = from_top[p - 1]; sym_top[p] = sym_top[p - 1]; p--;
                    }
                    ns_top[p] = ns; nk_top[p] = nk; from_top[p] = b; sym_top[p] = x;
                }
            }
        }
        for (int b = 0; b < top; b++) {
            cur_key[b] = nk_top[b]; cur_score[b] = ns_top[b];
            g_kp_sym[t][b]  = (unsigned char) sym_top[b];
            g_kp_back[t][b] = (unsigned char) from_top[b];
        }
        cur_cnt = top;
    }

    int best = 0; for (int b = 1; b < cur_cnt; b++) if (cur_score[b] > cur_score[best]) best = b;
    int b = best;
    for (int t = nl - 1; t >= 0; t--) { out[t] = g_kp_sym[t][b]; b = g_kp_back[t][b]; }
    return nl;
}

// ---------------------------------------------------------------------
//  CipherModel hooks
// ---------------------------------------------------------------------

static int kp_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = KEYPHRASE_LEN; out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Frequency warm start: give each group its highest-frequency plaintext letter first (so every
// group is non-empty), then distribute the rest rank-proportionally to ciphertext frequency, so
// the assigned plaintext mass roughly matches each ct letter's frequency. Jitter for diversity.
static void kp_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) cc;
    const KeyphraseScratch *a = (const KeyphraseScratch *) ctx->model_scratch;
    int K = a->K;
    for (int r = 0; r < K; r++) st->key[a->eng_rank[r]] = a->ctf_rank[r];
    for (int r = K; r < KEYPHRASE_LEN; r++) {
        int gi = (int) ((long) r * K / KEYPHRASE_LEN); if (gi >= K) gi = K - 1;
        st->key[a->eng_rank[r]] = a->ctf_rank[gi];
    }
    st->key_len = KEYPHRASE_LEN;
    int kicks = rand_int(0, 7);
    for (int k = 0; k < kicks; k++) keyphrase_move(st->key, K);
}

static void kp_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st, bool *force_primary) {
    (void) cc; (void) force_primary;
    const KeyphraseScratch *a = (const KeyphraseScratch *) ctx->model_scratch;
    keyphrase_move(st->key, a->K);
    if (frand() < 0.10) keyphrase_move(st->key, a->K);
}

static void kp_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < KEYPHRASE_LEN; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void kp_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                            int *out, double *score_adjust) {
    (void) cc;
    const KeyphraseScratch *a = (const KeyphraseScratch *) ctx->model_scratch;
    kp_decode(a, st->key, a->beam, out);
    *score_adjust = 0.0;
}

// ---------------------------------------------------------------------
//  Reporting
// ---------------------------------------------------------------------

// Walk the original layout, emitting each sentinel verbatim and each letter from `letters[]`
// (restoring the ciphertext's word divisions around the recovered plaintext / ciphertext).
static void kp_print_spaced(const int orig[], int orig_len, const int letters[]) {
    int c = 0;
    for (int i = 0; i < orig_len; i++)
        putchar(orig[i] >= 0 ? index_to_char(letters[c++]) : index_to_char(orig[i]));
}

static void kp_report_verbose(const SolverCtx *ctx, const SolverConfig *cc, const SolverState *st,
        double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) st;
    const KeyphraseScratch *a = (const KeyphraseScratch *) ctx->model_scratch;
    char params[32];
    snprintf(params, sizeof params, "K=%d", a->K);
    report_transposition_verbose(ctx, score, decrypted, stats, params);
}

static void kp_report(const SolverCtx *ctx, const SolverConfig *cc, const SolverState *st,
                      double score, int *decrypted) {
    (void) cc; (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    const KeyphraseScratch *a = (const KeyphraseScratch *) ctx->model_scratch;
    int nl = a->nl;

    kp_decode(a, st->key, KP_FINAL_BEAM, g_kp_dec);              // re-decode with the wide final beam

    int n_words_found = 0;
    static char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < nl; i++) plaintext_string[i] = index_to_char(g_kp_dec[i]);
    plaintext_string[nl] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    // Recovered key phrase (plaintext a..z -> ciphertext letter = obs[group]).
    char phrase[KEYPHRASE_LEN + 1];
    for (int p = 0; p < KEYPHRASE_LEN; p++) phrase[p] = index_to_char(a->obs[st->key[p]]);
    phrase[KEYPHRASE_LEN] = '\0';

    printf("\nResult Score: %.2f | Words: %d | K=%d | keyphrase=%s\n",
        score, n_words_found, a->K, phrase);
    kp_print_spaced(a->orig, a->orig_len, a->cipher); printf("\n");
    kp_print_spaced(a->orig, a->orig_len, g_kp_dec);  printf("\n");
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(g_kp_dec, ctx->result->decrypted, nl);
        ctx->result->decrypted_len = nl;
    }

    // >>> score, [words,] type, keyphrase=, K=, file, CIPHER, PLAINTEXT (letters-only).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, keyphrase=%s, K=%d, ", score, n_words_found, cfg->cipher_type, phrase, a->K);
    else
        printf(">>> %.2f, %d, keyphrase=%s, K=%d, ", score, cfg->cipher_type, phrase, a->K);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(a->cipher, nl);
    printf(", ");
    print_text(g_kp_dec, nl);
    printf("\n");
}

static const CipherModel KEYPHRASE_MODEL = {
    .name = "keyphrase", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = kp_enumerate, .key_len = NULL,
    .seed = kp_seed, .perturb = kp_perturb, .copy_state = kp_copy,
    .decrypt = kp_decrypt_hook, .report = kp_report,
    .report_verbose = kp_report_verbose,
};

// ---------------------------------------------------------------------
//  Entry point
// ---------------------------------------------------------------------

// Letters-only search stream + the original layout (spaces/punct), off the stack.
static int g_kp_letters[MAX_CIPHER_LENGTH];
static int g_kp_orig[MAX_CIPHER_LENGTH];

void solve_keyphrase(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != DEFAULT_ALPHABET_SIZE) {
        printf("\n\nERROR: Key Phrase needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    int nl = 0;
    for (int i = 0; i < cipher_len; i++) {
        g_kp_orig[i] = cipher_indices[i];
        if (cipher_indices[i] >= 0) g_kp_letters[nl++] = cipher_indices[i];
    }
    if (nl < 8) {
        printf("\n\nERROR: ciphertext too short for a Key Phrase solve (%d letters).\n\n", nl);
        return;
    }

    KeyphraseScratch scratch;
    memset(&scratch, 0, sizeof scratch);

    // Observed ciphertext letters -> groups (label = the ct letter, read off the ciphertext).
    int hist[ALPHABET_SIZE];
    for (int c = 0; c < ALPHABET_SIZE; c++) { hist[c] = 0; scratch.ct_group[c] = -1; }
    for (int i = 0; i < nl; i++) hist[g_kp_letters[i]]++;
    int K = 0;
    for (int c = 0; c < ALPHABET_SIZE; c++)
        if (hist[c] > 0) { scratch.obs[K] = c; scratch.ct_group[c] = K; K++; }
    if (K < 2) {
        printf("\n\nERROR: ciphertext uses only %d distinct letters -- too degenerate.\n\n", K);
        return;
    }

    scratch.cipher = g_kp_letters; scratch.nl = nl;
    scratch.ngram = shared->ngram_data; scratch.ngram_size = cfg->ngram_size;
    scratch.beam = KP_BEAM;
    scratch.orig = g_kp_orig; scratch.orig_len = cipher_len;
    scratch.K = K;

    // Warm-start ranks: plaintext letters by descending English monogram; groups by descending
    // ciphertext frequency (selection sort -- 26 / K elements).
    for (int i = 0; i < KEYPHRASE_LEN; i++) scratch.eng_rank[i] = i;
    for (int x = 0; x < KEYPHRASE_LEN; x++)
        for (int y = x + 1; y < KEYPHRASE_LEN; y++)
            if (g_monograms[scratch.eng_rank[y]] > g_monograms[scratch.eng_rank[x]]) {
                int t = scratch.eng_rank[x]; scratch.eng_rank[x] = scratch.eng_rank[y]; scratch.eng_rank[y] = t;
            }
    for (int g = 0; g < K; g++) scratch.ctf_rank[g] = g;
    for (int x = 0; x < K; x++)
        for (int y = x + 1; y < K; y++)
            if (hist[scratch.obs[scratch.ctf_rank[y]]] > hist[scratch.obs[scratch.ctf_rank[x]]]) {
                int t = scratch.ctf_rank[x]; scratch.ctf_rank[x] = scratch.ctf_rank[y]; scratch.ctf_rank[y] = t;
            }

    if (cfg->verbose)
        printf("\nkeyphrase: %d letters, %d distinct ciphertext letters, beam-Viterbi partition anneal\n", nl, K);

    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_kp_letters, nl, crib_indices, crib_positions, 0);
    ctx.model_scratch = &scratch;
    ctx.result = result;
    run_solver(&KEYPHRASE_MODEL, &ctx);
}

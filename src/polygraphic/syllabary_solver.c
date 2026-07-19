#include "syllabary_solver.h"
#include "syllabary.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// =====================================================================
//  Syllabary solver (TYPE syllabary) -- see syllabary_solver.h
// =====================================================================
//
// The key is a permutation of the 100 syllabary tokens over the 100 codes (indexed by code
// 0..99); decode is code -> token -> its letters, concatenated. Because tokens are 1-3 letters
// the decode length varies with the key, so -- like Fractionated Morse -- the decrypt hook
// cyclically TILES the variable decode into a fixed-length scoring buffer C = SYL_MAXTOK*ncodes
// (>= any decode, so tiling repeats rather than truncates), making the mean n-gram length-fair.
// The search is a plain SHAPE_ANNEAL two-code swap on the generic engine path (no incremental
// fast path: the length change breaks the per-position window bookkeeping). Seed is a frequency
// warm start (frequent codes -> common tokens). Cribs are not used (positions do not align).

#define SYL_MAX_CODES 100

// Parsed ciphertext + fixed scoring buffer, filled by solve_syllabary (read-only in the search).
static int g_syl_codes[MAX_CIPHER_LENGTH];             // the ncodes CT codes (0..99)
static int g_syl_cipher[MAX_CIPHER_LENGTH];            // dummy engine "cipher" of length C (zeros)
static _Thread_local int g_syl_decode[MAX_CIPHER_LENGTH];  // per-worker clean decode scratch

typedef struct {
    int    ncodes;                       // number of CT codes (positions)
    int    C;                            // fixed scoring length = SYL_MAXTOK * ncodes (<= MAX)
    int    code_rank[100];               // codes sorted by descending observed frequency
    int    tok_rank[SYL_NTOKENS];        // tokens sorted by descending heuristic commonness
} SyllabaryScratch;

// Concatenate the decode of the current key into g_syl_decode[]; return its letter length.
static int syl_decode_key(const SyllabaryScratch *s, const int *key) {
    int L = 0;
    for (int i = 0; i < s->ncodes; i++) {
        const char *tok = syllabary_tokens[key[g_syl_codes[i]]];
        for (int k = 0; tok[k] && L < MAX_CIPHER_LENGTH; k++)
            if (tok[k] >= 'A' && tok[k] <= 'Z') g_syl_decode[L++] = tok[k] - 'A';
    }
    return L;
}

static int syl_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = SYL_NTOKENS; out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

// Frequency warm start: assign the most-frequent codes the most-common tokens (ranked once in
// solve), then a few random swaps so shotgun restarts diversify.
static void syl_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const SyllabaryScratch *s = (const SyllabaryScratch *) ctx->model_scratch;
    (void) cc;
    for (int i = 0; i < SYL_NTOKENS; i++) st->key[s->code_rank[i]] = s->tok_rank[i];
    st->key_len = SYL_NTOKENS;
    int kicks = rand_int(0, 40);
    for (int k = 0; k < kicks; k++) {
        int a = rand_int(0, SYL_NTOKENS), b = rand_int(0, SYL_NTOKENS);
        int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
    }
}

static void syl_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                        SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    int a = rand_int(0, SYL_NTOKENS), b = rand_int(0, SYL_NTOKENS);
    int t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
    if (frand() < 0.10) {
        a = rand_int(0, SYL_NTOKENS); b = rand_int(0, SYL_NTOKENS);
        t = st->key[a]; st->key[a] = st->key[b]; st->key[b] = t;
    }
}

static void syl_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < SYL_NTOKENS; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

static void syl_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    const SyllabaryScratch *s = (const SyllabaryScratch *) ctx->model_scratch;
    int C = ctx->cipher_len;
    int L = syl_decode_key(s, st->key);
    if (L <= 0) { for (int i = 0; i < C; i++) out[i] = 0; *score_adjust = 0.0; return; }
    for (int i = 0; i < C; i++) out[i] = g_syl_decode[i % L];   // tile to the fixed length C
    *score_adjust = 0.0;
}

static void syl_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) cc; (void) st; (void) decrypted;
    const SyllabaryScratch *s = (const SyllabaryScratch *) ctx->model_scratch;
    char params[48];
    snprintf(params, sizeof(params), "codes=%d", s->ncodes);
    report_transposition_verbose(ctx, score, decrypted, stats, params);
}

// Lay the recovered code -> token map into a 10x10 by code (row = code/10, col = code%10). The
// label order is not identifiable, so this is a canonical arrangement, not the true keysquare.
static void syl_print_square(const int *key) {
    printf("     0    1    2    3    4    5    6    7    8    9\n");
    for (int r = 0; r < 10; r++) {
        printf("  %d ", r);
        for (int c = 0; c < 10; c++) {
            const char *tok = syllabary_tokens[key[r * 10 + c]];
            printf(" %-4s", tok[0] ? tok : ".");
        }
        printf("\n");
    }
}

static void syl_report(const SolverCtx *ctx, const SolverConfig *cc,
                       const SolverState *st, double score, int *decrypted) {
    (void) cc; (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    const SyllabaryScratch *s = (const SyllabaryScratch *) ctx->model_scratch;
    int L = syl_decode_key(s, st->key);
    if (L > MAX_CIPHER_LENGTH) L = MAX_CIPHER_LENGTH;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < L; i++) plaintext_string[i] = index_to_char(g_syl_decode[i]);
    plaintext_string[L] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    printf("\nResult Score: %.2f | Words: %d | codes=%d decode=%d letters\n",
        score, n_words_found, s->ncodes, L);
    print_text(g_syl_decode, L);
    printf("\n%s\n", ctx->cribtext);
    printf("\nrecovered code -> token map (by code; label order not identifiable):\n");
    syl_print_square(st->key);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;
        vec_copy(g_syl_decode, ctx->result->decrypted, L);
        ctx->result->decrypted_len = L;
    }

    // One-liner: >>> score, [words,] type, codes=N, file, PLAINTEXT (the clean decode).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, codes=%d, ", score, n_words_found, cfg->cipher_type, s->ncodes);
    else
        printf(">>> %.2f, %d, codes=%d, ", score, cfg->cipher_type, s->ncodes);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(g_syl_decode, L);
    printf("\n");
}

static const CipherModel SYLLABARY_MODEL = {
    .name = "syllabary", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = syl_enumerate, .key_len = NULL,
    .seed = syl_seed, .perturb = syl_perturb, .copy_state = syl_copy,
    .decrypt = syl_decrypt_hook, .report = syl_report, .report_verbose = syl_report_verbose,
};

// Heuristic token commonness: single letters by monogram, syllables by summed letter monogram
// log (frequent-letter syllables rank high), digits/null rare. Ranks the warm-start target.
static double syl_token_commonness(int t) {
    const char *s = syllabary_tokens[t];
    int L = (int) strlen(s);
    if (L == 0) return -1e18;                          // null: rarest
    double v = 0.0;
    for (int k = 0; k < L; k++) {
        if (s[k] < 'A' || s[k] > 'Z') return -1e17;    // digit token: rare
        double f = g_monograms[s[k] - 'A'];
        v += log(f > 1e-9 ? f : 1e-9);
    }
    return v / L;                                      // per-letter, so lengths compare fairly
}

static int syl_parse_codes(const char *str, int out[], int cap) {
    int digits[2 * MAX_CIPHER_LENGTH], nd = 0;
    for (int i = 0; str[i] && nd < 2 * MAX_CIPHER_LENGTH; i++)
        if (str[i] >= '0' && str[i] <= '9') digits[nd++] = str[i] - '0';
    int n = 0;
    for (int i = 0; i + 1 < nd && n < cap; i += 2) out[n++] = digits[i] * 10 + digits[i + 1];
    return n;
}

void solve_syllabary(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_indices; (void) cipher_len;
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != DEFAULT_ALPHABET_SIZE) {
        printf("\n\nERROR: Syllabary needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    int ncodes = syl_parse_codes(ciphertext_str, g_syl_codes, MAX_CIPHER_LENGTH);
    if (ncodes < 4) {
        printf("\n\nERROR: parsed only %d Syllabary codes; need >= 4. The ciphertext must be a "
               "stream of 2-digit numbers.\n\n", ncodes);
        return;
    }
    // Fixed scoring length: SYL_MAXTOK * ncodes bounds any decode; clamp to the buffer.
    int C = SYL_MAXTOK * ncodes;
    if (C > MAX_CIPHER_LENGTH) C = MAX_CIPHER_LENGTH;

    SyllabaryScratch scratch;
    memset(&scratch, 0, sizeof scratch);
    scratch.ncodes = ncodes;
    scratch.C = C;

    // Rank codes by observed frequency (descending); unobserved codes keep their natural order.
    int freq[100];
    for (int i = 0; i < 100; i++) { freq[i] = 0; scratch.code_rank[i] = i; }
    for (int i = 0; i < ncodes; i++) freq[g_syl_codes[i]]++;
    for (int x = 0; x < 100; x++)
        for (int y = x + 1; y < 100; y++)
            if (freq[scratch.code_rank[y]] > freq[scratch.code_rank[x]]) {
                int t = scratch.code_rank[x]; scratch.code_rank[x] = scratch.code_rank[y]; scratch.code_rank[y] = t;
            }
    // Rank tokens by heuristic commonness (descending).
    double comm[SYL_NTOKENS];
    for (int i = 0; i < SYL_NTOKENS; i++) { comm[i] = syl_token_commonness(i); scratch.tok_rank[i] = i; }
    for (int x = 0; x < SYL_NTOKENS; x++)
        for (int y = x + 1; y < SYL_NTOKENS; y++)
            if (comm[scratch.tok_rank[y]] > comm[scratch.tok_rank[x]]) {
                int t = scratch.tok_rank[x]; scratch.tok_rank[x] = scratch.tok_rank[y]; scratch.tok_rank[y] = t;
            }

    if (cfg->verbose)
        printf("\nsyllabary: %d codes, scoring length C=%d (tiled), 100-token composite-map anneal\n",
               ncodes, C);

    for (int i = 0; i < C; i++) g_syl_cipher[i] = 0;      // dummy engine cipher (scoring reads out[])
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_syl_cipher, C, crib_indices, crib_positions, 0);
    ctx.result = result;
    ctx.model_scratch = &scratch;
    run_solver(&SYLLABARY_MODEL, &ctx);
}

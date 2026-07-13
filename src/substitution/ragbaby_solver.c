#include "ragbaby_solver.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Ragbaby solver (TYPE ragbaby)
// =====================================================================
//
// Ragbaby (ragbaby.c) shifts each plaintext LETTER forward by its word-position number (mod 24)
// in a keyed 24-letter alphabet KA (I/J and W/X paired). The per-letter numbering is a FIXED
// function of the (spaced) ciphertext, so the ONLY unknown is the keyed alphabet -- there is NO
// period (a single engine config).
//
// KEYED-ALPHABET SEARCH (not a free 24! permutation), exactly like Fractionated Morse / Digrafid.
// An ACA Ragbaby alphabet is a KEYED alphabet -- a keyword (duplicates dropped) then the rest of
// the 24-letter alphabet ascending. The state is a keyed-alphabet SEQUENCE (KA = st->key[0..23])
// maintained as "keyword prefix of length kw + ascending tail" and searched with the structure-
// preserving ragbaby_move_seq / ragbaby_canonicalize (the 24-cell twin of fracmorse_move_seq): ~4%
// grow/shrink kw, ~48% keyword<->tail swap (re-sort the tail), ~48% in-keyword reorder. kw lives in
// st->aux[0], sampled per restart from [RAG_KW_MIN..RAG_KW_MAX] -- RESTARTS cover the keyword length.
// This tracks the keyword the way an ACA solver does and cracks the short ACA (80-150-letter) range.
// (Since the per-letter shift is KNOWN, this is far more constrained than a free monoalphabetic
// solve; the KA is recoverable only up to a cyclic rotation, à la Playfair, but the plaintext is
// unique.) score_adjust stays 0 (every KA is a bijection). Cribs are NOT used (the letters-only /
// spaced positional mapping is fragile; Ragbaby is not crib-driven).

#define RAG_ALPHA   24               // 24-letter Ragbaby keyed alphabet
#define RAG_KW_MIN  3                // keyword-length search range (as fracmorse/digrafid)
#define RAG_KW_MAX  12               //   tail >= 12 cells stays sorted

// Parsed problem instance (single-threaded), off the stack: the folded cipher LETTERS and the
// parallel per-letter shift numbers, both length g_rag_len, filled by solve_ragbaby and read by
// the decrypt hook (ctx->cipher points at g_rag_cipher, ctx->cipher_len == g_rag_len).
static int g_rag_cipher[MAX_CIPHER_LENGTH];
static int g_rag_num[MAX_CIPHER_LENGTH];
static int g_rag_len;

// ===================================================================
//  Keyed-alphabet move (24-cell twin of fracmorse_canonicalize / _move_seq)
// ===================================================================

void ragbaby_canonicalize(int *seq, int kw) {
    char used[RAG_ALPHA];
    for (int i = 0; i < RAG_ALPHA; i++) used[i] = 0;
    for (int i = 0; i < kw; i++) used[seq[i]] = 1;
    int m = kw;
    for (int s = 0; s < RAG_ALPHA && m < RAG_ALPHA; s++) if (!used[s]) seq[m++] = s;
}

void ragbaby_move_seq(int *seq, int *kw) {
    double r = frand();
    if (r < 0.04) {                                   // grow / shrink the keyword length
        int k = *kw;
        if (k <= RAG_KW_MIN) k++;
        else if (k >= RAG_KW_MAX) k--;
        else k += (frand() < 0.5) ? 1 : -1;
        *kw = k;
        ragbaby_canonicalize(seq, k);
    } else if (r < 0.52) {                            // keyword <-> tail (coarse set search)
        int i = rand_int(0, *kw), j = rand_int(*kw, RAG_ALPHA);
        int t = seq[i]; seq[i] = seq[j]; seq[j] = t;
        ragbaby_canonicalize(seq, *kw);              // keep the tail sorted
    } else {                                          // in-keyword reorder (smooth)
        int i = rand_int(0, *kw), j = rand_int(0, *kw);
        int t = seq[i]; seq[i] = seq[j]; seq[j] = t;
    }
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int ragbaby_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = 0; out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

static void ragbaby_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx; (void) cc;
    int kw = rand_int(RAG_KW_MIN, RAG_KW_MAX + 1);
    random_keyword(st->key, RAG_ALPHA, kw);           // keyword prefix + ascending tail (= KA)
    st->aux[0] = kw;
    st->key_len = RAG_ALPHA;
}

static void ragbaby_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                            SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    ragbaby_move_seq(st->key, &st->aux[0]);
}

static void ragbaby_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < RAG_ALPHA; i++) dst->key[i] = src->key[i];
    dst->aux[0] = src->aux[0];
    dst->key_len = src->key_len;
}

// Decrypt: KA = st->key; ka_inv[KA[p]] = p; pt[i] = KA[(ka_inv[ct[i]] - num[i]) mod 24].
static void ragbaby_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                 SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    int ka_inv[RAG_ALPHA];
    for (int p = 0; p < RAG_ALPHA; p++) ka_inv[st->key[p]] = p;
    ragbaby_decrypt(ctx->cipher, g_rag_num, ctx->cipher_len, st->key, ka_inv, RAG_ALPHA, out);
    *score_adjust = 0.0;
}

// ===================================================================
//  Reporting
// ===================================================================

static void ragbaby_alpha_string(const int ka[], char out[]) {
    for (int i = 0; i < RAG_ALPHA; i++) out[i] = index_to_char(ka[i]);
    out[RAG_ALPHA] = '\0';
}

static void ragbaby_report(const SolverCtx *ctx, const SolverConfig *cc,
                           const SolverState *st, double score, int *decrypted) {
    (void) cc;
    ColossusConfig *cfg = ctx->cfg;
    int L = ctx->cipher_len;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < L; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[L] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char alpha[RAG_ALPHA + 1]; ragbaby_alpha_string(st->key, alpha);

    printf("\nResult Score: %.2f | Words: %d | alphabet=%s\n", score, n_words_found, alpha);

    print_cipher(ctx->cipher, L, NULL);
    printf("\n");
    print_text(decrypted, L);
    printf("\n");
    print_spaces_line(g_spaces_table, decrypted, L);
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;              // no period for Ragbaby
        vec_copy(decrypted, ctx->result->decrypted, L);
        ctx->result->decrypted_len = L;
    }

    // One-liner summary: >>> score, [words,] type, alphabet=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, alphabet=%s, ", score, n_words_found, cfg->cipher_type, alpha);
    else
        printf(">>> %.2f, %d, alphabet=%s, ", score, cfg->cipher_type, alpha);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, L, NULL);
    printf(", ");
    print_text(decrypted, L);
    printf("\n");
}

static const CipherModel RAGBABY_MODEL = {
    .name = "ragbaby", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = ragbaby_enumerate, .key_len = NULL,
    .seed = ragbaby_seed, .perturb = ragbaby_perturb, .copy_state = ragbaby_copy,
    .decrypt = ragbaby_decrypt_hook, .report = ragbaby_report,
};

// ===================================================================
//  Entry point
// ===================================================================

void solve_ragbaby(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_indices; (void) cipher_len; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != RAG_ALPHA) {
        printf("\n\nERROR: Ragbaby needs the 24-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    // Parse the spaced ciphertext directly (not the A..Z decode): the word divisions drive the
    // numbering. Folds J->I / X->W through g_char_to_idx; strips other punctuation.
    ragbaby_number_stream(ciphertext_str, RAG_ALPHA, g_rag_cipher, g_rag_num, &g_rag_len);
    if (g_rag_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Ragbaby solve (%d letters).\n\n", g_rag_len);
        return;
    }
    if (g_rag_len > MAX_CIPHER_LENGTH) g_rag_len = MAX_CIPHER_LENGTH;

    if (cfg->verbose)
        printf("\nragbaby: %d ciphertext letters, keyed-alphabet anneal (kw %d..%d)\n",
            g_rag_len, RAG_KW_MIN, RAG_KW_MAX);

    // Cribs are not used (the letters-only / spaced positional mapping is fragile).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_rag_cipher, g_rag_len, crib_indices, crib_positions, 0);
    ctx.result = result;

    run_solver(&RAGBABY_MODEL, &ctx);
}

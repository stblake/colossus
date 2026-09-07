#include "compressocrat_solver.h"
#include "compressocrat.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Compressocrat solver (TYPE compressocrat) -- ACA type 92
// =====================================================================
//
// Compressocrat (compressocrat.c) is the fractionation twin of Fractionated Morse: it writes the
// plaintext in a FIXED prefix-free code over the digits {1,2,3} (the "irregular" compression
// alphabet), pads the concatenated stream with the digit 1 to a multiple of 3, groups it into
// trigraphs, and maps each trigraph (rank 9(a-1)+3(b-1)+(c-1), 333 impossible) to a ciphertext
// letter through a keyed alphabet sigma (rank -> letter). The ONLY unknown is that 26-letter
// alphabet -- there is NO period (the trigraph grouping is fixed at 3), exactly as fracmorse.
//
// KEYED-ALPHABET SEARCH (not a free 26! permutation). An ACA Compressocrat alphabet is a KEYED
// ALPHABET -- a keyword (duplicates dropped) then the rest of the alphabet ascending -- as the
// generator builds it (random_keyword). So, following fracmorse/digrafid, the state is a keyed-
// alphabet SEQUENCE (sigma = st->key[0..25]) maintained as "keyword prefix of length kw +
// ascending tail" and searched with the SAME structure-preserving moves (compressocrat_move_seq
// / _canonicalize, the length-26 twin of fracmorse_move_seq): ~4% grow/shrink kw, ~48% keyword
// <-> tail swap (coarse set search, re-sort the tail), ~48% in-keyword reorder (smooth). kw lives
// in st->aux[0], sampled per restart from [CR_KW_MIN..CR_KW_MAX] -- RESTARTS cover the keyword
// length. This collapses the keyspace and makes the sorted tail a strong prior, so it tracks the
// keyword the way an ACA solver does and cracks the short (~110-150-letter) ACA puzzles.
//
// LENGTH HANDLING (the decode length varies per key, but the engine scores a fixed length). We
// decrypt to plaintext and score over a FIXED length C = the ciphertext length: compressocrat_
// decrypt yields nt letters (one per token, a filler for an invalid token), which are CYCLICALLY
// TILED to fill C so the mean n-gram score is length-fair across keys. (The decode EXPANDS here
// -- m ~ N >= C -- so the tile scores the first C ~ plaintext letters.) A structural VALIDITY
// reward (fraction of tokens that are legal codewords) is folded into score_adjust (the analogue
// of ADFGVX's IoC term): the true key decodes ~100% valid, so this gives the anneal a gradient
// toward clean-parse keys on short text. It effectively needs -logprob, like the fractionation
// family. Cribs are NOT used (the length change + tiling break the positional crib mapping).

#define CR_ALPHA        ALPHABET_SIZE    // the 26-letter keyed alphabet (== the 26 trigraphs)
#define CR_KW_MIN       3                // keyword-length search range (as fracmorse/Digrafid:
#define CR_KW_MAX       13               //   realistic ACA keyword lengths, tail >= 13 stays sorted)
#define CR_FILLER       23               // letter emitted for an invalid token ('X')
#define CR_VALID_WEIGHT 2.0              // validity reward weight (calibrated vs the solver test)

// Decoded-plaintext scratch (single-threaded), off the stack. A decode yields at most ~3C/2
// tokens (every token the 2-digit minimum), so 2*MAX_CIPHER_LENGTH is a safe bound for C <= MAX.
static _Thread_local int g_cr_decode[2 * MAX_CIPHER_LENGTH];

// ===================================================================
//  Keyed-alphabet move (length-26 twin of fracmorse_move_seq / _canonicalize)
// ===================================================================

// Restore the keyed-alphabet invariant after a keyword/tail boundary change: keep seq[0..kw-1]
// as the keyword (any order), then refill seq[kw..25] with the remaining letters ascending.
void compressocrat_canonicalize(int *seq, int kw) {
    char used[CR_ALPHA];
    for (int i = 0; i < CR_ALPHA; i++) used[i] = 0;
    for (int i = 0; i < kw; i++) used[seq[i]] = 1;
    int m = kw;
    for (int s = 0; s < CR_ALPHA && m < CR_ALPHA; s++) if (!used[s]) seq[m++] = s;
}

// One structure-preserving move on the 26-cell keyed alphabet `seq` (keyword prefix length *kw).
// The invariant (seq a permutation of 0..25, seq[kw..25] strictly ascending, kw in range) holds
// on exit. Same three classes / weights as fracmorse_move_seq (coarse fractionation landscape).
void compressocrat_move_seq(int *seq, int *kw) {
    double r = frand();
    if (r < 0.04) {                                   // grow / shrink the keyword length
        int k = *kw;
        if (k <= CR_KW_MIN) k++;
        else if (k >= CR_KW_MAX) k--;
        else k += (frand() < 0.5) ? 1 : -1;
        *kw = k;
        compressocrat_canonicalize(seq, k);
    } else if (r < 0.52) {                            // keyword <-> tail (coarse set search)
        int i = rand_int(0, *kw), j = rand_int(*kw, CR_ALPHA);
        int t = seq[i]; seq[i] = seq[j]; seq[j] = t;
        compressocrat_canonicalize(seq, *kw);         // keep the tail sorted
    } else {                                          // in-keyword reorder (smooth)
        int i = rand_int(0, *kw), j = rand_int(0, *kw);
        int t = seq[i]; seq[i] = seq[j]; seq[j] = t;
    }
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

// One config (no period, no starter). enumerate must return >= 1.
static int compressocrat_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    (void) ctx;
    if (cap < 1) return 0;
    out[0].period = 0; out[0].j = 0; out[0].k = 0; out[0].aux[0] = 0; out[0].aux[1] = 0;
    return 1;
}

static void compressocrat_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx; (void) cc;
    int kw = rand_int(CR_KW_MIN, CR_KW_MAX + 1);      // [MIN .. MAX]
    random_keyword(st->key, CR_ALPHA, kw);            // keyword prefix + ascending tail (= sigma)
    st->aux[0] = kw;
    st->key_len = CR_ALPHA;
}

static void compressocrat_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                                  SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    compressocrat_move_seq(st->key, &st->aux[0]);
}

static void compressocrat_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < CR_ALPHA; i++) dst->key[i] = src->key[i];
    dst->aux[0] = src->aux[0];
    dst->key_len = src->key_len;
}

// Decrypt into out[0..C-1]: decode to nt plaintext letters (filler for invalid tokens), then
// cyclically tile to fill C so the scored length is fixed and key-independent. score_adjust =
// the validity reward (fraction of tokens that are legal codewords).
static void compressocrat_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                       SolverState *st, int *out, double *score_adjust) {
    (void) cc;
    int C = ctx->cipher_len, nt = 0, nv = 0;
    int m = compressocrat_decrypt(ctx->cipher, C, st->key, g_cr_decode, CR_FILLER, &nt, &nv);
    if (m <= 0) {
        for (int i = 0; i < C; i++) out[i] = CR_FILLER;
        *score_adjust = 0.0;
        return;
    }
    for (int i = 0; i < C; i++) out[i] = g_cr_decode[i % m];     // tile to the fixed length C
    *score_adjust = CR_VALID_WEIGHT * ((nt > 0) ? (double) nv / nt : 0.0);
}

// ===================================================================
//  Reporting
// ===================================================================

static void compressocrat_alpha_string(const int sigma[], char out[]) {
    for (int i = 0; i < CR_ALPHA; i++) out[i] = index_to_char(sigma[i]);
    out[CR_ALPHA] = '\0';
}

// The engine's `decrypted` is the TILED length-C scoring buffer; the report re-decodes cleanly
// (at the solution every token is valid, so this is the real N-letter plaintext) and reports THAT
// -- the >>> CSV plaintext (which run_tests.sh compares to the .solution) and ctx->result carry
// the clean decode, not the tiled buffer.
static void compressocrat_report(const SolverCtx *ctx, const SolverConfig *cc,
                                 const SolverState *st, double score, int *decrypted) {
    (void) cc; (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    int C = ctx->cipher_len, nt = 0, nv = 0;
    int m = compressocrat_decrypt(ctx->cipher, C, st->key, g_cr_decode, CR_FILLER, &nt, &nv);
    if (m > MAX_CIPHER_LENGTH) m = MAX_CIPHER_LENGTH;    // clamp for the fixed-size result buffer

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < m; i++) plaintext_string[i] = index_to_char(g_cr_decode[i]);
    plaintext_string[m] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    char alpha[CR_ALPHA + 1]; compressocrat_alpha_string(st->key, alpha);

    printf("\nResult Score: %.2f | Words: %d | valid=%d/%d | alphabet=%s\n",
        score, n_words_found, nv, nt, alpha);

    print_cipher(ctx->cipher, C, NULL);
    printf("\n");
    print_text(g_cr_decode, m);
    printf("\n");
    print_solution_check(g_cr_decode, m);
    print_spaces_line(g_spaces_table, g_cr_decode, m);
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;              // no period for Compressocrat
        vec_copy(g_cr_decode, ctx->result->decrypted, m);
        ctx->result->decrypted_len = m;
    }

    // One-liner summary: >>> score, [words,] type, alphabet=, valid=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, alphabet=%s, valid=%d/%d, ",
            score, n_words_found, cfg->cipher_type, alpha, nv, nt);
    else
        printf(">>> %.2f, %d, alphabet=%s, valid=%d/%d, ",
            score, cfg->cipher_type, alpha, nv, nt);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(ctx->cipher, C, NULL);
    printf(", ");
    print_text(g_cr_decode, m);
    printf("\n");
}

static const CipherModel COMPRESSOCRAT_MODEL = {
    .name = "compressocrat", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = compressocrat_enumerate, .key_len = NULL,
    .seed = compressocrat_seed, .perturb = compressocrat_perturb, .copy_state = compressocrat_copy,
    .decrypt = compressocrat_decrypt_hook, .report = compressocrat_report,
};

// ===================================================================
//  Entry point
// ===================================================================

void solve_compressocrat(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Compressocrat needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 4) {
        printf("\n\nERROR: ciphertext too short for a Compressocrat solve.\n\n");
        return;
    }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= g_alpha) {
            printf("\n\nERROR: Compressocrat ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    if (cfg->verbose)
        printf("\ncompressocrat: %d ciphertext letters, keyed-alphabet anneal (kw %d..%d)\n",
            cipher_len, CR_KW_MIN, CR_KW_MAX);

    // Score over the fixed ciphertext length C; the decrypt hook tiles the variable decode to C.
    // Cribs are not used (the length change breaks the positional mapping).
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, 0);
    ctx.result = result;

    run_solver(&COMPRESSOCRAT_MODEL, &ctx);
}

#include "baconian_solver.h"
#include "baconian.h"
#include "engine.h"
#include "scoring.h"

// =====================================================================
//  Baconian solver (TYPE baconian)
// =====================================================================
//
// The plaintext is biliteral-encoded (baconian.c: each letter -> five a/b symbols via
// the FIXED 24-letter table, I=J/U=V) and then CONCEALED in normal English cover text.
// A hidden CLASSIFIER labels each cover unit a or b; two grouping modes appear in the
// ACA examples:
//   - PER-LETTER: every cover letter is one a/b symbol.
//   - PER-WORD:   the FIRST letter of each word is one a/b symbol.
// Both ACA examples use the classifier A-M='a', N-Z='b'; polarity is a second unknown.
//
// KEY = the classifier, a 26-letter a/b LABELLING (st->key[letter] in {BAC_A,BAC_B}).
// The 5-bit->letter decode is fixed, so this labelling is the ONLY unknown. The grouping
// mode is carried in the outer config (cc->aux[0]); polarity folds into the labelling.
//
// SEARCH. Per selected mode the model enumerates:
//   - CANONICAL SWEEP configs (key_len==0: a single seed+decode+score) for the standard
//     rules {A-M/N-Z, vowel/consonant} x both polarities -- so a canonical ACA Baconian
//     (both PDF examples) is decoded EXACTLY and instantly; plus
//   - a FREE-CLIMB config (key_len==26, random seed) whose neighbour move flips one
//     letter's a/b label -- SHAPE_ANNEAL, so all three -method schedules apply. This
//     generalises to non-canonical classifiers and gives the engine something to
//     calibrate (tests/test_baconian_solver.c).
//
// LENGTH HANDLING (the decode length varies per classifier; the engine scores a fixed
// length). We decode to nplain letters and CYCLICALLY TILE them to fill C = the cover-
// letter count so the mean n-gram is length-fair (the Fractionated-Morse pattern). A
// biliteral-VALIDITY reward (fraction of 5-bit groups that are legal codes, v<=23) is
// folded into score_adjust: the true classifier decodes ~100% valid, so this gives the
// anneal a decoupled gradient toward clean classifiers. Effectively needs -logprob.
// Cribs are NOT used (the length change + tiling break the positional mapping).
//
// LIMITATION. Blind classifier recovery at the ACA <=25-letter maximum is gaming-limited
// (n-gram is weak at that length; a fluent-but-wrong classifier can out-score the truth)
// -- the reliable case is a canonical classifier, caught exactly by the sweep configs.
// Phase offset is fixed at 0 (ACA cover text is exact; a trailing partial group drops).

// Grouping mode lives in cc->aux[0] (BAC_MODE_LETTER / BAC_MODE_WORD, from colossus.h).
// Seed selector lives in cc->aux[1]:
enum {
    BAC_SEED_AMNZ_A = 0,   // A-M -> a, N-Z -> b   (canonical, sweep)
    BAC_SEED_AMNZ_B = 1,   //   the opposite polarity
    BAC_SEED_VOWEL_A = 2,  // vowels -> a, consonants -> b
    BAC_SEED_VOWEL_B = 3,  //   the opposite polarity
    BAC_SEED_RANDOM = 4,   // random labelling -> the free anneal (climb)
    BAC_N_SEEDS = 5
};

// Parsed cover text, built ONCE in solve_baconian on the main thread and only READ during
// the search (thread-safe shared, per the CLAUDE.md rule). letters[] are the alphabetic
// cover letters (0..25) in order; word_first[i] marks the first alpha letter of each word.
static struct {
    int  letters[MAX_CIPHER_LENGTH];
    int  word_first[MAX_CIPHER_LENGTH];
    int  nletters;
    int  nwords;
} g_bacon;

// The tiled mean n-gram is PERIODIC (period = the decoded length), so scoring more than a few
// cycles is pure waste. We cap the scored/tiled length at BAC_SCORE_CAP -- this makes per-word
// grouping (whose cover text is ~5x longer than per-letter for the same plaintext) as cheap to
// score as per-letter, with no loss of signal (BAC_SCORE_CAP >> any realistic decoded length).
#define BAC_SCORE_CAP 500

// Per-thread scratch, off the stack.
static _Thread_local int g_bac_bits[MAX_CIPHER_LENGTH];       // the a/b stream (<= nletters)
static _Thread_local int g_bac_decode[MAX_CIPHER_LENGTH];     // the decoded plaintext (nbits/5)

static inline int bac_is_vowel(int l) {
    return (l == 0 || l == 4 || l == 8 || l == 14 || l == 20);   // A E I O U
}

// Extract the a/b bit stream for `mode` from the cover text under classifier key[].
// Returns the bit count.
static int bac_build_bits(int mode, const int key[], int out_bits[]) {
    int nb = 0;
    for (int i = 0; i < g_bacon.nletters; i++) {
        if (mode == BAC_MODE_WORD && !g_bacon.word_first[i]) continue;
        out_bits[nb++] = key[g_bacon.letters[i]] ? BAC_B : BAC_A;
    }
    return nb;
}

// ===================================================================
//  CipherModel hooks
// ===================================================================

static int bac_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    int mode = ctx->cfg->bacon_mode;
    int modes[2], nm = 0;
    if (mode == BAC_MODE_LETTER)      modes[nm++] = BAC_MODE_LETTER;
    else if (mode == BAC_MODE_WORD)   modes[nm++] = BAC_MODE_WORD;
    else { modes[nm++] = BAC_MODE_LETTER; modes[nm++] = BAC_MODE_WORD; }  // auto: sweep both

    int n = 0;
    for (int mi = 0; mi < nm; mi++)
        for (int s = 0; s < BAC_N_SEEDS && n < cap; s++) {
            out[n].period = 0; out[n].j = 0; out[n].k = 0;
            out[n].aux[0] = modes[mi]; out[n].aux[1] = s;
            n++;
        }
    return n;
}

// Canonical seeds are single SWEEP cells; the random seed is a full climb.
static int bac_key_len(const SolverCtx *ctx, const SolverConfig *cc) {
    (void) ctx;
    return (cc->aux[1] == BAC_SEED_RANDOM) ? ALPHABET_SIZE : 0;
}

static void bac_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    (void) ctx;
    for (int l = 0; l < ALPHABET_SIZE; l++) {
        int label;
        switch (cc->aux[1]) {
            case BAC_SEED_AMNZ_A:  label = (l >= 13) ? BAC_B : BAC_A; break;
            case BAC_SEED_AMNZ_B:  label = (l >= 13) ? BAC_A : BAC_B; break;
            case BAC_SEED_VOWEL_A: label = bac_is_vowel(l) ? BAC_A : BAC_B; break;
            case BAC_SEED_VOWEL_B: label = bac_is_vowel(l) ? BAC_B : BAC_A; break;
            default:               label = rand_int(0, 2); break;      // BAC_SEED_RANDOM
        }
        st->key[l] = label;
    }
    st->key_len = ALPHABET_SIZE;
}

static void bac_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                        SolverState *st, bool *force_primary) {
    (void) ctx; (void) cc; (void) force_primary;
    st->key[rand_int(0, ALPHABET_SIZE)] ^= 1;             // flip one letter's a/b label
    if (frand() < 0.15) st->key[rand_int(0, ALPHABET_SIZE)] ^= 1;   // occasional second flip
}

static void bac_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    for (int i = 0; i < ALPHABET_SIZE; i++) dst->key[i] = src->key[i];
    dst->key_len = src->key_len;
}

// Decrypt into out[0..C-1]: build the bit stream for this mode + classifier, decode to
// nplain letters (filler for an invalid group), then cyclically tile to the fixed length
// C. score_adjust = the biliteral-validity reward (fraction of valid 5-bit groups).
static void bac_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                             SolverState *st, int *out, double *score_adjust) {
    int C = ctx->cipher_len, nt = 0, nv = 0;
    int nb = bac_build_bits(cc->aux[0], st->key, g_bac_bits);
    int m = baconian_decode(g_bac_bits, nb, g_bac_decode, BAC_FILLER, &nt, &nv);
    if (m <= 0) {
        for (int i = 0; i < C; i++) out[i] = BAC_FILLER;
        *score_adjust = 0.0;
        return;
    }
    for (int i = 0; i < C; i++) out[i] = g_bac_decode[i % m];    // tile to the fixed length C
    *score_adjust = BAC_VALID_WEIGHT * ((nt > 0) ? (double) nv / nt : 0.0);
}

// ===================================================================
//  Reporting
// ===================================================================

static void bac_classifier_string(const int key[], char out[]) {
    for (int l = 0; l < ALPHABET_SIZE; l++) out[l] = key[l] ? 'b' : 'a';
    out[ALPHABET_SIZE] = '\0';
}

static void bac_report(const SolverCtx *ctx, const SolverConfig *cc,
                       const SolverState *st, double score, int *decrypted) {
    (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    int nt = 0, nv = 0;
    int nb = bac_build_bits(cc->aux[0], st->key, g_bac_bits);
    int m = baconian_decode(g_bac_bits, nb, g_bac_decode, BAC_FILLER, &nt, &nv);
    if (m > MAX_CIPHER_LENGTH) m = MAX_CIPHER_LENGTH;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < m; i++) plaintext_string[i] = index_to_char(g_bac_decode[i]);
    plaintext_string[m] = '\0';
    if (cfg->dictionary_present && ctx->shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, ctx->shared->dict,
            ctx->shared->n_dict_words, ctx->shared->max_dict_word_len);

    const char *mode = (cc->aux[0] == BAC_MODE_WORD) ? "per-word" : "per-letter";
    char cls[ALPHABET_SIZE + 1]; bac_classifier_string(st->key, cls);

    printf("\nResult Score: %.2f | Words: %d | valid=%d/%d | mode=%s | classifier(a/b)=%s\n",
        score, n_words_found, nv, nt, mode, cls);

    print_cipher(g_bacon.letters, g_bacon.nletters, NULL);
    printf("\n");
    print_text(g_bac_decode, m);
    printf("\n");
    print_solution_check(g_bac_decode, m);
    print_spaces_line(g_spaces_table, g_bac_decode, m);
    printf("%s\n", ctx->cribtext);

    if (ctx->result) {
        ctx->result->solved = true;
        ctx->result->cipher_type = cfg->cipher_type;
        ctx->result->score = score;
        ctx->result->n_words = n_words_found;
        ctx->result->cycleword_len = 0;              // no period for Baconian
        vec_copy(g_bac_decode, ctx->result->decrypted, m);
        ctx->result->decrypted_len = m;
    }

    // One-liner: >>> score, [words,] type, mode=, classifier=, valid=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, mode=%s, classifier=%s, valid=%d/%d, ",
            score, n_words_found, cfg->cipher_type, mode, cls, nv, nt);
    else
        printf(">>> %.2f, %d, mode=%s, classifier=%s, valid=%d/%d, ",
            score, cfg->cipher_type, mode, cls, nv, nt);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(g_bacon.letters, g_bacon.nletters, NULL);
    printf(", ");
    print_text(g_bac_decode, m);
    printf("\n");
}

static const CipherModel BACONIAN_MODEL = {
    .name = "baconian", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = bac_enumerate, .key_len = bac_key_len,
    .seed = bac_seed, .perturb = bac_perturb, .copy_state = bac_copy,
    .decrypt = bac_decrypt_hook, .report = bac_report,
};

// ===================================================================
//  Entry point
// ===================================================================

void solve_baconian(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_indices; (void) cipher_len;
    (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Baconian needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    // Parse the cover text (space-significant): the alphabetic cover letters + word starts.
    // A word starts at the first alpha letter after any run of non-alpha characters.
    g_bacon.nletters = 0; g_bacon.nwords = 0;
    int in_word = 0;
    for (const char *p = ciphertext_str; *p && g_bacon.nletters < MAX_CIPHER_LENGTH; p++) {
        unsigned char c = (unsigned char) *p;
        if (isalpha(c)) {
            int wf = !in_word;
            if (wf) g_bacon.nwords++;
            g_bacon.word_first[g_bacon.nletters] = wf;
            g_bacon.letters[g_bacon.nletters] = g_char_to_idx[toupper(c)];
            g_bacon.nletters++;
            in_word = 1;
        } else {
            in_word = 0;
        }
    }

    if (g_bacon.nletters < 5) {
        printf("\n\nERROR: ciphertext too short for a Baconian solve (need >= 5 cover letters).\n\n");
        return;
    }
    if (cfg->bacon_mode == BAC_MODE_WORD && g_bacon.nwords < 5) {
        printf("\n\nERROR: -baconmode word needs >= 5 words of cover text (got %d).\n\n", g_bacon.nwords);
        return;
    }

    if (cfg->verbose)
        printf("\nbaconian: %d cover letters, %d words; classifier anneal (mode=%s)\n",
            g_bacon.nletters, g_bacon.nwords,
            cfg->bacon_mode == BAC_MODE_LETTER ? "per-letter" :
            cfg->bacon_mode == BAC_MODE_WORD ? "per-word" : "auto (both)");

    // Score over a fixed length C (the cover-letter count, capped at BAC_SCORE_CAP); the decrypt
    // hook tiles the variable decode to C. The full cover text stays in g_bacon (the hooks read
    // it directly). Cribs are not used (the length change breaks the positional mapping).
    int C = (g_bacon.nletters < BAC_SCORE_CAP) ? g_bacon.nletters : BAC_SCORE_CAP;
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        g_bacon.letters, C, crib_indices, crib_positions, 0);
    ctx.result = result;

    run_solver(&BACONIAN_MODEL, &ctx);
}

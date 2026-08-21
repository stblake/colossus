#include "affine_solver.h"
#include "affine.h"
#include "scoring.h"
#include <time.h>

// =====================================================================
//  Affine solver -- deterministic exhaustive 12x26 = 312-key search. See affine_solver.h.
// =====================================================================

// Letters-only search stream + the original layout (spaces/punct sentinels), off the stack.
static int g_affine_letters[MAX_CIPHER_LENGTH];
static int g_affine_orig[MAX_CIPHER_LENGTH];

// Aristocrat-style spaced print: walk the original layout, emitting each sentinel verbatim and
// each letter as its recovered value (restoring the ciphertext's word divisions).
static void affine_print_spaced(const int orig[], int orig_len, const int letters[]) {
    int c = 0;
    for (int i = 0; i < orig_len; i++)
        putchar(orig[i] >= 0 ? index_to_char(letters[c++]) : index_to_char(orig[i]));
}

void solve_affine(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str; (void) crib_indices; (void) crib_positions; (void) n_cribs;

    if (g_alpha != DEFAULT_ALPHABET_SIZE) {
        printf("\n\nERROR: Affine needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    // Derive the letters-only stream (spaces/punctuation dropped from the search, restored in
    // the report); every position is a letter so scoring is a plain n-gram walk.
    int nl = 0;
    for (int i = 0; i < cipher_len; i++) {
        g_affine_orig[i] = cipher_indices[i];
        if (cipher_indices[i] >= 0) g_affine_letters[nl++] = cipher_indices[i];
    }
    if (nl < 4) {
        printf("\n\nERROR: ciphertext too short for an Affine solve (%d letters).\n\n", nl);
        return;
    }

    if (cfg->verbose)
        printf("\naffine: %d letters, deterministic exhaustive 12x26 = 312-key search\n", nl);

    static int dec_stream[MAX_CIPHER_LENGTH];
    static int best_pt[MAX_CIPHER_LENGTH];
    int best_a = 1, best_b = 0, have = 0;
    double best = -1e300;

    clock_t t0 = clock();
    for (int mi = 0; mi < AFFINE_N_MULT; mi++) {
        int ainv = affine_inverses[mi];
        for (int b = 0; b < AFFINE_MOD; b++) {
            int dec[AFFINE_MOD];
            for (int c = 0; c < AFFINE_MOD; c++)
                dec[c] = (ainv * ((c - b + AFFINE_MOD) % AFFINE_MOD)) % AFFINE_MOD;
            for (int i = 0; i < nl; i++) dec_stream[i] = dec[g_affine_letters[i]];
            double s = state_score(dec_stream, nl, NULL, NULL, 0,
                                   shared->ngram_data, cfg->ngram_size,
                                   cfg->weight_ngram, 0.f, 0.f, 0.f);
            if (!have || s > best) {
                best = s; have = 1;
                best_a = affine_multipliers[mi]; best_b = b;
                for (int i = 0; i < nl; i++) best_pt[i] = dec_stream[i];
            }
        }
    }
    double elapsed = ((double) clock() - t0) / CLOCKS_PER_SEC;

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < nl; i++) plaintext_string[i] = index_to_char(best_pt[i]);
    plaintext_string[nl] = '\0';
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);

    printf("\nResult Score: %.2f | Words: %d | a=%d b=%d | %.3f sec\n",
        best, n_words_found, best_a, best_b, elapsed);
    affine_print_spaced(g_affine_orig, cipher_len, g_affine_letters); printf("\n");
    affine_print_spaced(g_affine_orig, cipher_len, best_pt);          printf("\n");
    printf("%s\n", cribtext_str ? cribtext_str : "");

    if (result) {
        result->solved = true;
        result->cipher_type = cfg->cipher_type;
        result->score = best;
        result->n_words = n_words_found;
        result->cycleword_len = 0;
        vec_copy(best_pt, result->decrypted, nl);
        result->decrypted_len = nl;
    }

    // >>> score, [words,] type, a=, b=, file, CIPHER, PLAINTEXT (letters-only, whitespace-strippable).
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, a=%d, b=%d, ", best, n_words_found, cfg->cipher_type, best_a, best_b);
    else
        printf(">>> %.2f, %d, a=%d, b=%d, ", best, cfg->cipher_type, best_a, best_b);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_text(g_affine_letters, nl);
    printf(", ");
    print_text(best_pt, nl);
    printf("\n");
}

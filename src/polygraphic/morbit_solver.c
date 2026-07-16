#include "morbit_solver.h"
#include "morbit.h"
#include "scoring.h"
#include <time.h>

// =====================================================================
//  Morbit solver -- deterministic exhaustive 9! search. See morbit_solver.h.
// =====================================================================

// Parse a stream of single decimal digits 1..9 (any other character, including 0, is
// a separator) into out[]; returns the count. A Morbit ciphertext is a run of digits
// 1..9, one per Morse-symbol pair (there is no 0 -- the key numbers 1..9 come from a
// 9-letter keyword's alphabetical ranks).
static int mb_parse_digits(const char *s, int out[], int cap) {
    int n = 0;
    for (int i = 0; s[i] && n < cap; i++)
        if (s[i] >= '1' && s[i] <= '9') out[n++] = s[i] - '0';
    return n;
}

double morbit_search(const int *digits, int clen,
    float *ngram_data, int ngram_size, float weight_ngram, double valid_weight,
    int best_key[10], int best_pt[], int *best_n, int *best_nt, int *best_nv) {

    static int pt[MAX_CIPHER_LENGTH];
    int perm[9], key[10], c[9] = {0};
    double best = -1e300;
    int have = 0;
    *best_n = 0; *best_nt = 0; *best_nv = 0;

    key[0] = 0;                              // digit 0 never occurs in a Morbit cipher
    for (int p = 0; p < 9; p++) perm[p] = p; // perm[d-1] = pair assigned to digit d

    // Heap's algorithm: enumerate all 9! = 362,880 digit(1..9) -> pair(0..8)
    // bijections. The first pass scores the identity perm; each later pass advances by
    // a single swap and scores again, for 9! evaluations total.
    int i = 0, started = 0;
    for (;;) {
        if (started) {
            while (i < 9 && c[i] >= i) { c[i] = 0; i++; }
            if (i >= 9) break;
            int a = (i & 1) ? c[i] : 0;      // swap index per Heap's algorithm
            int t = perm[a]; perm[a] = perm[i]; perm[i] = t;
            c[i]++;
            i = 0;
        }
        started = 1;

        for (int d = 1; d <= 9; d++) key[d] = perm[d - 1];
        int nt = 0, nv = 0;
        int n = morbit_decrypt(digits, clen, key, pt, MORBIT_FILLER, &nt, &nv);
        double s = state_score(pt, n, NULL, NULL, 0, ngram_data, ngram_size,
                               weight_ngram, 0.f, 0.f, 0.f);
        if (nt > 0) s += valid_weight * ((double) nv / (double) nt);

        if (!have || s > best) {
            best = s; have = 1;
            for (int d = 0; d <= 9; d++) best_key[d] = key[d];
            for (int k = 0; k < n; k++) best_pt[k] = pt[k];
            *best_n = n; *best_nt = nt; *best_nv = nv;
        }
    }
    return best;
}

void solve_morbit(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_indices; (void) cipher_len;            // digit input is parsed here, not decoded
    (void) crib_indices; (void) crib_positions; (void) n_cribs;   // cribs unused (length change)

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Morbit needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    static int digits[MAX_CIPHER_LENGTH];
    int clen = mb_parse_digits(ciphertext_str, digits, MAX_CIPHER_LENGTH);
    if (clen < 8) {
        printf("\n\nERROR: parsed only %d ciphertext digits; need >= 8. A Morbit "
               "ciphertext is a stream of digits 1-9.\n\n", clen);
        return;
    }

    if (cfg->verbose)
        printf("\nmorbit: %d ciphertext digits, exhaustive 9! = 362880-key search "
               "(validity weight %.2f)\n", clen, MORBIT_VALID_WEIGHT);

    static int best_pt[MAX_CIPHER_LENGTH];
    int best_key[10], n = 0, nt = 0, nv = 0;
    clock_t t0 = clock();
    double score = morbit_search(digits, clen, shared->ngram_data, cfg->ngram_size,
        cfg->weight_ngram, MORBIT_VALID_WEIGHT, best_key, best_pt, &n, &nt, &nv);
    double elapsed = ((double) clock() - t0) / CLOCKS_PER_SEC;

    // Recovered key as a 9-char pair->digit string (index = pair 0..8, value = digit
    // 1..9) -- the ACA array's "number" row, e.g. "958427136" for keyword WISECRACK.
    char keystr[10];
    for (int d = 1; d <= 9; d++) keystr[best_key[d]] = (char) ('0' + d);
    keystr[9] = '\0';

    int n_words_found = 0;
    char plaintext_string[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < n; i++) plaintext_string[i] = index_to_char(best_pt[i]);
    plaintext_string[n] = '\0';
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);

    // Digit cipher as a printable string (for the human block and the >>> CSV).
    static char cipherstr[MAX_CIPHER_LENGTH + 1];
    for (int i = 0; i < clen; i++) cipherstr[i] = (char) ('0' + digits[i]);
    cipherstr[clen] = '\0';

    printf("\nResult Score: %.2f | Words: %d | valid=%d/%d | key=%s | %.2f sec\n",
        score, n_words_found, nv, nt, keystr, elapsed);
    printf("%s\n", cipherstr);
    print_text(best_pt, n);
    printf("\n");
    print_solution_check(best_pt, n);
    print_spaces_line(g_spaces_table, best_pt, n);
    printf("%s\n", cribtext_str ? cribtext_str : "");

    if (result) {
        result->solved = true;
        result->cipher_type = cfg->cipher_type;
        result->score = score;
        result->n_words = n_words_found;
        result->cycleword_len = 0;                       // no period for Morbit
        vec_copy(best_pt, result->decrypted, n);
        result->decrypted_len = n;
    }

    // One-liner summary: >>> score, [words,] type, key=, valid=, file, CIPHER, PLAINTEXT
    if (cfg->dictionary_present)
        printf(">>> %.2f, %d, %d, key=%s, valid=%d/%d, ",
            score, n_words_found, cfg->cipher_type, keystr, nv, nt);
    else
        printf(">>> %.2f, %d, key=%s, valid=%d/%d, ",
            score, cfg->cipher_type, keystr, nv, nt);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    printf("%s, ", cipherstr);
    print_text(best_pt, n);
    printf("\n");
}

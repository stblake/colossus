#include "pollux_solver.h"
#include "pollux.h"
#include "scoring.h"
#include <time.h>

// =====================================================================
//  Pollux solver -- deterministic exhaustive 3^10 search. See pollux_solver.h.
// =====================================================================

// Parse a stream of single decimal digits (any non-digit is a separator) into
// out[]; returns the count. Unlike Nihilist Substitution's number parser, each digit
// is one symbol (a Pollux ciphertext is a run of digits 0..9, one per Morse symbol).
static int px_parse_digits(const char *s, int out[], int cap) {
    int n = 0;
    for (int i = 0; s[i] && n < cap; i++)
        if (s[i] >= '0' && s[i] <= '9') out[n++] = s[i] - '0';
    return n;
}

double pollux_search(const int *digits, int clen,
    float *ngram_data, int ngram_size, float weight_ngram, double valid_weight,
    int best_key[10], int best_pt[], int *best_n, int *best_nt, int *best_nv) {

    static int pt[MAX_CIPHER_LENGTH];
    int key[10];
    double best = -1e300;
    int have = 0;
    *best_n = 0; *best_nt = 0; *best_nv = 0;

    // Every map digit -> {DOT,DASH,X} is a base-3 number code in [0, 3^10). Skip keys
    // that cannot segment (no X) or cannot form both dot- and dash-letters (the true
    // key always has all three), a cheap prune that also avoids degenerate optima.
    const int total = 59049;                     // 3^10
    for (int code = 0; code < total; code++) {
        int c = code, has_x = 0, has_dot = 0, has_dash = 0;
        for (int d = 0; d < 10; d++) {
            int e = c % 3; c /= 3;
            key[d] = e;
            if      (e == PX_X)   has_x = 1;
            else if (e == PX_DOT) has_dot = 1;
            else                  has_dash = 1;
        }
        if (!has_x || !has_dot || !has_dash) continue;

        int nt = 0, nv = 0;
        int n = pollux_decrypt(digits, clen, key, pt, POLLUX_FILLER, &nt, &nv);
        double s = state_score(pt, n, NULL, NULL, 0, ngram_data, ngram_size,
                               weight_ngram, 0.f, 0.f, 0.f);
        if (nt > 0) s += valid_weight * ((double) nv / (double) nt);

        if (!have || s > best) {
            best = s; have = 1;
            for (int d = 0; d < 10; d++) best_key[d] = key[d];
            for (int i = 0; i < n; i++) best_pt[i] = pt[i];
            *best_n = n; *best_nt = nt; *best_nv = nv;
        }
    }
    return best;
}

void solve_pollux(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) cipher_indices; (void) cipher_len;            // digit input is parsed here, not decoded
    (void) crib_indices; (void) crib_positions; (void) n_cribs;   // cribs unused (length change)

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Pollux needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }

    static int digits[MAX_CIPHER_LENGTH];
    int clen = px_parse_digits(ciphertext_str, digits, MAX_CIPHER_LENGTH);
    if (clen < 8) {
        printf("\n\nERROR: parsed only %d ciphertext digits; need >= 8. A Pollux "
               "ciphertext is a stream of digits 0-9.\n\n", clen);
        return;
    }

    if (cfg->verbose)
        printf("\npollux: %d ciphertext digits, exhaustive 3^10 = 59049-key search "
               "(validity weight %.2f)\n", clen, POLLUX_VALID_WEIGHT);

    static int best_pt[MAX_CIPHER_LENGTH];
    int best_key[10], n = 0, nt = 0, nv = 0;
    clock_t t0 = clock();
    double score = pollux_search(digits, clen, shared->ngram_data, cfg->ngram_size,
        cfg->weight_ngram, POLLUX_VALID_WEIGHT, best_key, best_pt, &n, &nt, &nv);
    double elapsed = ((double) clock() - t0) / CLOCKS_PER_SEC;

    // Recovered key as a 10-char .-x string indexed by digit 0..9.
    char keystr[11];
    for (int d = 0; d < 10; d++)
        keystr[d] = best_key[d] == PX_DOT ? '.' : (best_key[d] == PX_DASH ? '-' : 'x');
    keystr[10] = '\0';

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
    printf("\n%s\n", cribtext_str ? cribtext_str : "");

    if (result) {
        result->solved = true;
        result->cipher_type = cfg->cipher_type;
        result->score = score;
        result->n_words = n_words_found;
        result->cycleword_len = 0;                       // no period for Pollux
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

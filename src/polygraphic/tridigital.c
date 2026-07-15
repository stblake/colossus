//
//  Tridigital cipher primitives (see tridigital.h for the full convention). A 3x10 keyed
//  block over the full 26-letter alphabet enciphers plaintext into a DIGIT stream (0..9):
//  each letter -> the digit heading its column, a word break -> the separator digit (the
//  digit above the blank 10th column). Because a letter's row does not affect its digit the
//  map is 3-to-1 (deliberately ambiguous). Everything is derived once by tridigital_build_grid();
//  encrypt shares it with the solver, the generator and the tests, so the convention can never drift.
//

#include "tridigital.h"
#include <ctype.h>

// Validate that p[0..n-1] is a permutation of 0..n-1 (n <= TRI_NALPHA).
static int tri_is_perm(const int p[], int n) {
    int seen[TRI_NALPHA];
    if (n > TRI_NALPHA) return 0;
    for (int i = 0; i < n; i++) seen[i] = 0;
    for (int i = 0; i < n; i++) {
        if (p[i] < 0 || p[i] >= n || seen[p[i]]) return 0;
        seen[p[i]] = 1;
    }
    return 1;
}

// The letter-column of reading-order cell k: row0 -> cols 0..8, row1 -> cols 0..8,
// row2 -> cols 0..7 (the 26-letter row-major fill of the 9 letter-columns).
static int tri_col_of_cell(int k) {
    return (k < 9) ? k : (k < 18) ? (k - 9) : (k - 18);
}

int tridigital_build_grid(TridigitalGrid *g, const int keyed_alpha[TRI_NALPHA],
                          const int col_digit[TRI_COLS]) {
    if (!tri_is_perm(keyed_alpha, TRI_NALPHA)) return -1;
    // col_digit must be a permutation of 0..9.
    int seen[10]; for (int d = 0; d < 10; d++) seen[d] = 0;
    for (int c = 0; c < TRI_COLS; c++) {
        int d = col_digit[c];
        if (d < 0 || d > 9 || seen[d]) return -1;
        seen[d] = 1;
    }

    for (int c = 0; c < TRI_COLS; c++) { g->digit_of_col[c] = col_digit[c]; g->col_of_digit[col_digit[c]] = c; }
    g->sep_digit = col_digit[TRI_LETTER_COLS];   // the 10th physical column (index 9) is the separator

    for (int c = 0; c < TRI_LETTER_COLS; c++) g->ncol[c] = 0;
    for (int s = 0; s < TRI_NALPHA; s++) { g->col_of_letter[s] = -1; g->digit_of_letter[s] = -1; }

    for (int k = 0; k < TRI_NALPHA; k++) {
        int col = tri_col_of_cell(k);
        int L = keyed_alpha[k];
        g->col_of_letter[L] = col;
        g->digit_of_letter[L] = col_digit[col];
        g->letters_of_col[col][g->ncol[col]++] = L;
    }
    return 0;
}

int tridigital_build_from_keywords(TridigitalGrid *g, const char *kw_cols, const char *kw_alpha) {
    // Column digits: 10 DISTINCT keyword letters -> digit by alphabetical rank (1..10, 10 -> 0).
    int cols[TRI_COLS], m = 0, seen[26] = {0};
    for (const char *p = kw_cols; *p && m < TRI_COLS; p++) {
        int c = toupper((unsigned char) *p);
        if (c < 'A' || c > 'Z') continue;
        int L = c - 'A';
        if (seen[L]) continue;                 // must be distinct letters
        seen[L] = 1; cols[m++] = L;
    }
    if (m != TRI_COLS) return -1;
    int col_digit[TRI_COLS];
    for (int i = 0; i < TRI_COLS; i++) {
        int rank = 1;                          // 1-based alphabetical rank (ties by position)
        for (int j = 0; j < TRI_COLS; j++)
            if (cols[j] < cols[i] || (cols[j] == cols[i] && j < i)) rank++;
        col_digit[i] = rank % 10;              // rank 10 -> digit 0
    }

    // Keyed 26-letter alphabet: keyword deduped, then the remaining letters A..Z in order.
    int keyed[TRI_NALPHA], used[26] = {0}, k = 0;
    for (const char *p = kw_alpha; *p; p++) {
        int c = toupper((unsigned char) *p);
        if (c < 'A' || c > 'Z') continue;
        int L = c - 'A';
        if (!used[L]) { used[L] = 1; keyed[k++] = L; }
    }
    for (int L = 0; L < 26; L++) if (!used[L]) keyed[k++] = L;

    return tridigital_build_grid(g, keyed, col_digit);
}

int tridigital_encrypt(const int plain[], int n, const TridigitalGrid *g, int out[]) {
    for (int i = 0; i < n; i++) {
        int x = plain[i];
        out[i] = (x >= 0 && x < TRI_NALPHA) ? g->digit_of_letter[x] : g->sep_digit;
    }
    return n;
}

int tridigital_candidates(int digit, const TridigitalGrid *g, int out_letters[TRI_MAXCOLLET]) {
    if (digit < 0 || digit > 9 || digit == g->sep_digit) return 0;
    int col = g->col_of_digit[digit];
    int nc = g->ncol[col];
    for (int i = 0; i < nc; i++) out_letters[i] = g->letters_of_col[col][i];
    return nc;
}

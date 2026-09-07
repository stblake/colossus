#include "grandpre.h"

// =====================================================================
//  Grandpre cipher primitive -- see grandpre.h
// =====================================================================

// Fill the reverse index (letter -> its cells' codes) from g->cell / g->N. Returns 0 if all
// 26 letters occur at least once, -1 otherwise.
static int grandpre_index_cells(GrandpreSquare *g) {
    int N = g->N;
    for (int c = 0; c < ALPHABET_SIZE; c++) g->ncodes[c] = 0;
    for (int r = 0; r < N; r++)
        for (int c = 0; c < N; c++) {
            int L = g->cell[r * N + c];
            if (L < 0 || L >= ALPHABET_SIZE) return -1;
            g->codes[L][g->ncodes[L]++] = grandpre_code(r, c, N);
        }
    for (int c = 0; c < ALPHABET_SIZE; c++) if (g->ncodes[c] == 0) return -1;
    return 0;
}

int grandpre_build(GrandpreSquare *g, int N, const int letters[]) {
    if (N < GRANDPRE_MIN_SIDE || N > GRANDPRE_MAX_SIDE) return -1;
    g->N = N;
    for (int i = 0; i < N * N; i++) {
        if (letters[i] < 0 || letters[i] >= ALPHABET_SIZE) return -1;
        g->cell[i] = letters[i];
    }
    return grandpre_index_cells(g);
}

int grandpre_build_from_words(GrandpreSquare *g, int N, const char *const words[]) {
    if (N < GRANDPRE_MIN_SIDE || N > GRANDPRE_MAX_SIDE) return -1;
    int letters[GRANDPRE_MAX_GRID];
    for (int r = 0; r < N; r++) {
        const char *w = words[r];
        int c = 0;
        for (; c < N; c++) {
            char ch = w[c];
            if (ch >= 'A' && ch <= 'Z') letters[r * N + c] = ch - 'A';
            else if (ch >= 'a' && ch <= 'z') letters[r * N + c] = ch - 'a';
            else return -1;                       // too short / non-letter
        }
        if (w[c] != '\0') return -1;              // longer than N letters
    }
    return grandpre_build(g, N, letters);
}

int grandpre_encrypt(const int plain[], int n, const GrandpreSquare *g, int out_codes[]) {
    int m = 0;
    for (int i = 0; i < n; i++) {
        int L = plain[i];
        if (L < 0 || L >= ALPHABET_SIZE) continue;            // skip non-letters (no separators)
        int k = rand_int(0, g->ncodes[L]);                   // uniform homophone choice (isolog)
        out_codes[m++] = g->codes[L][k];
    }
    return m;
}

int grandpre_decrypt(const int codes[], int ncodes, const GrandpreSquare *g, int out_pt[]) {
    int N = g->N, m = 0;
    for (int i = 0; i < ncodes; i++) {
        int code = codes[i];
        int r = grandpre_index(code / 10, N), c = grandpre_index(code % 10, N);
        if (r < 0 || c < 0) continue;                        // label outside the square
        out_pt[m++] = g->cell[r * N + c];
    }
    return m;
}

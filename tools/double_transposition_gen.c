// Double columnar transposition generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), UPPERCASES letters and carries any
// spaces / punctuation as reversible sentinels (like colossus's ord()), and enciphers it
// with a DOUBLE columnar transposition under two keywords: C = colenc(colenc(P, K1), K2).
// The non-letter cells ride through both transposition stages as real grid cells, so the
// solver -- which scores letters only and reverses the sentinels on output -- recovers the
// spaces/punctuation in place. (A letters-only plaintext is bit-identical to before.)
// Each stage is a standard columnar transposition in colossus's convention (plaintext
// row-major, leftmost `len % K` columns one cell taller, columns read in key order,
// top-to-bottom). The keyword's letters (ranked alphabetically, ties left-to-right)
// give the column read order -- the same numeric-key convention as the paper.
//
//   make double_transposition_gen
//   ./double_transposition_gen plaintext.txt KEYWORD SECRET >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <keyword1> <keyword2> [tb|bt]
// stdout: the ciphertext (one line, bare A..Z).
// stderr: the plaintext solution (bare A..Z -- exactly what the solver recovers).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

// Column read order from a keyword: order[j] = grid column read j-th (the column with
// the j-th smallest keyword letter; equal letters keep left-to-right order).
static int keyword_order(const char *kw, int order[]) {
    int lets[MAX_COLS], K = 0;
    for (const char *k = kw; *k && K < MAX_COLS; k++) {
        int c = toupper((unsigned char) *k);
        if (c >= 'A' && c <= 'Z') lets[K++] = c - 'A';
    }
    // stable selection sort of column indices by letter
    int used[MAX_COLS];
    for (int i = 0; i < K; i++) used[i] = 0;
    for (int j = 0; j < K; j++) {
        int best = -1;
        for (int c = 0; c < K; c++) {
            if (used[c]) continue;
            if (best < 0 || lets[c] < lets[best]) best = c;
        }
        used[best] = 1;
        order[j] = best;
    }
    return K;
}

// One columnar transposition stage (encryption): write plaintext row-major into a grid
// of K columns (leftmost len%K columns one cell taller), read columns in `order` (TB).
static void columnar_encrypt(const int in[], int len, int K, const int order[], int out[]) {
    static int grid[MAXLEN];
    int R = (len + K - 1) / K, rem = len % K;
    int pos = 0;
    for (int r = 0; r < R; r++)
        for (int c = 0; c < K; c++) {
            int h = (rem == 0 || c < rem) ? R : R - 1;
            if (r < h) grid[r * K + c] = in[pos++];
        }
    int o = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j];
        int h = (rem == 0 || c < rem) ? R : R - 1;
        for (int r = 0; r < h; r++) out[o++] = grid[r * K + c];
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <plaintext|-> <keyword1> <keyword2> [tb|bt]\n", argv[0]);
        return 1;
    }
    init_alphabet(NULL);

    int order1[MAX_COLS], order2[MAX_COLS];
    int K1 = keyword_order(argv[2], order1);
    int K2 = keyword_order(argv[3], order2);
    if (K1 < 2 || K2 < 2) { fprintf(stderr, "keywords must have >= 2 distinct-position letters\n"); return 1; }

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        if (ch == '\r') continue;
        int c = toupper(ch);
        // Letters -> 0..25; spaces / punctuation -> reversible negative sentinel over the
        // ORIGINAL byte (matches ord() / index_to_char), carried as real grid cells.
        int v = (c >= 'A' && c <= 'Z') ? (c - 'A') : (-(int) (unsigned char) ch - 1);
        if (n < MAXLEN) raw[n++] = v;
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int mid[MAXLEN], cipher[MAXLEN];
    columnar_encrypt(raw, n, K1, order1, mid);
    columnar_encrypt(mid, n, K2, order2, cipher);

    for (int i = 0; i < n; i++) putchar(index_to_char(cipher[i]));
    putchar('\n');
    fprintf(stderr, "[double-transposition: K1=%s (%d cols), K2=%s (%d cols), %d letters]\n",
        argv[2], K1, argv[3], K2, n);
    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

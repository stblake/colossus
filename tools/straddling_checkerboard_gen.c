// Straddling Checkerboard cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z and the digits 0..9
// (digits exercise the figure-shift), and enciphers with the REAL straddling_encrypt
// (straddling_checkerboard.c + utils.c) so the generator and the solver can never drift.
//
//   make straddling_checkerboard_gen
//   ./straddling_checkerboard_gen plaintext.txt SECRET - 2 6 >cipher.txt 2>solution.txt
//   ./straddling_checkerboard_gen plaintext.txt SECRET 7290518364 0 4 >c.txt 2>s.txt
//
// argv: <plaintext|-> <arrangement-keyword> <label-key: 10 digits or -> <blankcol0> <blankcol1>
//   arrangement-keyword : letters; the 26-letter keyed alphabet (keyword then A..Z ascending)
//                         fills the board, with the FIG marker appended as the 27th symbol.
//   label-key           : a permutation of 0..9 (the column headings), or '-' for 0..9 in order.
//   blankcol0/1         : the two physical columns (0..9) left blank on the top row; their
//                         labels become the two row indicators.
// stdout: the ciphertext (one line of digits 0-9; length differs from N)
// stderr: the cleaned plaintext (the solution: A..Z + 0..9, what the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "straddling_checkerboard.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <plaintext|-> <arrangement-keyword> <label-key|-> "
                        "<blankcol0> <blankcol1>\n", argv[0]);
        return 1;
    }

    init_alphabet_adfgvx();                          // 36-symbol A..Z + 0..9 (matches the solver)
    if (g_alpha != 36) { fprintf(stderr, "alphabet is %d symbols, need 36\n", g_alpha); return 1; }

    // Board arrangement: 26-letter keyed alphabet (keyword deduped, then A..Z ascending), then FIG.
    int seq[SC_NSYM], used[26] = {0}, m = 0;
    for (const char *k = argv[2]; *k; k++) {
        int c = toupper((unsigned char) *k);
        if (c >= 'A' && c <= 'Z' && !used[c - 'A']) { used[c - 'A'] = 1; seq[m++] = c - 'A'; }
    }
    for (int c = 0; c < 26; c++) if (!used[c]) seq[m++] = c;
    seq[26] = SC_FIG;                                // the figure-shift marker (27th board symbol)

    // Column labels: a permutation of 0..9, or '-' for the identity 0..9.
    int labels[10];
    if (strcmp(argv[3], "-") == 0) {
        for (int c = 0; c < 10; c++) labels[c] = c;
    } else {
        if ((int) strlen(argv[3]) != 10) { fprintf(stderr, "label-key must be 10 digits or -\n"); return 1; }
        int seen[10] = {0};
        for (int c = 0; c < 10; c++) {
            int d = argv[3][c] - '0';
            if (d < 0 || d > 9 || seen[d]) { fprintf(stderr, "label-key must be a permutation of 0..9\n"); return 1; }
            seen[d] = 1; labels[c] = d;
        }
    }

    int b0 = atoi(argv[4]), b1 = atoi(argv[5]);
    if (b0 < 0 || b0 > 9 || b1 < 0 || b1 > 9 || b0 == b1) {
        fprintf(stderr, "blank columns must be two distinct values in 0..9\n"); return 1;
    }

    StraddlingBoard board;
    if (straddling_build_board(&board, seq, labels, labels[b0], labels[b1]) != 0) {
        fprintf(stderr, "board build failed\n"); return 1;
    }

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int c = toupper(ch);
        int idx = (c >= 0 && c < 128) ? g_char_to_idx[c] : -1;   // A..Z -> 0..25, 0..9 -> 26..35
        if (idx >= 0 && idx < 36 && n < MAXLEN) raw[n++] = idx;
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[3 * MAXLEN];
    int clen = straddling_encrypt(raw, n, &board, cipher);
    if (clen < 0) { fprintf(stderr, "encrypt failed (board cannot encipher a plaintext digit)\n"); return 1; }

    for (int i = 0; i < clen; i++) putchar('0' + cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

// Monome-Dinome cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps letters only (J->I, Z->Y so the
// 24-letter box can represent them), and enciphers with the REAL monome_dinome_encrypt
// (monome_dinome.c + utils.c) so the generator and the solver can never drift.
//
//   make monome_dinome_gen
//   ./monome_dinome_gen plaintext.txt SECRET 18927054 6 3 >cipher.txt 2>solution.txt
//   ./monome_dinome_gen plaintext.txt GROUCHOMARX - 3 7   >c.txt 2>s.txt
//
// argv: <plaintext|-> <box-keyword> <col-label-key: 8 digits or -> <indicator0> <indicator1>
//   box-keyword   : letters; the 24-letter keyed alphabet (keyword deduped, then the unused
//                   letters A..Y skipping J and Z, in order) fills the 3x8 box in reading order.
//   col-label-key : the 8 COLUMN-LABEL digits (a permutation of the ten digits minus the two
//                   indicators), or '-' for the eight non-indicator digits in ascending order.
//   indicator0/1  : the two ROW-INDICATOR digits (0..9, distinct, heading rows 2 and 3).
// stdout: the ciphertext (one line of digits 0-9; length differs from N)
// stderr: the cleaned plaintext (the solution: 24-letter set, what the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "monome_dinome.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <plaintext|-> <box-keyword> <col-label-key|-> "
                        "<indicator0> <indicator1>\n", argv[0]);
        return 1;
    }

    init_alphabet_monome_dinome();                   // 24-letter A..Z with J->I, Z->Y (matches solver)
    if (g_alpha != MD_NALPHA) { fprintf(stderr, "alphabet is %d letters, need %d\n", g_alpha, MD_NALPHA); return 1; }

    int r0 = atoi(argv[4]), r1 = atoi(argv[5]);
    if (r0 < 0 || r0 > 9 || r1 < 0 || r1 > 9 || r0 == r1) {
        fprintf(stderr, "indicators must be two distinct digits 0..9\n"); return 1;
    }

    // Box arrangement: 24-letter keyed alphabet (keyword deduped, then the unused 24-set letters).
    int letters[MD_NALPHA], used[MD_NALPHA] = {0}, m = 0;
    for (const char *k = argv[2]; *k; k++) {
        int c = toupper((unsigned char) *k);
        int idx = (c >= 0 && c < 128) ? g_char_to_idx[c] : -1;     // J->I, Z->Y already folded
        if (idx >= 0 && idx < MD_NALPHA && !used[idx]) { used[idx] = 1; letters[m++] = idx; }
    }
    for (int idx = 0; idx < MD_NALPHA; idx++) if (!used[idx]) letters[m++] = idx;

    // Column labels: the 8 non-indicator digits, keyed or ascending.
    int col_label[MD_COLS];
    if (strcmp(argv[3], "-") == 0) {
        int c = 0;
        for (int d = 0; d < 10; d++) if (d != r0 && d != r1) col_label[c++] = d;
    } else {
        if ((int) strlen(argv[3]) != MD_COLS) { fprintf(stderr, "col-label-key must be 8 digits or -\n"); return 1; }
        for (int c = 0; c < MD_COLS; c++) col_label[c] = argv[3][c] - '0';
    }

    MonomeDinomeBoard board;
    if (monome_dinome_build_board(&board, letters, col_label, r0, r1) != 0) {
        fprintf(stderr, "board build failed (col labels must be the 8 digits other than the indicators)\n");
        return 1;
    }

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int c = toupper(ch);
        int idx = (c >= 0 && c < 128) ? g_char_to_idx[c] : -1;      // A..Y -> 0..23 (J->I, Z->Y)
        if (idx >= 0 && idx < MD_NALPHA && n < MAXLEN) raw[n++] = idx;
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[2 * MAXLEN];
    int clen = monome_dinome_encrypt(raw, n, &board, cipher);
    if (clen < 0) { fprintf(stderr, "encrypt failed\n"); return 1; }

    for (int i = 0; i < clen; i++) putchar('0' + cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

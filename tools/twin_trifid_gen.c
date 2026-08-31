// Twin Trifid cipher generator (test-data tool, not part of the solver).
//
// A Twin Trifid is TWO Trifid messages enciphered under the SAME keyed 3x3x3 cube (27
// symbols, A..Z + '+') but with DIFFERENT periods; the two plaintexts share a common phrase.
// This tool reads two plaintext files, keeps A..Z, builds one cube from a keyword, and
// enciphers message 1 at period1 and message 2 at period2. It links the REAL cipher code
// (twin_trifid.c + trifid.c + utils.c) so the generator and solver can never drift.
//
//   make twin_trifid_gen
//   ./twin_trifid_gen pt1.txt pt2.txt KEYWORD 7 8 >ciphers.txt 2>twin.solution
//
// argv: <pt1-file|-> <pt2-file> <keyword> <period1> <period2>
// stdout: TWO lines -- ciphertext 1 then ciphertext 2 (27-symbol alphabet, may contain '+')
// stderr: the cleaned plaintext 1 immediately followed by cleaned plaintext 2 (the .solution)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

// Read the first line of a file (or stdin for "-") into alphabet indices (A..Z kept; the
// cube's '+' is a cipher symbol, not used in plaintext). Returns the length.
static int read_plaintext(const char *path, int out[]) {
    FILE *fp = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int c = toupper(ch);
        if (c >= 'A' && c <= 'Z' && n < MAXLEN) out[n++] = g_char_to_idx[c];
    }
    if (fp != stdin) fclose(fp);
    return n;
}

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <pt1|-> <pt2> <keyword> <period1> <period2>\n", argv[0]);
        return 1;
    }
    const char *keyword = argv[3];
    int period1 = atoi(argv[4]);
    int period2 = atoi(argv[5]);
    if (period1 < 1 || period2 < 1) { fprintf(stderr, "periods must be >= 1\n"); return 1; }

    init_alphabet_trifid();                      // 27-symbol alphabet (A..Z + '+')
    if (g_alpha != TRIFID_CELLS) {
        fprintf(stderr, "alphabet is %d symbols, need %d\n", g_alpha, TRIFID_CELLS);
        return 1;
    }
    int side = TRIFID_SIDE;

    static int p1[MAXLEN], p2[MAXLEN];
    int n1 = read_plaintext(argv[1], p1);
    int n2 = read_plaintext(argv[2], p2);
    if (n1 == 0 || n2 == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int kw[256];
    int kwn = 0;
    for (int i = 0; keyword[i] && kwn < 256; i++) {
        int c = toupper((unsigned char) keyword[i]);
        int idx = (c < 128) ? g_char_to_idx[c] : -1;
        if (idx >= 0) kw[kwn++] = idx;
    }
    int cube[TRIFID_CELLS];
    trifid_cube_from_keyword(kw, kwn, cube, g_alpha);

    static int c1[MAXLEN], c2[MAXLEN];
    twin_trifid_encrypt(p1, n1, p2, n2, cube, side, period1, period2, c1, c2);

    for (int i = 0; i < n1; i++) putchar(index_to_char(c1[i]));
    putchar('\n');
    for (int i = 0; i < n2; i++) putchar(index_to_char(c2[i]));
    putchar('\n');

    for (int i = 0; i < n1; i++) fputc(index_to_char(p1[i]), stderr);
    for (int i = 0; i < n2; i++) fputc(index_to_char(p2[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

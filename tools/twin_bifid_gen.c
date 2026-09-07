// Twin Bifid cipher generator (test-data tool, not part of the solver).
//
// A Twin Bifid is TWO Bifid messages enciphered under the SAME keyed 5x5 square (J->I) but
// with DIFFERENT periods; the two plaintexts share a common phrase. This tool reads two
// plaintext files, keeps A..Z (J->I), builds one square from a keyword, and enciphers
// message 1 at period1 and message 2 at period2. It links the REAL cipher code
// (twin_bifid.c + bifid.c + utils.c) so the generator and solver can never drift.
//
//   make twin_bifid_gen
//   ./twin_bifid_gen pt1.txt pt2.txt KEYWORD 7 9 >ciphers.txt 2>twin.solution
//
// argv: <pt1-file|-> <pt2-file> <keyword> <period1> <period2> [omit=J]
// stdout: TWO lines -- ciphertext 1 then ciphertext 2 (bare A..Z over the 25-letter alphabet)
// stderr: the cleaned plaintext 1 immediately followed by cleaned plaintext 2 (the .solution
//         the solver recovers character-for-character; the two decrypts are concatenated)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

static int letter_to_index(int c, char omit) {
    c = toupper(c);
    if (c == toupper((unsigned char) omit)) c = (toupper((unsigned char) omit) == 'J') ? 'I' : 0;
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Read the first line of a file (or stdin for "-") into alphabet indices; returns the length.
static int read_plaintext(const char *path, char omit, int out[]) {
    FILE *fp = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int idx = letter_to_index(ch, omit);
        if (idx >= 0 && n < MAXLEN) out[n++] = idx;
    }
    if (fp != stdin) fclose(fp);
    return n;
}

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <pt1|-> <pt2> <keyword> <period1> <period2> [omit=J]\n", argv[0]);
        return 1;
    }
    const char *keyword = argv[3];
    int period1 = atoi(argv[4]);
    int period2 = atoi(argv[5]);
    char omit = (argc > 6 && argv[6][0]) ? (char) toupper((unsigned char) argv[6][0]) : 'J';
    if (period1 < 1 || period2 < 1) { fprintf(stderr, "periods must be >= 1\n"); return 1; }

    char omit_str[2] = { omit, '\0' };
    init_alphabet(omit_str);                     // 25-letter alphabet, base-25 indices
    if (g_alpha != PLAYFAIR_GRID) {
        fprintf(stderr, "alphabet is %d letters, need %d (one excluded)\n", g_alpha, PLAYFAIR_GRID);
        return 1;
    }
    int side = 5;

    static int p1[MAXLEN], p2[MAXLEN];
    int n1 = read_plaintext(argv[1], omit, p1);
    int n2 = read_plaintext(argv[2], omit, p2);
    if (n1 == 0 || n2 == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int kw[256];
    int kwn = 0;
    for (int i = 0; keyword[i] && kwn < 256; i++) {
        int idx = letter_to_index((unsigned char) keyword[i], omit);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int grid[PLAYFAIR_GRID];
    bifid_grid_from_keyword(kw, kwn, grid, g_alpha);

    static int c1[MAXLEN], c2[MAXLEN];
    twin_bifid_encrypt(p1, n1, p2, n2, grid, side, period1, period2, c1, c2);

    for (int i = 0; i < n1; i++) putchar(index_to_char(c1[i]));
    putchar('\n');
    for (int i = 0; i < n2; i++) putchar(index_to_char(c2[i]));
    putchar('\n');

    // The .solution is the two cleaned plaintexts concatenated (what the solver reports).
    for (int i = 0; i < n1; i++) fputc(index_to_char(p1[i]), stderr);
    for (int i = 0; i < n2; i++) fputc(index_to_char(p2[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

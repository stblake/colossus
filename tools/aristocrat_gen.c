// Aristocrat / Patristocrat cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), builds a 26-letter substitution map from a
// KEYWORD in an ACA keying arrangement (default K2 = keyed ciphertext alphabet), and enciphers it.
// It links the REAL cipher code (aristocrat.c + utils.c) so the generator and the solver can never
// drift in convention.
//
//   make aristocrat_gen
//   ./aristocrat_gen plaintext.txt KRYPTOS               >cipher.txt 2>solution.txt   # Aristocrat
//   ./aristocrat_gen plaintext.txt KRYPTOS patristocrat  >cipher.txt 2>solution.txt   # Patristocrat
//
// argv: <plaintext-file|-> <keyword> [aristocrat|patristocrat] [keying 1..4]
// stdout: the ciphertext -- Aristocrat: letters with single spaces at the word divisions;
//         Patristocrat: letters in 5-letter groups.
// stderr: the solution -- the bare uppercased plaintext letters (what the solver recovers).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLINE (1 << 20)

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <plaintext|-> <keyword> [aristocrat|patristocrat] [keying 1..4]\n", argv[0]);
        return 1;
    }
    bool patristocrat = (argc >= 4 && (strcasecmp(argv[3], "patristocrat") == 0 ||
                                       strcasecmp(argv[3], "patri") == 0 || strcasecmp(argv[3], "pat") == 0));
    int keying = (argc >= 5) ? atoi(argv[4]) : ARIST_K2;
    if (keying < 1 || keying > 4) keying = ARIST_K2;

    init_alphabet(NULL);                         // full 26-letter alphabet

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static char line[MAXLINE];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') if (n < MAXLINE - 1) line[n++] = (char) ch;
    line[n] = '\0';
    if (fp != stdin) fclose(fp);

    // Parse into a LETTER stream (indices 0..25); remember which letters started a new word so the
    // Aristocrat can reproduce the word divisions.
    static int letters[MAXLINE]; static char new_word[MAXLINE];
    int L = 0; bool prev_space = true;
    for (int i = 0; line[i]; i++) {
        unsigned char uc = (unsigned char) line[i];
        if (isspace(uc)) { prev_space = true; continue; }
        int idx = (uc < 128) ? g_char_to_idx[toupper(uc)] : -1;
        if (idx < 0 || idx >= 26) { prev_space = false; continue; }   // drop other punctuation
        new_word[L] = prev_space ? 1 : 0;
        letters[L++] = idx;
        prev_space = false;
    }
    if (L == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    int cmap[26];
    aristocrat_build_map(keying, argv[2], argv[2], 5, cmap);   // kw2==kw1 & shift 5 only matter for K3/K4
    static int ct[MAXLINE];
    aristocrat_apply(letters, L, cmap, ct);

    // Emit the ciphertext.
    if (patristocrat) {
        for (int i = 0; i < L; i++) {
            if (i > 0 && i % 5 == 0) putchar(' ');
            putchar(index_to_char(ct[i]));
        }
    } else {
        for (int i = 0; i < L; i++) {
            if (i > 0 && new_word[i]) putchar(' ');
            putchar(index_to_char(ct[i]));
        }
    }
    putchar('\n');

    // Solution: the bare plaintext letters (what the solver recovers and reports).
    for (int i = 0; i < L; i++) fputc(index_to_char(letters[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

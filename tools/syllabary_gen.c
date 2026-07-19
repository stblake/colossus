// Syllabary cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps LETTERS only, and enciphers with the
// REAL syllabary_encrypt (syllabary.c + utils.c) so the generator and the solver can never drift.
// The 100-token square is a reproducible random scramble (the general "Unknown Coordinates,
// Unknown Keysquare" case) from -sqseed. Tokenization mode selects the isolog spelling.
//
//   make syllabary_gen
//   ./syllabary_gen plaintext.txt -sqseed 5 -mode greedy >cipher.txt 2>solution.txt
//
// argv: <plaintext|->  [-sqseed <n>]  [-mode single|greedy|random]
//   -sqseed : seeds the random square (token arrangement + label permutations). Default 1.
//   -mode   : single (every letter its own token) | greedy (longest tokens) | random (isolog).
//             Default greedy. random additionally uses -sqseed to seed the isolog RNG.
// stdout: the ciphertext (one line of space-separated 2-digit codes).
// stderr: the cleaned plaintext (letters only -- the solution the solver recovers).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "syllabary.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <plaintext|-> [-sqseed <n>] [-mode single|greedy|random]\n", argv[0]);
        return 1;
    }
    unsigned sqseed = 1, mode = SYL_TOK_GREEDY;
    for (int a = 2; a + 1 < argc; a += 2) {
        if (strcmp(argv[a], "-sqseed") == 0) sqseed = (unsigned) atoi(argv[a + 1]);
        else if (strcmp(argv[a], "-mode") == 0) {
            if (strcmp(argv[a + 1], "single") == 0) mode = SYL_TOK_SINGLE;
            else if (strcmp(argv[a + 1], "random") == 0) mode = SYL_TOK_RANDOM;
            else mode = SYL_TOK_GREEDY;
        }
    }
    init_alphabet(NULL);

    SyllabarySquare sq;
    syllabary_build_random(&sq, sqseed);        // seeds the RNG; leaves it seeded for isolog picks

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }

    static int plain[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n')
        if (isalpha(ch) && n < MAXLEN) plain[n++] = toupper(ch) - 'A';
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[MAXLEN];
    int clen = syllabary_encrypt(plain, n, &sq, (int) mode, cipher);

    for (int i = 0; i < clen; i++) printf("%s%02d", i ? " " : "", cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(plain[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

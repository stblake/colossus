// Grandpre cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps LETTERS only, and enciphers with
// the REAL grandpre_encrypt (grandpre.c + utils.c) so the generator and the solver can never
// drift. The square is given as N words of N letters each (rows); ALL 26 letters must appear.
// Homophone choice is random (isolog frequency flattening) -- seed with -seed for reproducibility.
//
//   make grandpre_gen
//   ./grandpre_gen plaintext.txt LADYBUGS AZIMUTHS CALFSKIN QUACKISH UNJOVIAL EVULSION ROWDYISM SEXTUPLY >cipher.txt 2>solution.txt
//
// argv: <plaintext|-> <row-word-1> ... <row-word-N>   (N in 6..10; each word exactly N letters)
//   optional trailing  -seed <n>  fixes the homophone RNG (default 1).
// stdout: the ciphertext (one line of space-separated 2-digit codes).
// stderr: the cleaned plaintext (letters only -- the solution the solver recovers).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "grandpre.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    // Optional trailing "-seed <n>".
    unsigned seed = 1;
    if (argc >= 3 && strcmp(argv[argc - 2], "-seed") == 0) { seed = (unsigned) atoi(argv[argc - 1]); argc -= 2; }

    int N = argc - 2;                       // rows = words after <plaintext>
    if (N < GRANDPRE_MIN_SIDE || N > GRANDPRE_MAX_SIDE) {
        fprintf(stderr, "usage: %s <plaintext|-> <row-word-1> ... <row-word-N>  (N=6..10, each N letters)"
                        "  [-seed <n>]\n", argv[0]);
        return 1;
    }
    init_alphabet(NULL);                    // full 26-letter A..Z
    seed_rand(seed);

    const char *words[GRANDPRE_MAX_SIDE];
    for (int r = 0; r < N; r++) words[r] = argv[2 + r];

    GrandpreSquare g;
    if (grandpre_build_from_words(&g, N, words) != 0) {
        fprintf(stderr, "square build failed: need %d words of exactly %d letters and ALL 26 letters present\n", N, N);
        return 1;
    }

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }

    static int plain[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n')
        if (isalpha(ch) && n < MAXLEN) plain[n++] = toupper(ch) - 'A';
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[MAXLEN];
    int clen = grandpre_encrypt(plain, n, &g, cipher);

    for (int i = 0; i < clen; i++) printf("%s%02d", i ? " " : "", cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(plain[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

// Pollux cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and enciphers with the
// REAL pollux_encrypt (pollux.c + utils.c) so the generator and the solver can never drift.
// The key is the digit->{dot,dash,x} assignment given as a 10-char string of '.'/'-'/'x'
// indexed by digit 0..9. Encryption is polyphonic (each Morse symbol -> a random assigned
// digit), so a fixed RNG seed makes the ciphertext reproducible.
//
//   make pollux_gen
//   ./pollux_gen plaintext.txt ".x-..x.--x" 42 >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <assignment> [seed]
// stdout: the Pollux ciphertext (one line of digits 0-9; length differs from N)
// stderr: the cleaned plaintext (the solution: bare A..Z, what the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "pollux.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <plaintext|-> <assignment(10x ./-/x)> [seed]\n", argv[0]);
        return 1;
    }

    init_alphabet(NULL);                         // full 26-letter A..Z alphabet
    if (g_alpha != 26) {
        fprintf(stderr, "alphabet is %d symbols, need 26\n", g_alpha);
        return 1;
    }

    // Parse the digit->element assignment (indexed by digit 0..9).
    const char *a = argv[2];
    if ((int) strlen(a) != 10) {
        fprintf(stderr, "assignment must be exactly 10 chars of ./-/x (got \"%s\")\n", a);
        return 1;
    }
    int key[10], cnt[3] = {0, 0, 0};
    for (int d = 0; d < 10; d++) {
        char c = a[d];
        int e = (c == '.') ? PX_DOT : (c == '-') ? PX_DASH
              : (c == 'x' || c == 'X') ? PX_X : -1;
        if (e < 0) { fprintf(stderr, "bad assignment char '%c' (use ./-/x)\n", c); return 1; }
        key[d] = e; cnt[e]++;
    }
    if (cnt[PX_DOT] == 0 || cnt[PX_DASH] == 0 || cnt[PX_X] == 0) {
        fprintf(stderr, "assignment needs at least one of each of . - x\n");
        return 1;
    }

    seed_rand(argc > 3 ? (unsigned) strtoul(argv[3], NULL, 10) : 1u);

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int c = toupper(ch);
        if (c >= 'A' && c <= 'Z' && n < MAXLEN) raw[n++] = g_char_to_idx[c];
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[6 * MAXLEN];
    int clen = pollux_encrypt(raw, n, key, cipher);
    if (clen < 0) { fprintf(stderr, "encrypt failed (bad key)\n"); return 1; }

    for (int i = 0; i < clen; i++) putchar('0' + cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

// Morbit cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and enciphers with the
// REAL morbit_encrypt (morbit.c + utils.c) so the generator and the solver can never drift.
// The key is a 9-LETTER KEYWORD (ACA style): the digits 1..9 are assigned to the 9 base-3
// pairs by the keyword letters' alphabetical rank (stable for duplicate letters). Encryption
// is deterministic (a pair<->digit bijection), so no RNG is used; a [seed] argument is
// accepted for interface parity with the other generators but ignored.
//
//   make morbit_gen
//   ./morbit_gen plaintext.txt WISECRACK >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <9-letter keyword> [seed]
// stdout: the Morbit ciphertext (one line of digits 1-9; length differs from N)
// stderr: the cleaned plaintext (the solution: bare A..Z, what the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "morbit.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <plaintext|-> <9-letter keyword> [seed]\n", argv[0]);
        return 1;
    }

    init_alphabet(NULL);                         // full 26-letter A..Z alphabet
    if (g_alpha != 26) {
        fprintf(stderr, "alphabet is %d symbols, need 26\n", g_alpha);
        return 1;
    }

    // Parse the 9-letter keyword into alphabet indices.
    const char *kwarg = argv[2];
    if ((int) strlen(kwarg) != 9) {
        fprintf(stderr, "keyword must be exactly 9 letters (got \"%s\")\n", kwarg);
        return 1;
    }
    int kw[9];
    for (int i = 0; i < 9; i++) {
        int c = toupper((unsigned char) kwarg[i]);
        if (c < 'A' || c > 'Z') { fprintf(stderr, "bad keyword char '%c'\n", kwarg[i]); return 1; }
        kw[i] = g_char_to_idx[c];
    }

    // Derive pair -> digit by stable alphabetical rank (the ACA "number" row), then
    // invert to the digit -> pair key morbit_encrypt expects. Ranks of 9 items always
    // form a permutation of 1..9, so the key is always a valid bijection.
    int key[10]; key[0] = 0;
    for (int i = 0; i < 9; i++) {
        int rank = 1;
        for (int j = 0; j < 9; j++)
            if (kw[j] < kw[i] || (kw[j] == kw[i] && j < i)) rank++;
        key[rank] = i;                           // digit `rank` -> pair `i`
    }

    (void) (argc > 3);                           // [seed] accepted but unused (deterministic)

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
    int clen = morbit_encrypt(raw, n, key, cipher);
    if (clen < 0) { fprintf(stderr, "encrypt failed (bad key)\n"); return 1; }

    for (int i = 0; i < clen; i++) putchar('0' + cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

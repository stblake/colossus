// Condi cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and enciphers it with a Condi
// cipher: a plaintext-feedback substitution over a keyed alphabet sigma. The shift for each letter
// is the position of the preceding plaintext letter in sigma; the first letter uses a starter
// offset. It links the REAL cipher code (condi.c + utils.c), so the generator and the solver can
// never drift.
//
//   make condi_gen
//   ./condi_gen plaintext.txt STRANGE 25 >cipher.txt 2>solution.txt
//   ./condi_gen plaintext.txt STRANGE 25 21 >cipher.txt 2>solution.txt   # keyed alphabet shifted by 21
//
// argv: <plaintext-file|-> <keyword> <starter> [shift]
//   starter : the initial offset # (0..25)
//   shift   : optional cyclic left-rotation of the keyed alphabet (0..25, default 0; the ACA
//             "the keyed alphabet may or may not be shifted"). shift=21 turns the STRANGE keyed
//             alphabet STRANGEBCDF... into the worked-example VWXYZSTRANGE...
// stdout: the ciphertext (one line, bare A..Z over the 26-letter alphabet)
// stderr: line 1 = the plaintext (the solution: bare A..Z, what the solver recovers char-for-char)
//         line 2 = the keyed alphabet sigma used, and the starter

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <plaintext|-> <keyword> <starter> [shift]\n"
                        "  starter : initial offset # (0..25)\n"
                        "  shift   : cyclic left-rotation of the keyed alphabet (0..25, default 0)\n",
                        argv[0]);
        return 1;
    }
    const char *keyword_str = argv[2];
    int starter = ((atoi(argv[3]) % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    int shift = (argc >= 5) ? (((atoi(argv[4]) % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE) : 0;

    init_alphabet(NULL);                          // full 26-letter alphabet
    if (g_alpha != ALPHABET_SIZE) {
        fprintf(stderr, "alphabet is %d letters, need 26\n", g_alpha);
        return 1;
    }

    // Build the keyed alphabet (keyword + remaining letters), then apply the optional cyclic shift.
    int keyed[ALPHABET_SIZE], sigma[ALPHABET_SIZE], sigma_inv[ALPHABET_SIZE];
    make_keyed_alphabet((char *) keyword_str, keyed);
    for (int k = 0; k < ALPHABET_SIZE; k++) sigma[k] = keyed[(k + shift) % ALPHABET_SIZE];
    for (int k = 0; k < ALPHABET_SIZE; k++) sigma_inv[sigma[k]] = k;

    // Read the plaintext (first line of the file, or stdin), to alphabet indices.
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

    static int cipher[MAXLEN];
    condi_encrypt(raw, n, sigma, sigma_inv, starter, cipher);

    for (int i = 0; i < n; i++) putchar(index_to_char(cipher[i]));
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    fprintf(stderr, "sigma=");
    for (int k = 0; k < ALPHABET_SIZE; k++) fputc(index_to_char(sigma[k]), stderr);
    fprintf(stderr, " starter=%d\n", starter);
    return 0;
}

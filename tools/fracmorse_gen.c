// Fractionated Morse cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, builds the keyed 26-letter
// alphabet from KEYWORD (keyword letters, duplicates dropped, then the rest of the alphabet
// ascending -- exactly what the solver searches), and enciphers with fracmorse_encrypt. It
// links the REAL cipher code (fracmorse.c + utils.c) so the generator and the solver can never
// drift in convention.
//
//   make fracmorse_gen
//   ./fracmorse_gen plaintext.txt KEYWORD >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <keyword>
// stdout: the Fractionated Morse ciphertext (one line, bare A..Z; C differs from N)
// stderr: the cleaned plaintext (the solution: bare A..Z, what the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

// Build the keyed alphabet sigma (rank -> letter): the keyword letters in order with duplicates
// removed, then the remaining alphabet letters ascending.
static void build_keyed_alphabet(const char *keyword, int sigma[]) {
    char used[26];
    for (int i = 0; i < 26; i++) used[i] = 0;
    int m = 0;
    for (int i = 0; keyword[i] && m < 26; i++) {
        int c = toupper((unsigned char) keyword[i]);
        if (c < 'A' || c > 'Z') continue;
        int l = g_char_to_idx[c];
        if (used[l]) continue;
        used[l] = 1;
        sigma[m++] = l;
    }
    for (int l = 0; l < 26 && m < 26; l++)
        if (!used[l]) { used[l] = 1; sigma[m++] = l; }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <plaintext|-> <keyword>\n", argv[0]);
        return 1;
    }

    init_alphabet(NULL);                         // full 26-letter A..Z alphabet
    if (g_alpha != 26) {
        fprintf(stderr, "alphabet is %d symbols, need 26\n", g_alpha);
        return 1;
    }

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

    int sigma[26];
    build_keyed_alphabet(argv[2], sigma);

    static int cipher[2 * MAXLEN];
    int clen = fracmorse_encrypt(raw, n, sigma, cipher);

    for (int i = 0; i < clen; i++) putchar(index_to_char(cipher[i]));
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

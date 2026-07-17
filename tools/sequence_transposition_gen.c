// Sequence Transposition cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and enciphers it with a
// Sequence Transposition cipher. It links the REAL cipher code (sequence_transposition.c +
// gromark.c + utils.c), so the generator and the solver can never drift in convention.
//
//   make sequence_transposition_gen
//   ./sequence_transposition_gen plaintext.txt GUMMYBEARS 69315 >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <10-letter-keyword> <primer-digits>
//   The keyword's 10 letters are ranked 1..10 alphabetically (10 -> digit 0) to give the column
//   read order; the primer seeds the chain-addition digit sequence that assigns each plaintext
//   letter to a column.
// stdout: the ciphertext (one line, bare A..Z).
// stderr: the plaintext solution (bare A..Z) preceded by a parameter line.

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <plaintext|-> <10-letter-keyword> <primer-digits>\n", argv[0]);
        return 1;
    }
    const char *keyword = argv[2];
    const char *primerstr = argv[3];

    init_alphabet(NULL);                          // full 26-letter alphabet
    if (g_alpha != ALPHABET_SIZE) {
        fprintf(stderr, "alphabet is %d letters, need 26\n", g_alpha);
        return 1;
    }

    int primer[SEQ_TRANS_MAX_PRIMER], P = 0;
    for (int i = 0; primerstr[i] && P < SEQ_TRANS_MAX_PRIMER; i++)
        if (isdigit((unsigned char) primerstr[i])) primer[P++] = primerstr[i] - '0';
    if (P < 2) { fprintf(stderr, "primer needs at least 2 digits\n"); return 1; }

    int kwletters = 0;
    for (int i = 0; keyword[i]; i++) if (isalpha((unsigned char) keyword[i])) kwletters++;
    if (kwletters < SEQ_TRANS_BUCKETS)
        fprintf(stderr, "[warning: keyword has %d letters; Sequence Transposition uses 10]\n",
                kwletters);

    int pi[SEQ_TRANS_BUCKETS];
    sequence_transposition_pi_from_keyword(keyword, pi);

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
    sequence_transposition_encrypt(raw, n, primer, P, pi, cipher);

    for (int i = 0; i < n; i++) putchar(index_to_char(cipher[i]));
    putchar('\n');

    fprintf(stderr, "[sequence-transposition: keyword=%s primer=", keyword);
    for (int i = 0; i < P; i++) fprintf(stderr, "%d", primer[i]);
    fprintf(stderr, " read-order=[");
    for (int k = 0; k < SEQ_TRANS_BUCKETS; k++) fprintf(stderr, "%s%d", k ? " " : "", pi[k]);
    fprintf(stderr, "]]\n");
    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

// Baconian cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and biliteral-encodes it
// with the REAL baconian_encrypt (baconian.c + utils.c) so the generator and the solver can
// never drift. It then CONCEALS the a/b bit stream in cover text under a canonical classifier
// (a binary rule + polarity) in one of two grouping modes:
//   letter : each bit -> one random cover letter whose class matches the bit (spaced into
//            cosmetic 5-letter groups; the solver's per-letter mode ignores the spaces);
//   word   : each bit -> one random word whose FIRST letter's class matches the bit.
// Concealment is random (a fixed seed makes it reproducible). The decode folds J->I and V->U,
// so the reference SOLUTION written to stderr is the folded plaintext (what the solver recovers).
//
//   make baconian_gen
//   ./baconian_gen plaintext.txt word AM 0 42 >cover.txt 2>solution.txt
//
// argv: <plaintext-file|-> <mode:letter|word> <rule:AM|VOWEL> <polarity:0|1> [seed]
// stdout: the Baconian cover text
// stderr: the folded plaintext (the solution the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "baconian.h"

#define MAXLEN (1 << 16)

static int is_vowel(int l) { return (l == 0 || l == 4 || l == 8 || l == 14 || l == 20); }

// The a/b label of cover letter l under (rule, polarity): 0 = group a, 1 = group b.
static int classify(int l, int rule_am, int pol) {
    int c = rule_am ? ((l >= 13) ? 1 : 0) : (is_vowel(l) ? 0 : 1);
    return pol ? (c ^ 1) : c;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s <plaintext|-> <letter|word> <AM|VOWEL> <0|1> [seed]\n", argv[0]);
        return 1;
    }
    init_alphabet(NULL);
    if (g_alpha != 26) { fprintf(stderr, "alphabet is %d symbols, need 26\n", g_alpha); return 1; }

    int word_mode = (strcmp(argv[2], "word") == 0);
    if (!word_mode && strcmp(argv[2], "letter") != 0) {
        fprintf(stderr, "mode must be letter or word (got \"%s\")\n", argv[2]); return 1;
    }
    int rule_am = (strcasecmp(argv[3], "AM") == 0);
    if (!rule_am && strcasecmp(argv[3], "VOWEL") != 0) {
        fprintf(stderr, "rule must be AM or VOWEL (got \"%s\")\n", argv[3]); return 1;
    }
    int pol = (atoi(argv[4]) != 0);
    seed_rand(argc > 5 ? (unsigned) strtoul(argv[5], NULL, 10) : 1u);

    // Read + clean the plaintext.
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

    // Biliteral-encode with the real primitive.
    static int bits[5 * MAXLEN];
    int nbits = baconian_encrypt(raw, n, bits);

    // Precompute the two class lists (letters in group a / group b under this classifier).
    int grp[2][26], ng[2] = {0, 0};
    for (int l = 0; l < 26; l++) grp[classify(l, rule_am, pol)][ng[classify(l, rule_am, pol)]++] = l;
    if (ng[0] == 0 || ng[1] == 0) { fprintf(stderr, "degenerate classifier (empty group)\n"); return 1; }

    // Emit the cover text.
    for (int i = 0; i < nbits; i++) {
        int b = bits[i];
        if (word_mode) {
            int wlen = 2 + (int) rand_bounded(5);                 // 2..6 letters
            putchar(index_to_char(grp[b][rand_bounded(ng[b])]));  // first letter carries the bit
            for (int k = 1; k < wlen; k++) putchar(index_to_char(rand_bounded(26)));
            if (i + 1 < nbits) putchar(' ');
        } else {
            putchar(index_to_char(grp[b][rand_bounded(ng[b])]));
            if ((i + 1) % 5 == 0 && i + 1 < nbits) putchar(' ');  // cosmetic 5-letter groups
        }
    }
    putchar('\n');

    // The folded plaintext (J->I, V->U): exactly what the solver recovers.
    for (int i = 0; i < n; i++) {
        int l = raw[i];
        if (l == 9) l = 8;          // J -> I
        else if (l == 21) l = 20;   // V -> U
        fputc(index_to_char(l), stderr);
    }
    fputc('\n', stderr);
    return 0;
}

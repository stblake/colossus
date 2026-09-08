// Chaocipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps only A..Z, and enciphers it
// with a Chaocipher whose key is a pair of random STARTING alphabets derived from a seed.
// It links the REAL cipher code (chaocipher.c + utils.c), so the generator and the solver
// can never drift in convention.
//
//   make chaocipher_gen
//   ./chaocipher_gen plaintext.txt 1              >cipher.txt 2>solution.txt
//   ./chaocipher_gen plaintext.txt LEFTALPHA... RIGHTALPHA...   (explicit 26-letter keys)
//
// argv: <plaintext-file|-> [seed=1 | LEFT26 RIGHT26]
// stdout: the ciphertext (one line, bare A..Z).
// stderr: a "[chaocipher: left=... right=...]" key line, then the cleaned A..Z plaintext
//         solution (bare A..Z, what the solver recovers).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"
#include "chaocipher.h"

#define MAXLEN (1 << 20)

// Parse a 26-letter A..Z permutation string into alphabet indices; returns 1 on success.
static int parse_alphabet(const char *s, int a[CHAO_N]) {
    int seen[CHAO_N] = {0};
    if ((int) strlen(s) != CHAO_N) return 0;
    for (int i = 0; i < CHAO_N; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c < 'A' || c > 'Z') return 0;
        int idx = g_char_to_idx[c];
        if (idx < 0 || idx >= CHAO_N || seen[idx]) return 0;
        seen[idx] = 1;
        a[i] = idx;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <plaintext|-> [seed | LEFT26 RIGHT26]\n", argv[0]);
        return 1;
    }

    init_alphabet(NULL);                          // full 26-letter A..Z alphabet
    if (g_alpha != CHAO_N) {
        fprintf(stderr, "alphabet is %d letters, need %d\n", g_alpha, CHAO_N);
        return 1;
    }

    int left[CHAO_N], right[CHAO_N];
    if (argc >= 4) {                              // explicit LEFT / RIGHT starting alphabets
        if (!parse_alphabet(argv[2], left) || !parse_alphabet(argv[3], right)) {
            fprintf(stderr, "LEFT/RIGHT must each be a 26-letter A..Z permutation\n");
            return 1;
        }
    } else {                                      // random alphabets from a seed
        unsigned seed = (argc > 2) ? (unsigned) atoi(argv[2]) : 1u;
        seed_rand(seed);
        for (int i = 0; i < CHAO_N; i++) { left[i] = i; right[i] = i; }
        shuffle(left, CHAO_N);
        shuffle(right, CHAO_N);
    }

    // Read the plaintext (first line of the file, or stdin) to alphabet indices.
    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int c = toupper(ch);
        if (c >= 'A' && c <= 'Z' && n < MAXLEN) {
            int idx = g_char_to_idx[c];
            if (idx >= 0) raw[n++] = idx;
        }
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[MAXLEN];
    chaocipher_encrypt(raw, n, left, right, cipher);

    for (int i = 0; i < n; i++) putchar(index_to_char(cipher[i]));
    putchar('\n');

    fprintf(stderr, "[chaocipher: left=");
    for (int i = 0; i < CHAO_N; i++) fputc(index_to_char(left[i]), stderr);
    fprintf(stderr, " right=");
    for (int i = 0; i < CHAO_N; i++) fputc(index_to_char(right[i]), stderr);
    fprintf(stderr, "]\n");
    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

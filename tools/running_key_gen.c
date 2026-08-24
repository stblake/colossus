// Running Key cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and enciphers it as a
// Running Key. It links the REAL cipher code (running_key.c + utils.c), so the generator
// and the solver can never drift in convention.
//
//   make running_key_gen
//   # ACA self-keyed: the passage's first half keys its second half; CT is half as long.
//   ./running_key_gen passage.txt vig            >cipher.txt 2>solution.txt
//   # independent key: encipher the whole plaintext against a separate key TEXT.
//   ./running_key_gen plaintext.txt beau key.txt >cipher.txt 2>solution.txt
//
// argv: <plaintext|-> <family: vig|beau|variant|porta> [keyfile]
//   No keyfile  => ACA self-keyed. The plaintext is one passage of length 2N (trimmed to
//                  even); key = first N letters, plaintext half = last N; CT has N letters.
//                  The stderr solution is the WHOLE 2N passage (what the solver reports).
//   With keyfile => independent key. The key text (letters only, >= plaintext length)
//                  enciphers the whole plaintext; CT has plaintext length. The stderr
//                  solution is the plaintext.
// stdout: the ciphertext (one line, bare A..Z).
// stderr: a parameter line, then the solution plaintext (bare A..Z).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "running_key.h"

#define MAXLEN (1 << 20)

static int family_from_arg(const char *s) {
    if (!strcasecmp(s, "vig") || !strcasecmp(s, "vigenere")) return RK_VIGENERE;
    if (!strcasecmp(s, "beau") || !strcasecmp(s, "beaufort")) return RK_BEAUFORT;
    if (!strcasecmp(s, "variant") || !strcasecmp(s, "var")) return RK_VARIANT;
    if (!strcasecmp(s, "porta")) return RK_PORTA;
    return atoi(s);   // also accept 0..3
}

static int read_letters(const char *path, int out[], int cap) {
    FILE *fp = (strcmp(path, "-") == 0) ? stdin : fopen(path, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); return -1; }
    int n = 0, ch;
    // Read the whole first line (stops at newline) for the plaintext; a keyfile may span
    // multiple lines, so read to EOF when not stdin.
    while ((ch = fgetc(fp)) != EOF && n < cap) {
        if (ch == '\n' && fp == stdin) break;
        int c = toupper((unsigned char) ch);
        if (c >= 'A' && c <= 'Z') out[n++] = g_char_to_idx[c];
    }
    if (fp != stdin) fclose(fp);
    return n;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <plaintext|-> <family: vig|beau|variant|porta> [keyfile]\n",
                argv[0]);
        return 1;
    }
    init_alphabet(NULL);                          // full 26-letter alphabet
    if (g_alpha != ALPHABET_SIZE) {
        fprintf(stderr, "alphabet is %d letters, need 26\n", g_alpha);
        return 1;
    }
    int family = family_from_arg(argv[2]);
    if (family < 0 || family >= RK_N_FAMILIES) { fprintf(stderr, "bad family '%s'\n", argv[2]); return 1; }

    static int raw[MAXLEN];
    int L = read_letters(argv[1], raw, MAXLEN);
    if (L <= 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[MAXLEN];

    if (argc >= 4) {
        // Independent-key mode: encipher the whole plaintext against the key text.
        static int key[MAXLEN];
        int K = read_letters(argv[3], key, MAXLEN);
        if (K < L) { fprintf(stderr, "key text has %d letters, need >= %d\n", K, L); return 1; }
        running_key_encrypt(cipher, raw, L, key, family);
        for (int i = 0; i < L; i++) putchar(index_to_char(cipher[i]));
        putchar('\n');
        fprintf(stderr, "[running-key: family=%s mode=independent N=%d]\n", rk_family_name(family), L);
        for (int i = 0; i < L; i++) fputc(index_to_char(raw[i]), stderr);
        fputc('\n', stderr);
    } else {
        // ACA self-keyed mode: first half keys the second half.
        int N = L / 2;                            // trim to even
        if (N < 1) { fprintf(stderr, "passage too short\n"); return 1; }
        const int *key = raw;                     // first N letters
        const int *pt  = raw + N;                 // next N letters
        running_key_encrypt(cipher, pt, N, key, family);
        for (int i = 0; i < N; i++) putchar(index_to_char(cipher[i]));
        putchar('\n');
        fprintf(stderr, "[running-key: family=%s mode=self-keyed N=%d (passage=%d)]\n",
                rk_family_name(family), N, 2 * N);
        for (int i = 0; i < 2 * N; i++) fputc(index_to_char(raw[i]), stderr);  // full passage = solution
        fputc('\n', stderr);
    }
    return 0;
}

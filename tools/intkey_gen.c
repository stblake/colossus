// Interrupted Key cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, and enciphers it with an
// Interrupted Key cipher: a periodic base cipher (Vigenere / Variant / Beaufort) under a letter
// keyword whose key index RESETS to the first key letter at break points. It links the REAL
// cipher code (intkey.c + utils.c), so the generator and the solver can never drift.
//
//   make intkey_gen
//   ./intkey_gen plaintext.txt ORANGE vig ct N >cipher.txt 2>solution.txt
//   ./intkey_gen plaintext.txt ORANGE vig random 12345 >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <keyword> <vig|var|beau> <ct|pt|random> <param>
//   scheme ct/pt : <param> is the interruptor letter (reset AFTER that ct/pt letter)
//   scheme random: <param> is an integer RNG seed (group lengths random 1..P; the first group is
//                  length P so the entire keyword is used at least once, per the ACA rule)
// stdout: the ciphertext (one line, bare A..Z over the 26-letter alphabet)
// stderr: line 1 = the plaintext (the solution: bare A..Z, what the solver recovers char-for-char)
//         line 2 = the group-start break positions (space-separated 0-based), for -breaks

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s <plaintext|-> <keyword> <vig|var|beau> <ct|pt|random> <param>\n"
                        "  ct/pt : <param> = interruptor letter\n"
                        "  random: <param> = RNG seed\n", argv[0]);
        return 1;
    }
    const char *keyword_str = argv[2];
    const char *base_str = argv[3];
    const char *scheme_str = argv[4];
    const char *param = argv[5];

    int base = IK_BASE_VIG;
    if (strcmp(base_str, "var") == 0 || strcmp(base_str, "variant") == 0) base = IK_BASE_VAR;
    else if (strcmp(base_str, "beau") == 0 || strcmp(base_str, "beaufort") == 0) base = IK_BASE_BEAU;
    else if (strcmp(base_str, "vig") != 0 && strcmp(base_str, "vigenere") != 0) {
        fprintf(stderr, "unknown base '%s' (use vig|var|beau)\n", base_str);
        return 1;
    }

    init_alphabet(NULL);                          // full 26-letter alphabet
    if (g_alpha != ALPHABET_SIZE) {
        fprintf(stderr, "alphabet is %d letters, need 26\n", g_alpha);
        return 1;
    }

    // Parse the keyword into per-column base shifts 0..25 (each key letter is its own shift).
    static int keyword[MAX_CYCLEWORD_LEN];
    int P = 0;
    for (int i = 0; keyword_str[i] && P < MAX_CYCLEWORD_LEN; i++) {
        int c = toupper((unsigned char) keyword_str[i]);
        if (c >= 'A' && c <= 'Z') keyword[P++] = g_char_to_idx[c];
    }
    if (P == 0) { fprintf(stderr, "keyword has no letters\n"); return 1; }

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
    static int is_break[MAXLEN];
    for (int i = 0; i < n; i++) is_break[i] = 0;

    if (strcmp(scheme_str, "ct") == 0 || strcmp(scheme_str, "pt") == 0) {
        int c = toupper((unsigned char) param[0]);
        if (c < 'A' || c > 'Z') { fprintf(stderr, "interruptor must be a letter A..Z\n"); return 1; }
        int interruptor = g_char_to_idx[c];
        if (strcmp(scheme_str, "ct") == 0) {
            intkey_encrypt_ctint(cipher, raw, n, keyword, P, base, interruptor);
            intkey_build_mask_ct(is_break, cipher, n, interruptor);          // for the emitted breaks
        } else {
            intkey_encrypt_ptint(cipher, raw, n, keyword, P, base, interruptor);
            int k = 0;                                                       // recover pt-break mask
            for (int i = 0; i < n; i++) {
                is_break[i] = (i > 0 && raw[i - 1] == interruptor) ? 1 : 0;
                if (raw[i] == interruptor) k = 0; else if (++k >= P) k = 0;
            }
        }
    } else if (strcmp(scheme_str, "random") == 0) {
        seed_rand((uint32_t) strtoul(param, NULL, 10));
        int pos = 0, first = 1;
        while (pos < n) {
            is_break[pos] = 1;                                               // group start
            int L = first ? P : rand_int(1, P + 1);                         // first group = full keyword
            first = 0;
            pos += L;
        }
        intkey_encrypt_mask(cipher, raw, n, keyword, P, base, is_break);
    } else {
        fprintf(stderr, "unknown scheme '%s' (use ct|pt|random)\n", scheme_str);
        return 1;
    }

    for (int i = 0; i < n; i++) putchar(index_to_char(cipher[i]));
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    for (int i = 0; i < n; i++) if (is_break[i]) fprintf(stderr, "%d ", i);
    fputc('\n', stderr);
    return 0;
}

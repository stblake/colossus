// Ragbaby cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin) WITH its word divisions (spaces matter --
// they drive the Ragbaby numbering), builds the 24-letter keyed alphabet from KEYWORD (keyword
// letters, duplicates dropped, then the rest of the 24-letter alphabet ascending -- exactly what
// the solver searches), and enciphers with ragbaby_encrypt. Plaintext J/X fold to I/W (the ACA
// I/J, W/X pairing). It links the REAL cipher code (ragbaby.c + utils.c) so the generator and the
// solver can never drift in convention.
//
//   make ragbaby_gen
//   ./ragbaby_gen plaintext.txt GROSBEAK >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> <keyword>
// stdout: the Ragbaby ciphertext (one line, letters + single spaces at the word divisions)
// stderr: the solution -- the bare folded plaintext (letters only, J->I / X->W), what the
//         solver recovers and reports.

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLINE (1 << 20)
#define RAG_ALPHA 24

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <plaintext|-> <keyword>\n", argv[0]);
        return 1;
    }

    init_alphabet_ragbaby();                     // 24-letter alphabet (A..Z minus J,X; I/J, W/X paired)
    if (g_alpha != RAG_ALPHA) {
        fprintf(stderr, "alphabet is %d symbols, need 24\n", g_alpha);
        return 1;
    }

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static char line[MAXLINE];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') if (n < MAXLINE - 1) line[n++] = (char) ch;
    line[n] = '\0';
    if (fp != stdin) fclose(fp);

    // Parse the spaced line into the LETTER stream + per-letter numbers, then encrypt.
    static int letters[MAXLINE], num[MAXLINE], ct[MAXLINE];
    int L = 0;
    ragbaby_number_stream(line, RAG_ALPHA, letters, num, &L);
    if (L == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    int ka[RAG_ALPHA], ka_inv[RAG_ALPHA];
    ragbaby_build_keyed_alphabet(ka, argv[2], RAG_ALPHA);
    for (int p = 0; p < RAG_ALPHA; p++) ka_inv[ka[p]] = p;

    ragbaby_encrypt(letters, num, L, ka, ka_inv, RAG_ALPHA, ct);

    // Emit the ciphertext: walk the original line, replacing each letter with its cipher letter
    // (consuming ct[] in order) and collapsing each whitespace run to a single space (the word
    // division). Other punctuation is dropped (the solver skips it anyway).
    int idx = 0; bool prev_space = false, started = false;
    for (const char *q = line; *q; q++) {
        unsigned char uc = (unsigned char) *q;
        if (isspace(uc)) { prev_space = true; continue; }
        int a = (uc < 128) ? g_char_to_idx[toupper(uc)] : -1;
        if (a < 0 || a >= RAG_ALPHA) continue;    // non-letter: drop
        if (prev_space && started) putchar(' ');
        putchar(index_to_char(ct[idx++]));
        prev_space = false; started = true;
    }
    putchar('\n');

    // Solution: the bare folded plaintext (letters only) -- what the solver recovers.
    for (int i = 0; i < L; i++) fputc(index_to_char(letters[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

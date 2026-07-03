// Period column order transposition generator (test-data tool, not part of the solver).
//
// Applies one or more AZdecrypt "Period column order" stages (forward) to a plaintext
// and prints the resulting ciphertext. It links the REAL primitive (period_column.c),
// so the generator and the solver can never drift in convention. Spaces / punctuation
// are carried as genuine grid cells (the transposition permutes them like letters), so
// the plaintext may contain them and the ciphertext preserves the exact multiset.
//
//   make period_column_gen
//   ./period_column_gen plaintext.txt 4 3 tp 56 2 utp   >cipher.txt 2>solution.txt
//
// argv: <plaintext-file|-> ( <dx> <period> <tp|utp> )+
//   Each triple is one stage (grid width dx, period, transpose/untranspose), applied
//   left to right. To make cipher C such that the solver's forward decode recovers the
//   plaintext, pass the DECRYPT stages here in reverse with flipped directions -- or
//   simply pass any stages and let the solver find the inverse composition.
// stdout: the ciphertext (one line).
// stderr: the plaintext (what the solver should recover).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN MAX_CIPHER_LENGTH

int main(int argc, char **argv) {
    if (argc < 5 || (argc - 2) % 3 != 0) {
        fprintf(stderr, "usage: %s <plaintext|-> ( <dx> <period> <tp|utp> )+\n", argv[0]);
        return 1;
    }

    // Read the plaintext (first line of the file, or stdin), carrying every character
    // opaquely as its byte value (letters AND spaces/punctuation are real grid cells).
    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int a[MAXLEN], b[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n' && n < MAXLEN)
        a[n++] = ch;
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    // Remember the plaintext for the solution line.
    static int pt[MAXLEN];
    for (int i = 0; i < n; i++) pt[i] = a[i];

    // Apply each stage forward, ping-ponging between the two buffers.
    int *cur = a, *nxt = b;
    int nstages = (argc - 2) / 3;
    for (int s = 0; s < nstages; s++) {
        int dx  = atoi(argv[2 + 3 * s]);
        int p   = atoi(argv[3 + 3 * s]);
        int utp = (strcmp(argv[4 + 3 * s], "utp") == 0 || strcmp(argv[4 + 3 * s], "UTP") == 0) ? 1 : 0;
        if (dx < 1 || dx > n || p < 1) {
            fprintf(stderr, "stage %d: bad dx=%d p=%d (need 1<=dx<=%d, p>=1)\n", s + 1, dx, p, n);
            return 1;
        }
        period_column_transform(cur, nxt, n, dx, p, utp);
        int *t = cur; cur = nxt; nxt = t;
    }

    for (int i = 0; i < n; i++) putchar(cur[i]);
    putchar('\n');
    fprintf(stderr, "[period-column: %d stage(s), %d chars]\n", nstages, n);
    for (int i = 0; i < n; i++) fputc(pt[i], stderr);
    fputc('\n', stderr);
    return 0;
}

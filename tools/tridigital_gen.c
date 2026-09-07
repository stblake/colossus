// Tridigital cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps LETTERS and WORD BREAKS (any run
// of non-letters collapses to a single word separator; leading/trailing breaks are dropped),
// and enciphers with the REAL tridigital_encrypt (tridigital.c + utils.c) so the generator and
// the solver can never drift.
//
//   make tridigital_gen
//   ./tridigital_gen plaintext.txt NOVELCRAFT DRAGONFLY >cipher.txt 2>solution.txt
//
// argv: <plaintext|-> <column-keyword (10 distinct letters)> <alphabet-keyword> [--trailsep]
//   column-keyword   : 10 distinct letters -> the column digit-labels (alphabetical rank, 10->0).
//   alphabet-keyword : any keyword -> the 26-letter keyed alphabet (deduped keyword + A..Z tail).
//   --trailsep       : also emit a separator digit after the final word (default: internal only).
// stdout: the ciphertext (one line of digits 0-9; same length as the cleaned plaintext).
// stderr: the cleaned spaced plaintext (the solution the solver recovers).

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "tridigital.h"

#define MAXLEN (1 << 18)

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <plaintext|-> <column-keyword-10-letters> <alphabet-keyword> "
                        "[--trailsep]\n", argv[0]);
        return 1;
    }
    init_alphabet(NULL);                 // full 26-letter A..Z (for index_to_char on stderr)

    TridigitalGrid g;
    if (tridigital_build_from_keywords(&g, argv[2], argv[3]) != 0) {
        fprintf(stderr, "grid build failed: the column-keyword needs 10 DISTINCT letters\n");
        return 1;
    }
    int trailsep = (argc > 4 && strcmp(argv[4], "--trailsep") == 0);

    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }

    static int plain[MAXLEN];
    int n = 0, pending = 0, have = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        if (isalpha(ch)) {
            if (pending && have && n < MAXLEN) plain[n++] = TRI_SPACE;   // one break between words
            pending = 0;
            if (n < MAXLEN) { plain[n++] = toupper(ch) - 'A'; have = 1; }
        } else {
            pending = 1;                                                 // collapse a run of non-letters
        }
    }
    if (fp != stdin) fclose(fp);
    if (trailsep && have && n < MAXLEN) plain[n++] = TRI_SPACE;
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    static int cipher[MAXLEN];
    int clen = tridigital_encrypt(plain, n, &g, cipher);

    for (int i = 0; i < clen; i++) putchar('0' + cipher[i]);
    putchar('\n');

    for (int i = 0; i < n; i++) fputc(index_to_char(plain[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

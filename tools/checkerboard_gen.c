// Checkerboard cipher generator (test-data tool, not part of the solver).
//
// Reads a plaintext (first line of a file, or stdin), keeps A..Z, merges J into I (the ACA
// 25-letter convention), builds the keyed 5x5 square from a keyword along the ACA spiral route,
// and enciphers each plaintext letter to a (row label, column label) digraph. One keyword per
// axis is the SIMPLE case; a second keyword per axis (-rows2 / -cols2) is the COMPLEX case, whose
// 2x2 label choice is randomized per position (fixed seed -> reproducible). It links the REAL
// cipher code (checkerboard.c + bifid.c + utils.c), so the generator and the solver can never
// drift in convention.
//
//   make checkerboard_gen
//   ./tools/checkerboard_gen plaintext.txt KNIGHTS BLACK WHITE >cipher.txt 2>solution.txt
//   ./tools/checkerboard_gen - KNIGHTS HORSE GRAYS -rows2 BLACK -cols2 WHITE >c.txt 2>sol.txt
//
// argv: <plaintext|-> <square-keyword> <row-keyword> <col-keyword>
//       [-rows2 <kw>] [-cols2 <kw>] [-route rowmajor|spiral] [-seed N]
// stdout: the ciphertext (one line, bare A..Z, two letters per plaintext letter)
// stderr: the cleaned plaintext (the solution: bare A..Z, J->I -- what the solver recovers)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "colossus.h"

#define MAXLEN (1 << 20)

// Map an uppercase letter to its 25-letter (J->I) alphabet index; -1 otherwise.
static int letter_to_index(int c) {
    c = toupper(c);
    if (c == 'J') c = 'I';
    if (c < 'A' || c > 'Z') return -1;
    return g_char_to_idx[c];
}

// Parse a 5-letter keyword into `side` distinct alphabet indices (the labels of a line). Returns 1
// on success. A label letter colliding with another under J->I (e.g. both I and J) is rejected.
static int parse_labels(const char *s, int side, int lbl[]) {
    int n = 0;
    for (int i = 0; s[i] && n < side; i++) {
        int idx = letter_to_index((unsigned char) s[i]);
        if (idx >= 0) lbl[n++] = idx;
    }
    if (n != side) return 0;
    for (int i = 0; i < side; i++)
        for (int j = i + 1; j < side; j++)
            if (lbl[i] == lbl[j]) return 0;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s <plaintext|-> <square-keyword> <row-keyword> <col-keyword> "
                        "[-rows2 <kw>] [-cols2 <kw>] [-route rowmajor|spiral] [-seed N]\n", argv[0]);
        return 1;
    }
    const char *square_kw = argv[2];
    const char *row_kw1 = argv[3], *col_kw1 = argv[4];
    const char *row_kw2 = NULL, *col_kw2 = NULL;
    int route = CB_ROUTE_SPIRAL_CW;
    unsigned seed = 1u;

    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "-rows2") == 0 && i + 1 < argc)      row_kw2 = argv[++i];
        else if (strcmp(argv[i], "-cols2") == 0 && i + 1 < argc) col_kw2 = argv[++i];
        else if (strcmp(argv[i], "-route") == 0 && i + 1 < argc) {
            i++;
            if (strcasecmp(argv[i], "rowmajor") == 0)    route = CB_ROUTE_ROWMAJOR;
            else if (strcasecmp(argv[i], "spiral") == 0) route = CB_ROUTE_SPIRAL_CW;
            else { fprintf(stderr, "route must be rowmajor or spiral\n"); return 1; }
        } else if (strcmp(argv[i], "-seed") == 0 && i + 1 < argc) {
            seed = (unsigned) strtoul(argv[++i], NULL, 10);
        } else {
            fprintf(stderr, "unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    init_alphabet("J");
    if (g_alpha != CHECKERBOARD_GRID) {
        fprintf(stderr, "alphabet is %d letters, need %d\n", g_alpha, CHECKERBOARD_GRID);
        return 1;
    }
    seed_rand(seed);                               // reproducible complex-case label choices
    int side = CHECKERBOARD_SIDE;

    // Row/column labels, flat by line then choice: rowlbl[row*n_row_lbl + k].
    int n_row_lbl = row_kw2 ? 2 : 1, n_col_lbl = col_kw2 ? 2 : 1;
    int rlbl1[CHECKERBOARD_SIDE], rlbl2[CHECKERBOARD_SIDE];
    int clbl1[CHECKERBOARD_SIDE], clbl2[CHECKERBOARD_SIDE];
    if (!parse_labels(row_kw1, side, rlbl1) || !parse_labels(col_kw1, side, clbl1)) {
        fprintf(stderr, "row/col keyword must give %d distinct letters (J->I)\n", side);
        return 1;
    }
    if (row_kw2 && !parse_labels(row_kw2, side, rlbl2)) { fprintf(stderr, "bad -rows2\n"); return 1; }
    if (col_kw2 && !parse_labels(col_kw2, side, clbl2)) { fprintf(stderr, "bad -cols2\n"); return 1; }

    int rowlbl[CHECKERBOARD_SIDE * 2], collbl[CHECKERBOARD_SIDE * 2];
    for (int r = 0; r < side; r++) {
        rowlbl[r * n_row_lbl + 0] = rlbl1[r];
        if (n_row_lbl > 1) rowlbl[r * n_row_lbl + 1] = rlbl2[r];
    }
    for (int c = 0; c < side; c++) {
        collbl[c * n_col_lbl + 0] = clbl1[c];
        if (n_col_lbl > 1) collbl[c * n_col_lbl + 1] = clbl2[c];
    }
    // Within an axis the union of both label keywords must be distinct (10 labels, complex case).
    if (n_row_lbl > 1)
        for (int a = 0; a < side; a++) for (int b = 0; b < side; b++)
            if (rlbl1[a] == rlbl2[b]) { fprintf(stderr, "row labels collide under J->I\n"); return 1; }
    if (n_col_lbl > 1)
        for (int a = 0; a < side; a++) for (int b = 0; b < side; b++)
            if (clbl1[a] == clbl2[b]) { fprintf(stderr, "col labels collide under J->I\n"); return 1; }

    // Read the plaintext (first line of the file, or stdin), to alphabet indices.
    FILE *fp = (strcmp(argv[1], "-") == 0) ? stdin : fopen(argv[1], "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    static int raw[MAXLEN];
    int n = 0, ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        int idx = letter_to_index(ch);
        if (idx >= 0 && n < MAXLEN) raw[n++] = idx;
    }
    if (fp != stdin) fclose(fp);
    if (n == 0) { fprintf(stderr, "empty plaintext\n"); return 1; }

    // Square keyword -> keyed square along the route.
    int kw[256], kwn = 0;
    for (int i = 0; square_kw[i] && kwn < 256; i++) {
        int idx = letter_to_index((unsigned char) square_kw[i]);
        if (idx >= 0) kw[kwn++] = idx;
    }
    int grid[CHECKERBOARD_GRID];
    checkerboard_square_from_keyword(kw, kwn, route, grid, g_alpha);

    static int cipher[2 * MAXLEN];
    checkerboard_encrypt(raw, n, grid, side, rowlbl, n_row_lbl, collbl, n_col_lbl, cipher);

    for (int i = 0; i < 2 * n; i++) fputc(index_to_char(cipher[i]), stdout);
    fputc('\n', stdout);

    for (int i = 0; i < n; i++) fputc(index_to_char(raw[i]), stderr);
    fputc('\n', stderr);
    return 0;
}

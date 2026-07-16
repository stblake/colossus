//
//  Checkerboard cipher primitives (ACA Checkerboard: keyed 5x5 Polybius square, plaintext letter
//  -> (row label, column label) digraph, read from outside the square).
//
//  The square is a permutation of the active n = side*side letter alphabet (the binary forces
//  g_alpha == 25 for the default -type checkerboard: a 5x5 grid of 0..24, J merged into I). grid[p]
//  is the letter at cell p (row p/side, col p%side); the inverse pos[l] is the cell of letter l.
//  The square is built by laying the keyed-alphabet SEQUENCE (keyword + ascending tail, via
//  bifid_grid_from_keyword) into the cells along a route: CB_ROUTE_SPIRAL_CW (the ACA convention:
//  read the keyword+tail clockwise from the top-left inward) or CB_ROUTE_ROWMAJOR (plain fill).
//
//  Each plaintext letter is enciphered to a digraph. Its cell (row, col) is labelled from OUTSIDE:
//      out[2i]   = a row label of `row`
//      out[2i+1] = a column label of `col`
//  SIMPLE case (one label per line): the digraph is unique -- a bijection side^2 codes <-> letters.
//  COMPLEX case (two labels per line): the encipherer picks one of the 2 row labels and one of the
//  2 column labels FREELY per position, so a letter has 2x2 codes. This primitive randomizes the
//  choice with the shared RNG (like Tri-Square's polyphonic encode); a canonical pick would starve
//  the solver's gradient AND collapse the observed label set, defeating simple/complex detection.
//
//  Decryption maps each label letter back to its row/col line via label2row[]/label2col[] (indexed
//  by label letter, giving 0..side-1 or -1), looks the cell up in the square, and counts the LEGAL
//  positions (a digraph with an unknown label -- only possible under a wrong label map -- decrypts
//  to the sentinel letter 0). The label ORDER (which label denotes line 0 vs 1) is not identifiable
//  ciphertext-only: it is absorbed by a row/column permutation of the recovered square, so the
//  solver searches only the square (simple) or the square plus the label PAIRING (complex).
//

#include "colossus.h"

// Clockwise-from-top-left spiral cell visit order: order[k] is the grid cell (0..side*side-1) that
// receives the k-th sequence letter. For the ACA 5x5 example the sequence KNIGHT+tail lands exactly
// as the printed square. Side-generic (the tests exercise 6x6).
void checkerboard_spiral_order(int order[], int side) {
    int top = 0, bottom = side - 1, left = 0, right = side - 1, idx = 0;
    while (top <= bottom && left <= right) {
        for (int c = left; c <= right; c++) order[idx++] = top * side + c;
        top++;
        for (int r = top; r <= bottom; r++) order[idx++] = r * side + right;
        right--;
        if (top <= bottom) {
            for (int c = right; c >= left; c--) order[idx++] = bottom * side + c;
            bottom--;
        }
        if (left <= right) {
            for (int r = bottom; r >= top; r--) order[idx++] = r * side + left;
            left++;
        }
    }
}

// Build an n-cell (n == side*side) square from a keyword: the keyed-alphabet sequence (keyword
// letters in order, duplicates dropped, then the remaining alphabet ascending) laid into the cells
// along `route`. CB_ROUTE_ROWMAJOR fills cells 0,1,2,...; CB_ROUTE_SPIRAL_CW follows the spiral.
void checkerboard_square_from_keyword(const int keyword[], int kwlen, int route, int grid[], int n) {
    int seq[CHECKERBOARD_MAX_GRID];
    bifid_grid_from_keyword(keyword, kwlen, seq, n);       // the keyed-alphabet sequence
    if (route == CB_ROUTE_SPIRAL_CW) {
        int side = 0;
        while (side * side < n) side++;                    // n is a perfect square (side*side)
        int order[CHECKERBOARD_MAX_GRID];
        checkerboard_spiral_order(order, side);
        for (int k = 0; k < n; k++) grid[order[k]] = seq[k];
    } else {
        for (int k = 0; k < n; k++) grid[k] = seq[k];      // row-major
    }
}

// Encipher plen plaintext letters into 2*plen label letters. rowlbl[row*n_row_lbl + k] is row
// `row`'s k-th label; likewise collbl. The complex case (n_*_lbl > 1) picks a label uniformly at
// random per position (shared RNG); the simple case (n_*_lbl == 1) is deterministic.
void checkerboard_encrypt(const int plain[], int plen, const int grid[], int side,
        const int rowlbl[], int n_row_lbl, const int collbl[], int n_col_lbl, int out[]) {
    int ncell = side * side;
    int pos[CHECKERBOARD_MAX_GRID];
    for (int p = 0; p < ncell; p++) pos[grid[p]] = p;      // letter -> cell
    for (int i = 0; i < plen; i++) {
        int cell = pos[plain[i]];
        int row = cell / side, col = cell % side;
        int kr = (n_row_lbl > 1) ? rand_int(0, n_row_lbl) : 0;
        int kc = (n_col_lbl > 1) ? rand_int(0, n_col_lbl) : 0;
        out[2 * i]     = rowlbl[row * n_row_lbl + kr];
        out[2 * i + 1] = collbl[col * n_col_lbl + kc];
    }
}

// Decipher clen label letters (clen/2 digraphs) into clen/2 plaintext letters; returns the number
// of positions that decrypted to a legal cell (the rest get the sentinel letter 0). label2row[l]
// and label2col[l] give the line (0..side-1) for label letter l, or -1 if l is not a label.
int checkerboard_decrypt(const int cipher[], int clen, const int grid[], int side,
        const int label2row[], const int label2col[], int out_letters[]) {
    int n = clen / 2, n_valid = 0;
    for (int i = 0; i < n; i++) {
        int a = cipher[2 * i], b = cipher[2 * i + 1];
        int row = (a >= 0 && a < MAX_ALPHABET_SIZE) ? label2row[a] : -1;
        int col = (b >= 0 && b < MAX_ALPHABET_SIZE) ? label2col[b] : -1;
        if (row >= 0 && row < side && col >= 0 && col < side) {
            out_letters[i] = grid[row * side + col];
            n_valid++;
        } else {
            out_letters[i] = 0;                            // sentinel: scores as gibberish
        }
    }
    return n_valid;
}

// Count the distinct row labels (1st digraph position) and column labels (2nd position). A simple
// axis can never exceed `side` labels, so > side proves the axis is complex.
void checkerboard_detect(const int cipher[], int clen, int *n_row_lbl, int *n_col_lbl) {
    char seen_r[MAX_ALPHABET_SIZE] = {0}, seen_c[MAX_ALPHABET_SIZE] = {0};
    int nr = 0, nc = 0, n = clen / 2;
    for (int i = 0; i < n; i++) {
        int a = cipher[2 * i], b = cipher[2 * i + 1];
        if (a >= 0 && a < MAX_ALPHABET_SIZE && !seen_r[a]) { seen_r[a] = 1; nr++; }
        if (b >= 0 && b < MAX_ALPHABET_SIZE && !seen_c[b]) { seen_c[b] = 1; nc++; }
    }
    *n_row_lbl = nr;
    *n_col_lbl = nc;
}

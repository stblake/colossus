//
//  Phillips cipher primitives (8-square keyed-Polybius monographic substitution).
//
//  A base square is a permutation of the active n = side*side letter alphabet: the binary
//  forces g_alpha == 25 for -type phillips (J merged into I, ACA convention), so the base
//  is a 5x5 grid of 0..24 indices. base[p] is the letter at cell p (row p/side, col p%side).
//
//  From the base, nsq = 2*side - 2 squares are derived (8 for the 5x5). The plaintext is
//  split into blocks of `side` letters; block b is enciphered with square (b mod nsq).
//  A plaintext letter at cell (r, c) of its square is enciphered by the letter diagonally
//  DOWN-RIGHT with wrap: cipher = sq[(r+1)%side][(c+1)%side]. Decryption takes the UP-LEFT
//  neighbour. The overall period is nsq*side (40 for the 5x5). Each square is a bijection,
//  so the whole map is a bijection -- decryption is exact and needs no length padding.
//
//  The derived squares permute either rows or columns of the base, per `variant`:
//
//    PHILLIPS_ROW (the ACA standard, verified cell-for-cell against the ACA worked
//      example and two tutorials): the squares permute ROWS by a reinsertion sequence --
//        squares 0..side-1 : base row 0 reinserted at row positions 0..side-1
//                            ([0,1,2,3,4] -> [1,0,2,3,4] -> [1,2,0,3,4] -> ... -> [1,2,3,4,0]);
//        squares side..nsq-1: then base row 1 reinserted at positions 1..side-2 of the
//                            square-(side-1) ordering ([2,1,3,4,0] -> [2,3,1,4,0] -> [2,3,4,1,0]).
//
//    PHILLIPS_COL: the column-shift dual -- the SAME reinsertion sequence applied to
//      COLUMNS instead of rows.
//
//    PHILLIPS_ROWCOL: a row+column hybrid --
//        squares 0..side-1 : base row 0 reinserted (columns fixed), as in PHILLIPS_ROW;
//        squares side..nsq-1: base column 0 reinserted at positions 1..side-2 (rows fixed).
//
//  No external authority publishes the COL / ROWCOL conventions (they are the CryptoCrack
//  "Column" / "Row-Column" setup names); the definitions above are fixed here and the
//  generator + solver are round-trip tested against them, so they are internally consistent.
//

#include "colossus.h"

// Per-call scratch for the nsq derived squares and their inverse (letter -> cell) maps.
// Kept file-static so the per-iteration decrypt hook does not carry a large frame.
static _Thread_local int g_phil_squares[PHILLIPS_MAX_SQUARES * PHILLIPS_MAX_GRID];
static _Thread_local int g_phil_pos[PHILLIPS_MAX_SQUARES * PHILLIPS_MAX_GRID];

// Fill order[s*side + pos] with the row/col reinsertion sequence shared by ROW and COL
// (and reused by ROWCOL): the index of the base row/column placed at position `pos` of
// derived square `s`, for s in 0..2*side-3.
static void phillips_fill_orders(int side, int order[]) {
    int last = side - 1;
    // Phase 1 (squares 0..side-1): base index 0 reinserted at position s; the others
    // (1..side-1) keep their order around it.
    for (int s = 0; s < side; s++) {
        int idx = 0;
        for (int pos = 0; pos < side; pos++)
            order[s * side + pos] = (pos == s) ? 0 : (1 + idx++);
    }
    // Phase 2 (squares side..2*side-3): base index 1 reinserted at position j (1..side-2)
    // into the phase-1-final list with its leading 1 removed, i.e. [2,3,...,side-1,0].
    for (int j = 1; j <= side - 2; j++) {
        int s = last + j;
        int idx = 0;
        for (int pos = 0; pos < side; pos++) {
            if (pos == j) order[s * side + pos] = 1;
            else order[s * side + pos] = order[last * side + 1 + idx++];
        }
    }
}

// Build the nsq = 2*side-2 derived squares (cell -> letter) into squares[s*n + cell].
void phillips_build_squares(const int base[], int side, int variant, int squares[]) {
    int n = side * side;
    int nsq = 2 * side - 2;
    int order[PHILLIPS_MAX_SQUARES * PHILLIPS_MAX_SIDE];
    phillips_fill_orders(side, order);

    for (int s = 0; s < nsq; s++) {
        int *sq = squares + s * n;
        if (variant == PHILLIPS_COL) {
            // Permute columns by order[s]; rows unchanged.
            for (int c = 0; c < side; c++) {
                int sc = order[s * side + c];
                for (int r = 0; r < side; r++)
                    sq[r * side + c] = base[r * side + sc];
            }
        } else if (variant == PHILLIPS_ROWCOL && s >= side) {
            // Second phase: permute columns by base col 0 reinserted at position j (1..side-2),
            // which is exactly the phase-1 order[j]; rows unchanged.
            int j = s - (side - 1);
            for (int c = 0; c < side; c++) {
                int sc = order[j * side + c];
                for (int r = 0; r < side; r++)
                    sq[r * side + c] = base[r * side + sc];
            }
        } else {
            // PHILLIPS_ROW (any s), or PHILLIPS_ROWCOL first phase (s < side):
            // permute rows by order[s]; columns unchanged.
            for (int r = 0; r < side; r++) {
                int sr = order[s * side + r];
                for (int c = 0; c < side; c++)
                    sq[r * side + c] = base[sr * side + c];
            }
        }
    }
}

// Build the derived squares and their inverse maps into the file-static scratch.
static void phillips_prepare(const int base[], int side, int variant) {
    int n = side * side, nsq = 2 * side - 2;
    phillips_build_squares(base, side, variant, g_phil_squares);
    for (int s = 0; s < nsq; s++) {
        const int *sq = g_phil_squares + s * n;
        int *pos = g_phil_pos + s * n;
        for (int p = 0; p < n; p++) pos[sq[p]] = p;
    }
}

// Encipher plaintext (any length) under the base square: block b (of `side` letters) uses
// square (b mod nsq); each letter -> its down-right neighbour (with wrap) in that square.
void phillips_encrypt(const int plain[], int len, const int base[], int side, int variant, int out[]) {
    int n = side * side, nsq = 2 * side - 2;
    phillips_prepare(base, side, variant);
    for (int i = 0; i < len; i++) {
        int s = (i / side) % nsq;
        const int *sq = g_phil_squares + s * n;
        int cell = g_phil_pos[s * n + plain[i]];
        int r = cell / side, c = cell % side;
        out[i] = sq[((r + 1) % side) * side + (c + 1) % side];
    }
}

// Decipher ciphertext (any length): the inverse of phillips_encrypt -- each letter -> its
// up-left neighbour (with wrap) in the block's square.
void phillips_decrypt(const int cipher[], int len, const int base[], int side, int variant, int out[]) {
    int n = side * side, nsq = 2 * side - 2;
    phillips_prepare(base, side, variant);
    for (int i = 0; i < len; i++) {
        int s = (i / side) % nsq;
        const int *sq = g_phil_squares + s * n;
        int cell = g_phil_pos[s * n + cipher[i]];
        int r = cell / side, c = cell % side;
        out[i] = sq[((r - 1 + side) % side) * side + (c - 1 + side) % side];
    }
}

// Build a base square from a keyword (alphabet indices): keyword letters in order with
// duplicates removed, then the remaining alphabet letters in ascending order. n cells.
void phillips_grid_from_keyword(const int keyword[], int kwlen, int grid[], int n) {
    char used[PHILLIPS_MAX_GRID];
    for (int l = 0; l < n; l++) used[l] = 0;
    int m = 0;
    for (int i = 0; i < kwlen && m < n; i++) {
        int l = keyword[i];
        if (l < 0 || l >= n || used[l]) continue;
        used[l] = 1;
        grid[m++] = l;
    }
    for (int l = 0; l < n && m < n; l++)
        if (!used[l]) { used[l] = 1; grid[m++] = l; }
}

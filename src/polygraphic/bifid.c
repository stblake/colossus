//
//  Bifid cipher primitives (Delastelle fractionation over a keyed Polybius square).
//
//  The square is a permutation of the active n = side*side letter alphabet: the binary
//  forces g_alpha == 25 for the default -type bifid (J merged into I, ACA convention),
//  so the square is a 5x5 grid of 0..24 indices. grid[p] is the letter at cell p (row
//  p/side, col p%side) and the inverse pos[l] is the cell of letter l. Everything here
//  is side-generic so a 6x6 (36-cell) square works as soon as a 36-letter alphabet is
//  active.
//
//  Encryption of a block of L = period plaintext letters:
//    1. each letter -> its cell coords (r_i, c_i);
//    2. lay the coords out as one length-2L stream  r_0..r_{L-1} c_0..c_{L-1};
//    3. re-pair that stream consecutively: ciphertext letter k = grid[stream[2k]*side
//       + stream[2k+1]].
//  Decryption (all the solver needs) reverses it: expand each ciphertext letter to its
//  two coords to rebuild the stream, then the first L entries are the plaintext rows and
//  the last L the plaintext cols, so plaintext_i = grid[stream[i]*side + stream[L+i]].
//  A short final block (L < period) is enciphered/deciphered the same way in place, so
//  bifid -- unlike Playfair -- needs no length padding.
//

#include "colossus.h"

// Single-threaded scratch for the interleaved coordinate stream of the current block
// (length 2*L <= 2*period <= 2*MAX_CIPHER_LENGTH). Kept off the stack so the per-
// iteration decrypt hook does not carry an 80KB frame.
static _Thread_local int g_bifid_stream[2 * MAX_CIPHER_LENGTH];

// Build pos[] (letter -> cell) from grid[] (cell -> letter); n = side*side cells.
void bifid_build_inverse(const int grid[], int pos[], int n) {
    for (int p = 0; p < n; p++) pos[grid[p]] = p;
}

// Encipher plaintext (any length) under the keyed square, in blocks of `period`.
void bifid_encrypt(const int plain[], int len, const int grid[], int side, int period, int out[]) {
    int pos[BIFID_MAX_GRID];
    bifid_build_inverse(grid, pos, side * side);
    int *stream = g_bifid_stream;
    for (int off = 0; off < len; off += period) {
        int L = (off + period <= len) ? period : (len - off);
        // Lay out this block's coords: all rows first, then all columns.
        for (int i = 0; i < L; i++) {
            int cell = pos[plain[off + i]];
            stream[i]     = cell / side;       // row coordinate
            stream[L + i] = cell % side;       // column coordinate
        }
        // Re-pair the stream consecutively into ciphertext cells.
        for (int k = 0; k < L; k++)
            out[off + k] = grid[stream[2 * k] * side + stream[2 * k + 1]];
    }
}

// Decipher ciphertext (any length) into out[]. Inverse of bifid_encrypt.
void bifid_decrypt(const int cipher[], int len, const int grid[], int side, int period, int out[]) {
    int pos[BIFID_MAX_GRID];
    bifid_build_inverse(grid, pos, side * side);
    int *stream = g_bifid_stream;
    for (int off = 0; off < len; off += period) {
        int L = (off + period <= len) ? period : (len - off);
        // Expand each ciphertext letter to its two coords, rebuilding the 2L stream.
        for (int k = 0; k < L; k++) {
            int cell = pos[cipher[off + k]];
            stream[2 * k]     = cell / side;
            stream[2 * k + 1] = cell % side;
        }
        // First L stream entries are the plaintext rows, last L the plaintext cols.
        for (int i = 0; i < L; i++)
            out[off + i] = grid[stream[i] * side + stream[L + i]];
    }
}

// Build a square from a keyword (alphabet indices): the keyword letters in order with
// duplicates removed, then the remaining alphabet letters in ascending order. n cells.
void bifid_grid_from_keyword(const int keyword[], int kwlen, int grid[], int n) {
    char used[BIFID_MAX_GRID];
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

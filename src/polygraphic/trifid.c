//
//  Trifid cipher primitives (Delastelle fractionation over a keyed 3x3x3 cube).
//
//  The 3D generalization of Bifid. The cube is a permutation of the active n = side^3
//  letter alphabet: the binary forces g_alpha == 27 for the default -type trifid (the
//  full A..Z plus a 27th symbol '+'), so the cube is a 3x3x3 grid of 0..26 indices.
//  cube[p] is the letter at cell p, whose three coordinates are
//      c0 = p / (side*side)     ("layer")
//      c1 = (p / side) % side   ("row")
//      c2 = p % side            ("column")
//  each in 0..side-1, and the inverse pos[l] is the cell of letter l. Everything here is
//  side-generic so a smaller (2x2x2 = 8) or larger (4x4x4 = 64) cube also works.
//
//  Encryption of a block of L = period plaintext letters (matching the Wikipedia /
//  Practical-Cryptography worked example, alphabet "FELIXMARDSTBCGHJKNOPQUVWYZ+":
//  AIDET, period 5 -> FMJFV):
//    1. each letter -> its three coords (c0_i, c1_i, c2_i);
//    2. lay the coords out as one length-3L stream, ALL layers, then ALL rows, then ALL
//       columns:  c0_0..c0_{L-1} c1_0..c1_{L-1} c2_0..c2_{L-1};
//    3. re-group that stream into consecutive triples: ciphertext letter k =
//       cube[stream[3k]*side*side + stream[3k+1]*side + stream[3k+2]].
//  Decryption (all the solver needs) reverses it: expand each ciphertext letter to its
//  three coords to rebuild the stream, then the first L entries are the plaintext layers,
//  the next L the plaintext rows, and the last L the plaintext columns, so
//  plaintext_i = cube[stream[i]*side*side + stream[L+i]*side + stream[2L+i]].
//  A short final block (L < period) is enciphered/deciphered the same way in place, so
//  trifid -- like bifid -- needs no length padding.
//

#include "colossus.h"

// Single-threaded scratch for the interleaved coordinate stream of the current block
// (length 3*L <= 3*period <= 3*MAX_CIPHER_LENGTH). Kept off the stack so the per-
// iteration decrypt hook does not carry a 120KB frame.
static _Thread_local int g_trifid_stream[3 * MAX_CIPHER_LENGTH];

// Build pos[] (letter -> cell) from cube[] (cell -> letter); n = side^3 cells.
void trifid_build_inverse(const int cube[], int pos[], int n) {
    for (int p = 0; p < n; p++) pos[cube[p]] = p;
}

// Encipher plaintext (any length) under the keyed cube, in blocks of `period`.
void trifid_encrypt(const int plain[], int len, const int cube[], int side, int period, int out[]) {
    int pos[TRIFID_MAX_CELLS];
    int side2 = side * side;
    trifid_build_inverse(cube, pos, side2 * side);
    int *stream = g_trifid_stream;
    for (int off = 0; off < len; off += period) {
        int L = (off + period <= len) ? period : (len - off);
        // Lay out this block's coords: all layers, then all rows, then all columns.
        for (int i = 0; i < L; i++) {
            int cell = pos[plain[off + i]];
            stream[i]         = cell / side2;          // layer coordinate (c0)
            stream[L + i]     = (cell / side) % side;  // row coordinate   (c1)
            stream[2 * L + i] = cell % side;           // column coordinate (c2)
        }
        // Re-group the stream into consecutive triples indexing new cube cells.
        for (int k = 0; k < L; k++)
            out[off + k] = cube[stream[3 * k] * side2 + stream[3 * k + 1] * side + stream[3 * k + 2]];
    }
}

// Decipher ciphertext (any length) into out[]. Inverse of trifid_encrypt.
void trifid_decrypt(const int cipher[], int len, const int cube[], int side, int period, int out[]) {
    int pos[TRIFID_MAX_CELLS];
    int side2 = side * side;
    trifid_build_inverse(cube, pos, side2 * side);
    int *stream = g_trifid_stream;
    for (int off = 0; off < len; off += period) {
        int L = (off + period <= len) ? period : (len - off);
        // Expand each ciphertext letter to its three coords, rebuilding the 3L stream.
        for (int k = 0; k < L; k++) {
            int cell = pos[cipher[off + k]];
            stream[3 * k]     = cell / side2;
            stream[3 * k + 1] = (cell / side) % side;
            stream[3 * k + 2] = cell % side;
        }
        // First L stream entries are the plaintext layers, next L the rows, last L cols.
        for (int i = 0; i < L; i++)
            out[off + i] = cube[stream[i] * side2 + stream[L + i] * side + stream[2 * L + i]];
    }
}

// Build a cube from a keyword (alphabet indices): the keyword letters in order with
// duplicates removed, then the remaining alphabet letters in ascending order. n cells.
void trifid_cube_from_keyword(const int keyword[], int kwlen, int cube[], int n) {
    char used[TRIFID_MAX_CELLS];
    for (int l = 0; l < n; l++) used[l] = 0;
    int m = 0;
    for (int i = 0; i < kwlen && m < n; i++) {
        int l = keyword[i];
        if (l < 0 || l >= n || used[l]) continue;
        used[l] = 1;
        cube[m++] = l;
    }
    for (int l = 0; l < n && m < n; l++)
        if (!used[l]) { used[l] = 1; cube[m++] = l; }
}

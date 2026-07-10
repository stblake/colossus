//
//  ADFGVX / ADFGX cipher primitives (fractionation over a keyed Polybius square
//  composed with a keyed columnar transposition).
//
//  The classic WWI German field cipher. Two layers:
//    1. SUBSTITUTION/FRACTIONATION. A keyed side x side Polybius square holds the
//       n = side*side plaintext symbols (ADFGX: side 5, 25 letters with J->I;
//       ADFGVX: side 6, 36 symbols A..Z + 0..9). Each plaintext symbol is replaced
//       by its (row, col) cell coordinates, each a label drawn from a small set --
//       {A,D,F,G,X} for ADFGX, {A,D,F,G,V,X} for ADFGVX. This DOUBLES the length:
//       N plaintext symbols become a 2N stream of coordinates (each 0..side-1).
//    2. TRANSPOSITION. The 2N coordinate stream is written row-major into a grid of
//       K columns (K = the transposition keyword length) and the columns are read
//       off in keyed order `order[0..K-1]`, each top-to-bottom (COL_READ_TB) or
//       bottom-to-top (COL_READ_BT). The result is the ciphertext, a 2N string over
//       the label alphabet.
//
//  This module works in COORDINATE space: the ciphertext is carried as a length-2N
//  array of coordinates 0..side-1 (the solver maps the label characters to
//  coordinates once up front). The square is a permutation of 0..n-1 carried in
//  0..n-1 indices, exactly like Bifid -- so the keyed-square build and the cell
//  inverse are shared with bifid.c (bifid_grid_from_keyword / bifid_build_inverse).
//
//  Decryption (what the solver needs) inverts both layers: undo the columnar with
//  the shared decrypt_columnar() primitive to recover the row-major coordinate
//  stream, then pair consecutive coordinates and look the cell up in the square:
//    plain[i] = square[stream[2i]*side + stream[2i+1]].
//  Because the fractionation always emits exactly two coordinates per symbol, the
//  coordinate length 2N is always even and the columnar grid may be ragged
//  (2N % K != 0), which decrypt_columnar handles.
//

#include "colossus.h"

// The label alphabets. The coordinate value c (0..side-1) is displayed as labels[c];
// conversely the solver maps each ciphertext label character back to its coordinate.
// Shared by the solver and the test-data generator so they cannot drift.
const char ADFGX_LABELS[]  = "ADFGX";    // side 5
const char ADFGVX_LABELS[] = "ADFGVX";   // side 6

const char *adfgvx_labels(int side) {
    return (side >= 6) ? ADFGVX_LABELS : ADFGX_LABELS;
}

// Single-threaded scratch for the row-major coordinate stream (length 2N <=
// MAX_CIPHER_LENGTH). Kept off the stack so the per-iteration decrypt hook does not
// carry a large frame.
static _Thread_local int g_adfgvx_stream[2 * MAX_CIPHER_LENGTH];

// Columnar transposition ENCRYPT (the inverse of decrypt_columnar): write `in`
// row-major into K columns (the leftmost len % K columns one cell taller when the
// grid is ragged), then read the columns in read order `order[0..K-1]`, each
// top-to-bottom (COL_READ_TB) or bottom-to-top (COL_READ_BT). out[] must not alias
// in[]. Used only by adfgvx_encrypt; the solver inverts via the exposed
// decrypt_columnar.
static void adfgvx_columnar_encrypt(const int in[], int len, int K, const int order[],
                                    int dir, int out[]) {
    if (K <= 1 || K > len) {            // degenerate: a single column is the identity
        for (int i = 0; i < len; i++) out[i] = in[i];
        return;
    }
    int R = (len + K - 1) / K;          // number of rows (ceiling)
    int rem = len % K;                  // tall columns are 0..rem-1 (all K if rem==0)
    int pos = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j];
        int h = (rem == 0 || c < rem) ? R : R - 1;   // height of grid column c
        if (dir == COL_READ_BT)
            for (int r = h - 1; r >= 0; r--) out[pos++] = in[r * K + c];
        else
            for (int r = 0; r < h; r++)     out[pos++] = in[r * K + c];
    }
}

// Encipher n plaintext symbols (alphabet indices 0..side*side-1) into a 2n stream of
// ciphertext coordinates (each 0..side-1). square[p] is the symbol at cell p (row
// p/side, col p%side); side*side cells. out[] must not alias plain[].
void adfgvx_encrypt(const int plain[], int n, const int square[], int side, int K,
                    const int order[], int dir, int out[]) {
    int pos[SQUARE_MAX_GRID];                       // symbol -> cell
    bifid_build_inverse(square, pos, side * side);   // shared with Bifid
    int *frac = g_adfgvx_stream;
    for (int i = 0; i < n; i++) {
        int cell = pos[plain[i]];
        frac[2 * i]     = cell / side;              // row label
        frac[2 * i + 1] = cell % side;              // column label
    }
    adfgvx_columnar_encrypt(frac, 2 * n, K, order, dir, out);
}

// Decipher a 2N coordinate stream (each 0..side-1) back into N plaintext symbols.
// Inverse of adfgvx_encrypt: undo the columnar (decrypt_columnar) to recover the
// row-major coordinate stream, then pair consecutive coordinates and look up the
// square. len2 == 2N must be even; out[] (length N) must not alias cipher[].
void adfgvx_decrypt(const int cipher[], int len2, const int square[], int side, int K,
                    const int order[], int dir, int out[]) {
    int *frac = g_adfgvx_stream;
    // decrypt_columnar reads cipher/order and writes frac only; cast away const.
    decrypt_columnar((int *) cipher, len2, K, (int *) order, dir, frac);
    int n = len2 / 2;
    for (int i = 0; i < n; i++)
        out[i] = square[frac[2 * i] * side + frac[2 * i + 1]];
}

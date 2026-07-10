//
//  CM Bifid (Conjugated Matrix Bifid) cipher primitives.
//
//  CM Bifid is plain Bifid (Delastelle fractionation, block size = period) with TWO keyed
//  Polybius squares instead of one: the plaintext is fractionated into (row,col) coords
//  using square 1, and after the rows-then-cols reshape and re-pairing each coordinate
//  pair is mapped to a ciphertext letter using a DIFFERENT square 2. When sq1 == sq2 this
//  is exactly Bifid -- a property the unit tests pin.
//
//  Both squares are permutations of the active n = side*side letter alphabet (the binary
//  forces g_alpha == 25 for -type cm-bifid, J merged into I), so each is a 5x5 grid of
//  0..24 indices. sq[p] is the letter at cell p (row p/side, col p%side); the inverse
//  pos[l] is the cell of letter l. Everything is side-generic so a 6x6 (36-cell) pair
//  works once a 36-letter alphabet is active. Reuses bifid_build_inverse (bifid.c).
//
//  Encryption of a block of L = period plaintext letters:
//    1. each plaintext letter -> its cell coords in SQ1: (r_i, c_i);
//    2. lay the coords out as one length-2L stream  r_0..r_{L-1} c_0..c_{L-1};
//    3. re-pair that stream consecutively and look the pairs up in SQ2:
//       ciphertext letter k = sq2[stream[2k]*side + stream[2k+1]].
//  Decryption reverses it: expand each ciphertext letter to its two coords via SQ2 to
//  rebuild the stream, then the first L entries are the plaintext rows and the last L the
//  plaintext cols, recombined through SQ1: plaintext_i = sq1[stream[i]*side + stream[L+i]].
//  A short final block (L < period) is enciphered/deciphered the same way in place, so
//  CM Bifid -- like Bifid -- needs no length padding.
//

#include "colossus.h"

// Single-threaded scratch for the interleaved coordinate stream of the current block
// (length 2*L <= 2*period <= 2*MAX_CIPHER_LENGTH). Kept off the stack so the per-
// iteration decrypt hook does not carry an 80KB frame.
static _Thread_local int g_cm_bifid_stream[2 * MAX_CIPHER_LENGTH];

// Encipher plaintext (any length): fractionate with sq1, recombine the pairs through sq2.
void cm_bifid_encrypt(const int plain[], int len, const int sq1[], const int sq2[],
                      int side, int period, int out[]) {
    int pos1[BIFID_MAX_GRID];
    bifid_build_inverse(sq1, pos1, side * side);
    int *stream = g_cm_bifid_stream;
    for (int off = 0; off < len; off += period) {
        int L = (off + period <= len) ? period : (len - off);
        // Lay out this block's coords from SQ1: all rows first, then all columns.
        for (int i = 0; i < L; i++) {
            int cell = pos1[plain[off + i]];
            stream[i]     = cell / side;       // row coordinate
            stream[L + i] = cell % side;       // column coordinate
        }
        // Re-pair the stream consecutively and look the pairs up in SQ2.
        for (int k = 0; k < L; k++)
            out[off + k] = sq2[stream[2 * k] * side + stream[2 * k + 1]];
    }
}

// Decipher ciphertext (any length) into out[]. Inverse of cm_bifid_encrypt.
void cm_bifid_decrypt(const int cipher[], int len, const int sq1[], const int sq2[],
                      int side, int period, int out[]) {
    int pos2[BIFID_MAX_GRID];
    bifid_build_inverse(sq2, pos2, side * side);
    int *stream = g_cm_bifid_stream;
    for (int off = 0; off < len; off += period) {
        int L = (off + period <= len) ? period : (len - off);
        // Expand each ciphertext letter to its two coords via SQ2, rebuilding the 2L stream.
        for (int k = 0; k < L; k++) {
            int cell = pos2[cipher[off + k]];
            stream[2 * k]     = cell / side;
            stream[2 * k + 1] = cell % side;
        }
        // First L stream entries are the plaintext rows, last L the plaintext cols;
        // recombine through SQ1.
        for (int i = 0; i < L; i++)
            out[off + i] = sq1[stream[i] * side + stream[L + i]];
    }
}

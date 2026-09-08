// =====================================================================
//  Chaocipher primitive (TYPE chaocipher)
// =====================================================================
//
// Chaocipher (John F. Byrne, 1918) is a self-modifying substitution driven by two
// 26-letter alphabets ("disks"):
//
//   * LEFT  alphabet = the CIPHERTEXT disk.
//   * RIGHT alphabet = the PLAINTEXT  disk.
//
// Each alphabet is a permutation of A..Z. To encipher a plaintext letter you find it in
// the RIGHT alphabet; the letter at the SAME position in the LEFT alphabet is the
// ciphertext letter. Both alphabets are then PERMUTED, so the substitution map drifts
// after every character. The KEY is the pair of STARTING alphabets.
//
// Positions here are 0-indexed: the "zenith" is index 0 and the "nadir" is index 13 (the
// paper's positions 1 and 14). The permutations after enciphering one letter are:
//
//   LEFT  (ciphertext): rotate so the just-used ciphertext letter is at index 0; then
//     remove the element at index 1 and reinsert it at index 13 (rotate block [1..13] by 1).
//   RIGHT (plaintext) : rotate so the just-used plaintext letter is at index 0; rotate one
//     MORE (its successor to index 0); then remove the element at index 2 and reinsert it
//     at index 13 (rotate block [2..13] by 1).
//
// Deciphering is identical with the disks' roles swapped: find the ciphertext letter in the
// LEFT alphabet, read the plaintext letter from the RIGHT alphabet, then permute both disks
// exactly as above. A letter CAN encipher to itself (this is not a reciprocal cipher).
//
// Rotating BOTH starting alphabets by the same amount yields an identical cipher (after the
// first step the states converge), so the effective keyspace is 26! * 25!, not (26!)^2.
//
// Source: Moshe Rubin, "Chaocipher Revealed: The Algorithm" (2010); the worked example
// there is pinned as a KAT in tests/test_chaocipher.c. Letters are 0..25 (A=0).

#ifndef CHAOCIPHER_H
#define CHAOCIPHER_H

#include "colossus.h"

#define CHAO_N       26   // alphabet size (Chaocipher is always the full 26-letter A..Z set)
#define CHAO_ZENITH   0   // reading position (paper position 1)
#define CHAO_NADIR   13   // opposite position (paper position 14)

// Advance one disk by one step: cyclically rotate a[0..25] so a[pos] moves to the zenith
// (index 0), then apply the disk's permutation. These are the low-level disk steppers used
// by encrypt/decrypt (and by the Chaocipher solver's aggregate-displacement-error fitness,
// which advances each disk independently by its own known letter stream). LEFT is the
// ciphertext disk's permutation, RIGHT the plaintext disk's.
void chaocipher_step_left(int a[CHAO_N], int pos);
void chaocipher_step_right(int a[CHAO_N], int pos);

// Encipher plain[0..n-1] with the two STARTING alphabets left0[26]/right0[26] (each a
// permutation of 0..25); write ciphertext to out[0..n-1]. out[] must not alias left0/right0.
void chaocipher_encrypt(const int plain[], int n,
                        const int left0[CHAO_N], const int right0[CHAO_N], int out[]);

// Decipher cipher[0..n-1] with the same starting alphabets; write plaintext to out[0..n-1].
void chaocipher_decrypt(const int cipher[], int n,
                        const int left0[CHAO_N], const int right0[CHAO_N], int out[]);

#endif

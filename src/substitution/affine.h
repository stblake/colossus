#ifndef AFFINE_H
#define AFFINE_H
#include "colossus.h"

// =====================================================================
//  Affine cipher primitives (TYPE affine)
// =====================================================================
//
// A monoalphabetic substitution CT = (a*PT + b) mod 26 with gcd(a,26) = 1 (a "decimation"
// a composed with a Caesar shift b). Decode is PT = a_inv*(CT - b) mod 26. There are only
// 12 valid multipliers a in {1,3,5,7,9,11,15,17,19,21,23,25} (odd, not divisible by 13) and
// 26 shifts b -> 312 keys total, so the solver enumerates them exhaustively. a = 1 degenerates
// to a Caesar shift; a = 1, b = 0 is the identity.

#define AFFINE_MOD    26   // the classical 26-letter modulus (the phrase alphabet is fixed A..Z)
#define AFFINE_N_MULT 12   // count of a in [1,25] coprime to 26

// The 12 coprime multipliers and their inverses mod 26 (index-aligned: inverse of
// affine_multipliers[i] is affine_inverses[i]).
extern const int affine_multipliers[AFFINE_N_MULT];
extern const int affine_inverses[AFFINE_N_MULT];

// Modular inverse of a mod 26 (a coprime to 26). Returns -1 if a is not invertible.
int affine_mod_inverse(int a);

// Build the 26-entry encode / decode maps: enc[p] = (a*p + b) mod 26, dec[c] = a_inv*(c-b) mod 26.
void affine_encrypt_map(int a, int b, int enc[AFFINE_MOD]);
void affine_decrypt_map(int a, int b, int dec[AFFINE_MOD]);

// Encipher / decipher a letter-index stream in[0..len-1] (values 0..25). Non-letter sentinels
// (values < 0) pass through unchanged.
void affine_encrypt(const int in[], int len, int a, int b, int out[]);
void affine_decrypt(const int in[], int len, int a, int b, int out[]);
#endif

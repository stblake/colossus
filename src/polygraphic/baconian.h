#ifndef BACONIAN_H
#define BACONIAN_H
#include "colossus.h"

// =====================================================================
//  Baconian cipher primitives (biliteral 5-symbol substitution)
// =====================================================================
//
//  The Baconian cipher (ACA) replaces each plaintext letter with a group of five
//  a/b symbols through the FIXED 24-letter biliteral table (I=J and U=V share a
//  code), read in alphabetical order:
//
//      A=aaaaa B=aaaab C=aaaba D=aaabb E=aabaa F=aabab G=aabba H=aabbb I/J=abaaa
//      K=abaab L=ababa M=ababb N=abbaa O=abbab P=abbba Q=abbbb R=baaaa S=baaab
//      T=baaba U/V=baabb W=babaa X=babab Y=babba Z=babbb
//
//  With a=0 and b=1 the 5-bit code (most-significant symbol first) is exactly the
//  integer 0..23; the eight patterns 24..31 (11000..11111) are IMPOSSIBLE, which is
//  a strong VALIDITY signal the solver exploits (the Morse-validity analogue).
//
//  These primitives are the PURE biliteral layer only: letter <-> 5-bit a/b group.
//  The concealment (which cover-text letters/words stand for a vs b) is the solver's
//  job and lives in baconian_solver.c -- it produces the bit[] stream this decodes.

// The two biliteral stream symbols.
enum { BAC_A = 0, BAC_B = 1 };

// Build the fixed value<->letter tables (idempotent). Called by encrypt/decode.
void baconian_init(void);

// Encipher plaintext pt[0..n-1] (each 0..25; other values skipped) into the a/b bit
// stream out_bits[] -- five BAC_A/BAC_B symbols per letter, most-significant first.
// J folds to I's code and V to U's code (the 24-letter table). Returns the bit count
// (5 * #letters emitted). out_bits[] must hold up to 5*n entries. Used by the
// generator + unit tests (the solver only decodes).
int baconian_encrypt(const int pt[], int n, int out_bits[]);

// Decode the a/b bit stream bits[0..nbits-1] (each BAC_A/BAC_B) in fixed groups of
// five: each group is the 5-bit value v (0..31); v <= 23 -> its letter (a valid
// token), v >= 24 -> `filler` (an invalid token, an impossible biliteral code). A
// trailing partial group of < 5 bits is dropped. Returns the number of letters
// written (= nbits/5); *n_tokens / *n_valid receive the total-group / valid-group
// counts (either pointer may be NULL). out[] must hold up to nbits/5 letters. On the
// true bit stream this recovers exactly the plaintext (all groups valid).
int baconian_decode(const int bits[], int nbits, int out[],
                    int filler, int *n_tokens, int *n_valid);
#endif

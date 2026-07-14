#include "colossus.h"
#include <ctype.h>

// =====================================================================
//  Aristocrat / Patristocrat primitives (TYPES aristocrat, patristocrat)
// =====================================================================
//
// The ACA Aristocrat (see cryptogram.org .../cipher-types) is a simple MONOALPHABETIC
// substitution: one fixed 26-letter bijection cmap[] sends each plaintext letter p to the
// ciphertext letter cmap[p]. The defining feature is that the WORD DIVISIONS are PRESERVED --
// the ciphertext keeps its spaces, so a solver can read off word lengths. The Patristocrat is
// the identical cipher with the word divisions removed (ciphertext written in 5-letter groups).
// The enciphering math is the same for both; only the presentation differs, so both share the
// single bijection here.
//
// ACA aristocrats build cmap[] from a KEYED alphabet (keyword letters, duplicates dropped, then
// the rest ascending) in one of four "keying" arrangements (K1..K4). The solver does NOT need
// the keying -- it recovers the free bijection directly -- but aristocrat_build_map() lets the
// generator and the unit tests plant realistic ACA vectors. The math works on already-parsed
// LETTER streams (indices 0..alpha-1); spaces / non-letters are handled by the caller.
//
// Worked example (K2, keyword KRYPTOS): the keyed ciphertext alphabet is
//   KRYPTOSABCDEFGHIJLMNQUVWXZ, straight plaintext A..Z, so A->K B->R ... and
//   "THE QUICK BROWN FOX" enciphers to "NAT JQBYD RLHVG OHW".

// Build the substitution map cmap[] (cmap[pt] = ct). The plaintext row pt_alpha[] and the
// ciphertext row ct_alpha[] are two alphabets written one above the other; enciphering sends
// pt_alpha[i] -> ct_alpha[i], i.e. cmap[pt_alpha[i]] = ct_alpha[i].
void aristocrat_build_map(int keying, const char *kw1, const char *kw2, int shift, int cmap[]) {
    int pt_alpha[ALPHABET_SIZE], ct_alpha[ALPHABET_SIZE];
    int keyed1[ALPHABET_SIZE], keyed2[ALPHABET_SIZE];
    make_keyed_alphabet((char *) (kw1 ? kw1 : ""), keyed1);
    make_keyed_alphabet((char *) (kw2 ? kw2 : ""), keyed2);

    for (int i = 0; i < ALPHABET_SIZE; i++) { pt_alpha[i] = i; ct_alpha[i] = i; }

    switch (keying) {
        case ARIST_K1:                          // keyed plaintext, straight ciphertext
            for (int i = 0; i < ALPHABET_SIZE; i++) pt_alpha[i] = keyed1[i];
            break;
        case ARIST_K3:                          // both keyed with the same keyword, CT row rotated
            for (int i = 0; i < ALPHABET_SIZE; i++) {
                pt_alpha[i] = keyed1[i];
                ct_alpha[i] = keyed1[((i + shift) % ALPHABET_SIZE + ALPHABET_SIZE) % ALPHABET_SIZE];
            }
            break;
        case ARIST_K4:                          // both keyed, different keywords
            for (int i = 0; i < ALPHABET_SIZE; i++) { pt_alpha[i] = keyed1[i]; ct_alpha[i] = keyed2[i]; }
            break;
        case ARIST_K2:                          // straight plaintext, keyed ciphertext (default)
        default:
            for (int i = 0; i < ALPHABET_SIZE; i++) ct_alpha[i] = keyed1[i];
            break;
    }

    for (int i = 0; i < ALPHABET_SIZE; i++) cmap[pt_alpha[i]] = ct_alpha[i];
}

void aristocrat_invert(const int cmap[], int cinv[], int alpha) {
    for (int p = 0; p < alpha; p++) cinv[cmap[p]] = p;
}

void aristocrat_apply(const int in[], int len, const int map[], int out[]) {
    for (int i = 0; i < len; i++) out[i] = map[in[i]];
}

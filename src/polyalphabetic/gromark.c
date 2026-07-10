//
// Gromark / Periodic Gromark cipher primitives
//

#include "colossus.h"

/*
   Gromark ("GROnsfeld with Mixed Alphabet and Running Key", ACA, DUMBO 1969)

   A keyed 26-letter substitution composed with a chain-addition running key.
   Let sigma be the mixed (keyed) CIPHER alphabet -- a permutation of A..Z -- and let
   the plaintext alphabet be the straight A..Z. A P-digit primer seeds a numeric
   running key by chain addition:

       d[i] = primer[i]                              for i < P
       d[i] = (d[i-P] + d[i-P+1]) mod 10             for i >= P

   (the "add successive pairs of digits, dropping tens" rule: 1st+2nd give the
   (P+1)th, 2nd+3rd give the (P+2)th, ...). One digit per plaintext letter.

   Encipherment locates the plaintext letter in the straight alphabet, counts d[i]
   places to the right, and reads the letter below in the mixed cipher alphabet:

       BASIC GROMARK    C[i] = sigma[(p[i] + d[i]) mod 26]

   The standard ACA primer is 5 digits. Decryption inverts via sigma's inverse:
   p[i] = (sigma_inv[C[i]] - d[i]) mod 26.

   Periodic Gromark (DUMBO 1973) adds a per-GROUP offset. The plaintext is split into
   consecutive groups of P letters (P = keyword length); group g (cycling mod P) is
   headed by the g-th keyword letter, whose offset is that letter's position in the
   mixed alphabet (offset[g] = sigma_inv[keyword[g]]). Within a group the running key
   still advances per position:

       PERIODIC GROMARK C[i] = sigma[(p[i] + d[i] + offset[(i/P) mod P]) mod 26]

   The primer for Periodic Gromark is the alphabetical ranks of the keyword letters
   (a permutation of 1..P), but the primitive takes the primer/offsets explicitly --
   the solver searches them, the generator derives them from a keyword.

   The mixed alphabet is built by a K2M transposition block: the simple keyed
   alphabet (keyword's distinct letters, then the rest A..Z) written row-major at
   width W = #distinct keyword letters, read off by columns in the order of their
   heading letter (i.e. ascending alphabetical value of the keyword letters). Only the
   generator and the unit tests build sigma from a keyword; the solver hill-climbs it
   as a free permutation.

   Everything runs on the full 26-letter alphabet (mod base 26).
*/

// Chain-addition running key: d[0..n-1] from a primer of `primer_len` digits (0..9).
// d[i] = primer[i] (i < primer_len), else (d[i-P] + d[i-P+1]) mod 10. A 1-digit primer
// has no "pair" to add, so it degenerates to the constant primer digit.
void gromark_chain_key(const int primer[], int primer_len, int n, int out_digits[]) {
    int P = primer_len;
    for (int i = 0; i < n; i++) {
        if (i < P) {
            out_digits[i] = primer[i] % 10;
        } else if (P >= 2) {
            out_digits[i] = (out_digits[i - P] + out_digits[i - P + 1]) % 10;
        } else {
            out_digits[i] = out_digits[i - 1];   // P == 1: constant key
        }
    }
}

// Number of DISTINCT letters in a keyword (the K2M transposition-block width).
static int gromark_keyword_width(const char *keyword) {
    int seen[ALPHABET_SIZE] = {0}, w = 0;
    for (int i = 0; keyword[i]; i++) {
        int idx = g_char_to_idx[(unsigned char) toupper((unsigned char) keyword[i]) & 127];
        if (idx >= 0 && idx < ALPHABET_SIZE && !seen[idx]) { seen[idx] = 1; w++; }
    }
    return w;
}

// Build the K2M mixed cipher alphabet (a permutation of 0..25) from a keyword. The
// simple keyed alphabet is laid row-major at width W = #distinct keyword letters;
// the columns are read off in ascending order of their heading letter's value. For
// keyword "ENIGMA" this reproduces the ACA sigma AJRXEBKSYGFPVIDOUMHQWNCLTZ.
void gromark_mixed_alphabet(const char *keyword, int sigma[]) {
    int keyed[ALPHABET_SIZE];
    make_keyed_alphabet((char *) keyword, keyed);     // distinct keyword letters, then the rest
    int W = gromark_keyword_width(keyword);
    if (W < 1) W = 1;

    // Column read order: the W column indices sorted by their heading letter keyed[col]
    // (ascending alphabetical value). Insertion sort (W <= 26).
    int order[ALPHABET_SIZE];
    for (int c = 0; c < W; c++) order[c] = c;
    for (int a = 1; a < W; a++) {
        int v = order[a], b = a - 1;
        while (b >= 0 && keyed[order[b]] > keyed[v]) { order[b + 1] = order[b]; b--; }
        order[b + 1] = v;
    }

    int out = 0;
    for (int oc = 0; oc < W; oc++) {
        int col = order[oc];
        for (int row = 0; row * W + col < ALPHABET_SIZE; row++)
            sigma[out++] = keyed[row * W + col];
    }
}

// Inverse permutation of a 26-letter alphabet: inv[sigma[i]] = i.
static void gromark_inverse(const int sigma[], int inv[]) {
    for (int i = 0; i < ALPHABET_SIZE; i++) inv[sigma[i]] = i;
}

// Build the full Periodic Gromark key from a keyword given as P DISTINCT letter indices
// kw[0..P-1]: the K2M mixed alphabet sigma, the primer (the alphabetical ranks of the keyword
// letters, 1..P), and the per-group offsets (the keyword letters' positions in sigma). This is
// the index-space form of the keyword->everything derivation the Periodic Gromark solver climbs.
void gromark_build_from_keyword_idx(const int kw[], int P, int sigma[],
                                    int primer[], int offsets[]) {
    // Keyed alphabet: the P keyword letters, then the remaining letters in alphabetical order.
    int seen[ALPHABET_SIZE] = {0}, keyed[ALPHABET_SIZE], pos = 0;
    for (int i = 0; i < P; i++) { keyed[pos++] = kw[i]; seen[kw[i]] = 1; }
    for (int i = 0; i < ALPHABET_SIZE; i++) if (!seen[i]) keyed[pos++] = i;

    // K2M transposition block of width W = P: read columns in heading-letter (kw) order.
    int order[ALPHABET_SIZE];
    for (int c = 0; c < P; c++) order[c] = c;
    for (int a = 1; a < P; a++) {
        int v = order[a], b = a - 1;
        while (b >= 0 && keyed[order[b]] > keyed[v]) { order[b + 1] = order[b]; b--; }
        order[b + 1] = v;
    }
    int out = 0;
    for (int oc = 0; oc < P; oc++) {
        int col = order[oc];
        for (int row = 0; row * P + col < ALPHABET_SIZE; row++) sigma[out++] = keyed[row * P + col];
    }

    int sinv[ALPHABET_SIZE];
    gromark_inverse(sigma, sinv);
    for (int g = 0; g < P; g++) {
        int rank = 1;
        for (int h = 0; h < P; h++) if (kw[h] < kw[g]) rank++;
        primer[g] = rank;
        offsets[g] = sinv[kw[g]];
    }
}

// --- basic Gromark ------------------------------------------------------------

void gromark_encrypt(const int plain[], int len, const int sigma[],
                     const int primer[], int primer_len, int out[]) {
    static _Thread_local int d[MAX_CIPHER_LENGTH];
    gromark_chain_key(primer, primer_len, len, d);
    for (int i = 0; i < len; i++)
        out[i] = sigma[(plain[i] + d[i]) % ALPHABET_SIZE];
}

void gromark_decrypt(const int cipher[], int len, const int sigma[],
                     const int primer[], int primer_len, int out[]) {
    static _Thread_local int d[MAX_CIPHER_LENGTH];
    int inv[ALPHABET_SIZE];
    gromark_chain_key(primer, primer_len, len, d);
    gromark_inverse(sigma, inv);
    for (int i = 0; i < len; i++)
        out[i] = (inv[cipher[i]] - d[i] % ALPHABET_SIZE + ALPHABET_SIZE) % ALPHABET_SIZE;
}

// --- Periodic Gromark ---------------------------------------------------------

void gromark_periodic_encrypt(const int plain[], int len, const int sigma[],
                     const int primer[], int period, const int offsets[], int out[]) {
    static _Thread_local int d[MAX_CIPHER_LENGTH];
    gromark_chain_key(primer, period, len, d);
    for (int i = 0; i < len; i++) {
        int g = (i / period) % period;
        out[i] = sigma[(plain[i] + d[i] + offsets[g]) % ALPHABET_SIZE];
    }
}

void gromark_periodic_decrypt(const int cipher[], int len, const int sigma[],
                     const int primer[], int period, const int offsets[], int out[]) {
    static _Thread_local int d[MAX_CIPHER_LENGTH];
    int inv[ALPHABET_SIZE];
    gromark_chain_key(primer, period, len, d);
    gromark_inverse(sigma, inv);
    for (int i = 0; i < len; i++) {
        int g = (i / period) % period;
        out[i] = ((inv[cipher[i]] - d[i] - offsets[g]) % ALPHABET_SIZE + 2 * ALPHABET_SIZE)
                 % ALPHABET_SIZE;
    }
}

// Hot-path decrypt for the solver: precomputed running key d[] and inverse alphabet
// sigma_inv[]; offsets == NULL => basic Gromark, else periodic with the given period.
void gromark_decrypt_core(const int cipher[], int len, const int sigma_inv[],
                          const int d[], const int offsets[], int period, int out[]) {
    if (offsets == NULL) {
        for (int i = 0; i < len; i++)
            out[i] = (sigma_inv[cipher[i]] - d[i] + ALPHABET_SIZE) % ALPHABET_SIZE;
    } else {
        for (int i = 0; i < len; i++) {
            int g = (i / period) % period;
            out[i] = ((sigma_inv[cipher[i]] - d[i] - offsets[g]) % ALPHABET_SIZE
                      + 2 * ALPHABET_SIZE) % ALPHABET_SIZE;
        }
    }
}

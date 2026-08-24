//
//  Baconian cipher primitives -- see baconian.h for the full description.
//
//  The pure biliteral layer: a plaintext letter <-> a 5-symbol a/b group through the
//  FIXED 24-letter table (I=J, U=V). The concealment (cover-text -> a/b bits) is the
//  solver's job (baconian_solver.c); these primitives only map letters <-> bits and
//  serve the generator + unit tests.
//

#include "baconian.h"

// value 0..23 -> A..Z index of the 24-letter Baconian order (I=J, U=V merged): the
// alphabet skipping J and V. value 24..31 are impossible (no biliteral code).
static int  g_bac_val_to_letter[24];
// A..Z index (0..25) -> its 5-bit value 0..23 (J shares I's, V shares U's).
static int  g_bac_letter_to_val[26];
static bool g_bac_init = false;

void baconian_init(void) {
    if (g_bac_init) return;
    // The 24 letters in Baconian order, i.e. A..Z with J and V removed.
    static const int order24[24] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8,      /* A B C D E F G H I */
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,  /* K L M N O P Q R S T U */
        22, 23, 24, 25                  /* W X Y Z */
    };
    for (int v = 0; v < 24; v++) {
        g_bac_val_to_letter[v] = order24[v];
        g_bac_letter_to_val[order24[v]] = v;
    }
    g_bac_letter_to_val[9]  = g_bac_letter_to_val[8];   // J -> I's code
    g_bac_letter_to_val[21] = g_bac_letter_to_val[20];  // V -> U's code
    g_bac_init = true;
}

int baconian_encrypt(const int pt[], int n, int out_bits[]) {
    baconian_init();
    int o = 0;
    for (int i = 0; i < n; i++) {
        int l = pt[i];
        if (l < 0 || l >= 26) continue;                 // skip non-letters defensively
        int v = g_bac_letter_to_val[l];
        for (int b = 4; b >= 0; b--)                     // most-significant symbol first
            out_bits[o++] = (v >> b) & 1;               // 0 = BAC_A, 1 = BAC_B
    }
    return o;
}

int baconian_decode(const int bits[], int nbits, int out[],
                    int filler, int *n_tokens, int *n_valid) {
    baconian_init();
    int o = 0, ntok = 0, nval = 0;
    int ngroups = nbits / 5;                            // a trailing partial group is dropped
    for (int g = 0; g < ngroups; g++) {
        int v = 0;
        for (int b = 0; b < 5; b++)
            v = v * 2 + (bits[g * 5 + b] ? 1 : 0);
        if (v <= 23) { out[o++] = g_bac_val_to_letter[v]; nval++; }
        else         { out[o++] = filler; }
        ntok++;
    }
    if (n_tokens) *n_tokens = ntok;
    if (n_valid)  *n_valid  = nval;
    return o;
}

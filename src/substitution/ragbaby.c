#include "colossus.h"
#include <ctype.h>

// =====================================================================
//  Ragbaby primitives (TYPE ragbaby)
// =====================================================================
//
// The ACA Ragbaby (see cryptogram.org .../ciphers/Ragbaby.pdf) enciphers over a keyed alphabet.
// Historically the alphabet is 24 letters with I/J and W/X paired (so plaintext J -> I, X -> W);
// a 26- or 36-letter keyed alphabet could be used, so the math here is written for a generic
// `alpha`-letter alphabet (the solver forces 24 via init_alphabet_ragbaby()).
//
// NUMBERING: number each *letter* of the plaintext. The first letter of the first word is 1, the
// first letter of the second word is 2, ... (the per-word start increments by one); within a word
// the number increments by one per letter. The sequence runs 1..alpha and repeats (alpha+1 == 1),
// i.e. the shift used is (word_index + within_word_pos) mod alpha. Word divisions (whitespace) are
// structural; a hyphen/apostrophe keeps a word single (any non-letter passes through and is NOT
// numbered).
//
// ENCIPHER: move the plaintext letter `num` places to the RIGHT in the keyed alphabet (mod alpha):
//   CT = KA[(idx_KA(pt) + num) mod alpha].   DECIPHER: move LEFT: pt = KA[(idx_KA(ct) - num) mod alpha].
//
// Worked example (verified against the ACA PDF): KA = GROSBEAKCDFHILMNPQTUVWYZ (24 letters),
//   pt "word divisions are kept" -> CT "YBBL HNGQDUFGL DEF HFYR".

// Build a keyed alphabet: keyword letters (as alphabet indices, duplicates and out-of-range
// dropped) then the remaining indices ascending. `keyword` is ASCII; folded through g_char_to_idx
// (so a J/X in the keyword pairs to I/W under the Ragbaby alphabet).
void ragbaby_build_keyed_alphabet(int ka[], const char *keyword, int alpha) {
    char used[128];
    for (int i = 0; i < alpha; i++) used[i] = 0;
    int n = 0;
    if (keyword) {
        for (const char *p = keyword; *p; p++) {
            unsigned char uc = (unsigned char) toupper((unsigned char) *p);
            if (uc >= 128) continue;
            int idx = g_char_to_idx[uc];
            if (idx < 0 || idx >= alpha || used[idx]) continue;
            used[idx] = 1;
            ka[n++] = idx;
        }
    }
    for (int i = 0; i < alpha; i++)
        if (!used[i]) ka[n++] = i;
}

// Parse a spaced string into a LETTER stream + parallel per-letter shift number. Whitespace
// separates words (advancing the per-word start number); any other non-letter passes through
// (not numbered, not a separator -- so hyphenated / apostrophe words stay single). Only letters
// (chars whose g_char_to_idx is a valid alphabet index) are emitted into letters[]/num[].
void ragbaby_number_stream(const char *str, int alpha, int letters[], int num[], int *out_len) {
    int word_index = 0;      // becomes 1 at the first letter of the first word
    int within = 0;          // 0-based position within the current word (letters only)
    bool in_word = false;    // seen a letter since the last whitespace run
    int n = 0;
    for (const char *p = str; *p; p++) {
        unsigned char uc = (unsigned char) *p;
        if (isspace(uc)) { in_word = false; continue; }
        int idx = (uc < 128) ? g_char_to_idx[toupper(uc)] : -1;
        if (idx < 0 || idx >= alpha) continue;     // non-alphabet char: pass through, not numbered
        if (!in_word) { word_index++; within = 0; in_word = true; }   // first letter of a new word
        letters[n] = idx;
        num[n] = (word_index + within) % alpha;
        within++;
        n++;
    }
    *out_len = n;
}

void ragbaby_encrypt(const int in[], const int num[], int len, const int ka[], const int ka_inv[],
                     int alpha, int out[]) {
    for (int i = 0; i < len; i++) {
        int p = ka_inv[in[i]];                     // position of the plaintext letter in the KA
        out[i] = ka[((p + num[i]) % alpha + alpha) % alpha];
    }
}

void ragbaby_decrypt(const int in[], const int num[], int len, const int ka[], const int ka_inv[],
                     int alpha, int out[]) {
    for (int i = 0; i < len; i++) {
        int p = ka_inv[in[i]];                     // position of the cipher letter in the KA
        out[i] = ka[((p - num[i]) % alpha + alpha) % alpha];
    }
}

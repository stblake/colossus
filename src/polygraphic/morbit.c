//
//  Morbit cipher primitives -- see morbit.h for the full description.
//
//  The Morse table + reverse-token trick mirror fracmorse.c / pollux.c, kept as a
//  small self-contained copy so those stay byte-identical (the seriated_playfair
//  precedent). The solver needs only morbit_decrypt(); morbit_encrypt serves the
//  generator + unit tests.
//

#include "morbit.h"

// International Morse for A..Z (dot = '.', dash = '-'), longest code is 4 symbols.
static const char *const MB_MORSE[26] = {
    ".-",   "-...", "-.-.", "-..",  ".",    "..-.", "--.",  "....", "..",   ".---",
    "-.-",  ".-..", "--",   "-.",   "---",  ".--.", "--.-", ".-.",  "...",  "-",
    "..-",  "...-", ".--",  "-..-", "-.--", "--.."
};

// Forward table: per-letter DOT/DASH symbols + length. Reverse table: a codeword
// encoded as v = 1, then v = v*2 + (sym==DASH) for each symbol (so v in [2..31] for
// lengths 1..4, all distinct) maps back to a letter; unused v's stay -1 (invalid
// codeword). Built once.
static int  g_mb_sym[26][4];
static int  g_mb_len[26];
static int  g_mb_rev[32];
static bool g_mb_init = false;

void morbit_init(void) {
    if (g_mb_init) return;
    for (int v = 0; v < 32; v++) g_mb_rev[v] = -1;
    for (int l = 0; l < 26; l++) {
        const char *m = MB_MORSE[l];
        int len = 0, v = 1;
        for (int i = 0; m[i]; i++) {
            int s = (m[i] == '-') ? MB_DASH : MB_DOT;
            g_mb_sym[l][len++] = s;
            v = v * 2 + (s == MB_DASH ? 1 : 0);
        }
        g_mb_len[l] = len;
        g_mb_rev[v] = l;
    }
    g_mb_init = true;
}

// Symbol-stream scratch, kept off the stack (each letter is <= 4 code symbols plus
// up to 2 separators, and one trailing pad X).
#define MB_STREAM_MAX (6 * MAX_CIPHER_LENGTH + 8)
static int g_mb_stream[MB_STREAM_MAX];

int morbit_encrypt(const int pt[], int n, const int key[10], int out[]) {
    morbit_init();

    // Invert the key: pair -> digit. Require a valid bijection (each of the 9 pairs
    // owned by exactly one digit 1..9).
    int inv[9];
    for (int p = 0; p < 9; p++) inv[p] = -1;
    for (int d = 1; d <= 9; d++) {
        int p = key[d];
        if (p < 0 || p > 8 || inv[p] != -1) return -1;   // out of range or duplicate pair
        inv[p] = d;
    }
    for (int p = 0; p < 9; p++) if (inv[p] < 0) return -1;

    // Build the Morse symbol stream: letters joined by a single X, word breaks (a
    // negative sentinel) by XX, no leading/trailing X.
    int *s = g_mb_stream, slen = 0;
    int emitted = 0;        // have we emitted a letter yet?
    int seps = 0;           // pending separator X count before the next letter
    for (int i = 0; i < n; i++) {
        int l = pt[i];
        if (l < 0) { if (emitted) seps = 2; continue; }   // word divider -> XX
        if (l >= 26) continue;                            // skip non-letters defensively
        if (emitted) {
            if (seps == 0) seps = 1;                      // default inter-letter separator
            for (int k = 0; k < seps && slen < MB_STREAM_MAX; k++) s[slen++] = MB_X;
        }
        for (int k = 0; k < g_mb_len[l] && slen < MB_STREAM_MAX; k++) s[slen++] = g_mb_sym[l][k];
        emitted = 1;
        seps = 0;
    }

    // Morbit takes the stream in PAIRS; pad one trailing X if the length is odd so it
    // groups into whole pairs (the ACA convention -- verified against the worked
    // example, where "once upon a time" is 45 symbols padded to 46).
    if ((slen & 1) && slen < MB_STREAM_MAX) s[slen++] = MB_X;

    // Each pair (top, bottom) -> the digit the key assigns to pair 3*top+bottom.
    int c = 0;
    for (int j = 0; j + 1 < slen; j += 2)
        out[c++] = inv[3 * s[j] + s[j + 1]];
    return c;
}

int morbit_decrypt(const int cipher[], int clen, const int key[10], int out[],
                   int filler, int *n_tokens, int *n_valid) {
    morbit_init();

    int o = 0, ntok = 0, nval = 0;
    int runlen = 0, v = 1;     // current DOT/DASH codeword
    int xrun = 0;              // current consecutive-X count
    int total = 2 * clen;      // two Morse symbols per cipher digit
    for (int i = 0; i <= total; i++) {
        int sym;
        if (i < total) {
            int pair = key[cipher[i / 2]];               // the digit's pair index 0..8
            sym = (i & 1) ? (pair % 3) : (pair / 3);     // even = top, odd = bottom
        } else {
            sym = MB_X;                                  // sentinel X flushes the final run
        }
        if (sym == MB_X) {
            if (runlen > 0) {                            // flush a completed codeword token
                int letter = (runlen <= 4) ? g_mb_rev[v] : -1;
                if (letter >= 0) { out[o++] = letter; nval++; }
                else             { out[o++] = filler; }
                ntok++;
                runlen = 0; v = 1;
            }
            xrun++;
        } else {
            if (xrun >= 3) ntok++;                       // a closed illegal xxx+ run -> invalid token
            xrun = 0;
            if (runlen < 4) v = v * 2 + (sym == MB_DASH ? 1 : 0);
            runlen++;                                    // runlen may exceed 4 -> invalid codeword
        }
    }
    if (xrun - 1 >= 3) ntok++;                            // trailing run (drop the sentinel X)

    if (n_tokens) *n_tokens = ntok;
    if (n_valid)  *n_valid  = nval;
    return o;
}

//
//  Pollux cipher primitives -- see pollux.h for the full description.
//
//  The Morse table + reverse-token trick mirror fracmorse.c, kept as a small
//  self-contained copy so fracmorse.c stays byte-identical (the seriated_playfair
//  precedent). The solver needs only pollux_decrypt(); pollux_encrypt serves the
//  generator + unit tests.
//

#include "pollux.h"

// International Morse for A..Z (dot = '.', dash = '-'), longest code is 4 symbols.
static const char *const PX_MORSE[26] = {
    ".-",   "-...", "-.-.", "-..",  ".",    "..-.", "--.",  "....", "..",   ".---",
    "-.-",  ".-..", "--",   "-.",   "---",  ".--.", "--.-", ".-.",  "...",  "-",
    "..-",  "...-", ".--",  "-..-", "-.--", "--.."
};

// Forward table: per-letter DOT/DASH symbols + length. Reverse table: a codeword
// encoded as v = 1, then v = v*2 + (sym==DASH) for each symbol (so v in [2..31] for
// lengths 1..4, all distinct) maps back to a letter; unused v's stay -1 (invalid
// codeword). Built once.
static int  g_px_sym[26][4];
static int  g_px_len[26];
static int  g_px_rev[32];
static bool g_px_init = false;

void pollux_init(void) {
    if (g_px_init) return;
    for (int v = 0; v < 32; v++) g_px_rev[v] = -1;
    for (int l = 0; l < 26; l++) {
        const char *m = PX_MORSE[l];
        int len = 0, v = 1;
        for (int i = 0; m[i]; i++) {
            int s = (m[i] == '-') ? PX_DASH : PX_DOT;
            g_px_sym[l][len++] = s;
            v = v * 2 + (s == PX_DASH ? 1 : 0);
        }
        g_px_len[l] = len;
        g_px_rev[v] = l;
    }
    g_px_init = true;
}

// Symbol-stream scratch, kept off the stack (each letter is <= 4 code symbols plus
// up to 2 separators).
#define PX_STREAM_MAX (6 * MAX_CIPHER_LENGTH + 8)
static int g_px_stream[PX_STREAM_MAX];

int pollux_encrypt(const int pt[], int n, const int key[10], int out[]) {
    pollux_init();

    // Group the digits by the symbol each is assigned to.
    int dig[3][10], ndig[3] = {0, 0, 0};
    for (int d = 0; d < 10; d++) {
        int e = key[d];
        if (e < 0 || e > 2) return -1;
        dig[e][ndig[e]++] = d;
    }
    if (ndig[PX_DOT] == 0 || ndig[PX_DASH] == 0 || ndig[PX_X] == 0) return -1;

    // Build the Morse symbol stream: letters joined by a single X, word breaks (a
    // negative sentinel) by XX, no leading/trailing X.
    int *s = g_px_stream, slen = 0;
    int emitted = 0;        // have we emitted a letter yet?
    int seps = 0;           // pending separator X count before the next letter
    for (int i = 0; i < n; i++) {
        int l = pt[i];
        if (l < 0) { if (emitted) seps = 2; continue; }   // word divider -> XX
        if (l >= 26) continue;                            // skip non-letters defensively
        if (emitted) {
            if (seps == 0) seps = 1;                      // default inter-letter separator
            for (int k = 0; k < seps && slen < PX_STREAM_MAX; k++) s[slen++] = PX_X;
        }
        for (int k = 0; k < g_px_len[l] && slen < PX_STREAM_MAX; k++) s[slen++] = g_px_sym[l][k];
        emitted = 1;
        seps = 0;
    }

    // Map each Morse symbol to a random digit assigned to it (polyphonic).
    for (int j = 0; j < slen; j++) {
        int e = s[j];
        out[j] = dig[e][rand_int(0, ndig[e])];
    }
    return slen;
}

int pollux_decrypt(const int cipher[], int clen, const int key[10], int out[],
                   int filler, int *n_tokens, int *n_valid) {
    pollux_init();

    int o = 0, ntok = 0, nval = 0;
    int runlen = 0, v = 1;     // current DOT/DASH codeword
    int xrun = 0;              // current consecutive-X count
    for (int i = 0; i <= clen; i++) {
        int sym = (i < clen) ? key[cipher[i]] : PX_X;    // sentinel X flushes the final run
        if (sym == PX_X) {
            if (runlen > 0) {                            // flush a completed codeword token
                int letter = (runlen <= 4) ? g_px_rev[v] : -1;
                if (letter >= 0) { out[o++] = letter; nval++; }
                else             { out[o++] = filler; }
                ntok++;
                runlen = 0; v = 1;
            }
            xrun++;
        } else {
            if (xrun >= 3) ntok++;                       // a closed illegal xxx+ run -> invalid token
            xrun = 0;
            if (runlen < 4) v = v * 2 + (sym == PX_DASH ? 1 : 0);
            runlen++;                                    // runlen may exceed 4 -> invalid codeword
        }
    }
    if (xrun - 1 >= 3) ntok++;                            // trailing run (drop the sentinel X)

    if (n_tokens) *n_tokens = ntok;
    if (n_valid)  *n_valid  = nval;
    return o;
}

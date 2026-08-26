//
//  Compressocrat cipher primitives (ACA type 92) -- see compressocrat.h for the full description.
//
//  A FIXED prefix-free "irregular" alphabet maps each letter to a variable-length code over the
//  digits {1,2,3}; the plaintext's codes are concatenated, padded with the digit 1 to a multiple
//  of 3, grouped into trigraphs, and each trigraph (RANK r = 9(a-1)+3(b-1)+(c-1), 333 impossible
//  so r in 0..25) mapped to a ciphertext letter through the keyed alphabet sigma[r]. The unknown
//  is sigma; the compression table is fixed. This is the fractionation twin of Fractionated Morse
//  (fracmorse.c): the Morse table is replaced by this Huffman-style {1,2,3} code, and the split-
//  on-x parse by a GREEDY prefix parse (prefix-free => the first complete match is unique).
//

#include "colossus.h"
#include "compressocrat.h"

// The fixed Compressocrat compression alphabet: A..Z -> a prefix-free code over the digits 1,2,3
// (common letters short, rare letters long). Longest code is 6 digits (J, X, Z).
static const char *const CR_CODE[26] = {
    "13",     "32112",  "1112",   "213",    "31",     "3213",   "32113",  "113",
    "322",    "321112", "11112",  "212",    "2111",   "23",     "22",     "3212",
    "11113",  "323",    "112",    "12",     "1113",   "11111",  "2112",   "321111",
    "2113",   "321113"
};

#define CR_MAX_CODE_LEN 6

// Forward table: per-letter digit array (values 1..3) + length. Reverse table: a code hashed as
// v = 1, then v = v*4 + digit for each digit (digits are 1..3, so v is distinct across all codes
// and prefixes; v <= 8191 for <=6 digits) maps back to a letter; unused v's are -1 (not a code).
// Built once, thread-locally.
static _Thread_local int  g_cr_sym[26][CR_MAX_CODE_LEN];
static _Thread_local int  g_cr_len[26];
static _Thread_local int  g_cr_rev[8192];
static _Thread_local bool g_cr_init = false;

static void compressocrat_init(void) {
    if (g_cr_init) return;
    for (int v = 0; v < 8192; v++) g_cr_rev[v] = -1;
    for (int l = 0; l < 26; l++) {
        const char *c = CR_CODE[l];
        int len = 0, v = 1;
        for (int i = 0; c[i]; i++) {
            int d = c[i] - '0';                          // digit 1..3
            g_cr_sym[l][len++] = d;
            v = v * 4 + d;
        }
        g_cr_len[l] = len;
        g_cr_rev[v] = l;
    }
    g_cr_init = true;
}

// Symbol-stream scratch, kept off the stack. Encrypt needs up to 6 digits per plaintext letter
// plus <=2 padding; decrypt needs exactly 3 per ciphertext letter.
#define CR_STREAM_MAX (6 * MAX_CIPHER_LENGTH + 8)
static _Thread_local int g_cr_stream[CR_STREAM_MAX];

// Encipher plaintext (n letters, indices 0..25) under the keyed alphabet sigma (rank -> letter).
// Writes the ciphertext letters into out[] and returns the ciphertext length C. out[] must hold
// up to C = (sum of code lengths, rounded up to a multiple of 3) / 3 letters.
int compressocrat_encrypt(const int plain[], int n, const int sigma[], int out[]) {
    compressocrat_init();
    int *s = g_cr_stream;
    int slen = 0;
    for (int i = 0; i < n; i++) {
        int l = plain[i];
        if (l < 0 || l >= 26) continue;                 // skip any non-letter defensively
        for (int k = 0; k < g_cr_len[l] && slen < CR_STREAM_MAX; k++) s[slen++] = g_cr_sym[l][k];
    }
    while (slen % 3 != 0 && slen < CR_STREAM_MAX) s[slen++] = 1;   // pad with the digit 1
    int clen = slen / 3;
    for (int g = 0; g < clen; g++) {
        int r = 9 * (s[3 * g] - 1) + 3 * (s[3 * g + 1] - 1) + (s[3 * g + 2] - 1);  // 0..25 (333 impossible)
        out[g] = sigma[r];
    }
    return clen;
}

// Decipher ciphertext (clen letters) under the keyed alphabet sigma. Expands each letter to its
// trigraph's 3 digits, then greedily prefix-parses the {1,2,3} stream: a legal codeword -> its
// letter (a valid token), a 6-digit run with no match -> `filler` (an invalid token). Trailing
// padding (a run of <=2 ones completing no code) is dropped. Writes one letter per token into
// out[] and returns the token count; *n_tokens / *n_valid receive the total / legal-codeword
// token counts (either pointer may be NULL). out[] must hold up to ~3*clen/2 letters (worst case:
// every token the 2-digit minimum). On the true key this returns exactly the original plaintext.
int compressocrat_decrypt(const int cipher[], int clen, const int sigma[], int out[],
                          int filler, int *n_tokens, int *n_valid) {
    compressocrat_init();
    int inv[26];
    for (int r = 0; r < 26; r++) inv[sigma[r]] = r;      // letter -> trigraph rank

    int *s = g_cr_stream;
    for (int i = 0; i < clen; i++) {
        int r = inv[cipher[i]];
        s[3 * i]     = r / 9 + 1;                         // digits 1..3 (rank never encodes 333)
        s[3 * i + 1] = (r / 3) % 3 + 1;
        s[3 * i + 2] = r % 3 + 1;
    }
    int slen = 3 * clen;

    int o = 0, ntok = 0, nval = 0;
    int runlen = 0, v = 1, all_ones = 1;
    for (int i = 0; i < slen; i++) {
        int d = s[i];
        v = v * 4 + d; runlen++;
        if (d != 1) all_ones = 0;
        if (g_cr_rev[v] >= 0) {                          // a complete, legal codeword
            out[o++] = g_cr_rev[v]; nval++; ntok++;
            runlen = 0; v = 1; all_ones = 1;
        } else if (runlen == CR_MAX_CODE_LEN) {          // 6 digits with no match -> invalid token
            out[o++] = filler; ntok++;
            runlen = 0; v = 1; all_ones = 1;
        }
    }
    if (runlen > 0) {                                    // a trailing partial run
        if (!(all_ones && runlen <= 2)) {                // padding is <=2 ones; anything else is invalid
            out[o++] = filler; ntok++;
        }
    }
    if (n_tokens) *n_tokens = ntok;
    if (n_valid)  *n_valid  = nval;
    return o;
}

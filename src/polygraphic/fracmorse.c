//
//  Fractionated Morse cipher primitives (Morse fractionation over a keyed 26-letter alphabet).
//
//  The plaintext is first written in Morse code with a single separator between letters:
//
//     S = morse(p0) x morse(p1) x ... x morse(p_{N-1})
//
//  where the three symbols are DOT (.), DASH (-) and X (x, the separator). There is exactly
//  ONE x between consecutive letters and none at the ends (Colossus works on bare A..Z with
//  no word divisions, so the ACA "xx between words" convention does not apply -- see the
//  header comment in fracmorse_solver.c). S is then padded with trailing x until its length
//  is a multiple of 3 (0..2 x's; the body has only single x separators, so the run xxx can
//  never occur), grouped into trigraphs, and each trigraph mapped to one ciphertext letter.
//
//  The 26 trigraphs are all length-3 combinations of {DOT,DASH,X} in the order DOT<DASH<X
//  EXCEPT xxx: a trigraph (a,b,c) has RANK r = 9a + 3b + c (base 3, digit 0=DOT,1=DASH,2=X),
//  and because xxx = rank 26 never occurs every group has r in 0..25. The keyed alphabet
//  sigma (a permutation of A..Z, built keyword-then-ascending-tail like every ACA mixed
//  alphabet) is written across the trigraphs in rank order, so ciphertext letter = sigma[r].
//  The ciphertext is C = len(S_padded)/3 letters -- a length CHANGE from the N plaintext
//  letters.
//
//  Decryption inverts sigma (letter -> rank), expands each ciphertext letter to its 3 Morse
//  symbols (a fixed, key-independent 3C-symbol stream), splits the stream on x, and decodes
//  each maximal {DOT,DASH} run as a Morse codeword. On the true key this recovers exactly the
//  N plaintext letters (the trailing pad x's leave only empty tokens, dropped); a WRONG key
//  can leave runs that are not legal codewords (invalid) -- fracmorse_decrypt reports the
//  token / valid-token counts (the solver folds validity into its score) and emits a caller-
//  supplied filler letter for each invalid token.
//
//  The solver needs only fracmorse_decrypt(); fracmorse_encrypt serves the generator + tests.
//

#include "colossus.h"

// Morse symbols (also the base-3 digits of a trigraph rank).
enum { FM_DOT = 0, FM_DASH = 1, FM_X = 2 };

// International Morse for A..Z (dot = '.', dash = '-'), longest code is 4 symbols.
static const char *const FM_MORSE[26] = {
    ".-",   "-...", "-.-.", "-..",  ".",    "..-.", "--.",  "....", "..",   ".---",
    "-.-",  ".-..", "--",   "-.",   "---",  ".--.", "--.-", ".-.",  "...",  "-",
    "..-",  "...-", ".--",  "-..-", "-.--", "--.."
};

// Forward table: per-letter symbol array (DOT/DASH) + length. Reverse table: a codeword
// encoded as v = 1, then v = v*2 + (sym==DASH) for each symbol (so v in [2..31] for lengths
// 1..4, all distinct) maps back to a letter; unused v's are -1 (invalid codeword). Built once.
static int  g_fm_sym[26][4];
static int  g_fm_len[26];
static int  g_fm_rev[32];
static bool g_fm_init = false;

static void fracmorse_init(void) {
    if (g_fm_init) return;
    for (int v = 0; v < 32; v++) g_fm_rev[v] = -1;
    for (int l = 0; l < 26; l++) {
        const char *m = FM_MORSE[l];
        int len = 0, v = 1;
        for (int i = 0; m[i]; i++) {
            int s = (m[i] == '-') ? FM_DASH : FM_DOT;
            g_fm_sym[l][len++] = s;
            v = v * 2 + (s == FM_DASH ? 1 : 0);
        }
        g_fm_len[l] = len;
        g_fm_rev[v] = l;
    }
    g_fm_init = true;
}

// Symbol-stream scratch, kept off the stack (encrypt needs up to ~5 symbols per plaintext
// letter: <=4 code + 1 separator; decrypt needs exactly 3 per ciphertext letter).
#define FM_STREAM_MAX (5 * MAX_CIPHER_LENGTH + 8)
static int g_fm_stream[FM_STREAM_MAX];

// Encipher plaintext (n letters, indices 0..25) under the keyed alphabet sigma (rank -> letter).
// Writes the ciphertext letters into out[] and returns the ciphertext length C. out[] must hold
// up to C = (sum of Morse lengths + (n-1) separators, rounded up to a multiple of 3) / 3 letters.
int fracmorse_encrypt(const int plain[], int n, const int sigma[], int out[]) {
    fracmorse_init();
    int *s = g_fm_stream;
    int slen = 0;
    for (int i = 0; i < n; i++) {
        int l = plain[i];
        if (l < 0 || l >= 26) continue;                 // skip any non-letter defensively
        for (int k = 0; k < g_fm_len[l] && slen < FM_STREAM_MAX; k++) s[slen++] = g_fm_sym[l][k];
        if (i + 1 < n && slen < FM_STREAM_MAX) s[slen++] = FM_X;   // single separator, none at end
    }
    while (slen % 3 != 0 && slen < FM_STREAM_MAX) s[slen++] = FM_X; // pad to a whole trigraph
    int clen = slen / 3;
    for (int g = 0; g < clen; g++) {
        int r = 9 * s[3 * g] + 3 * s[3 * g + 1] + s[3 * g + 2];    // rank in 0..25 (xxx impossible)
        out[g] = sigma[r];
    }
    return clen;
}

// Decipher ciphertext (clen letters) under the keyed alphabet sigma. Expands each letter to
// its trigraph's 3 Morse symbols, splits on X, and decodes each maximal DOT/DASH run: a legal
// codeword -> its letter (a valid token), anything else -> `filler` (an invalid token). Empty
// runs (consecutive X) are not tokens. Writes one letter per token into out[] and returns the
// token count; *n_tokens / *n_valid receive the total / legal-codeword token counts (either
// pointer may be NULL). out[] must hold up to ~3*clen/2 letters (worst case: every token one
// symbol). On the true key this returns exactly the original plaintext (all tokens valid).
int fracmorse_decrypt(const int cipher[], int clen, const int sigma[], int out[],
                      int filler, int *n_tokens, int *n_valid) {
    fracmorse_init();
    int inv[26];
    for (int r = 0; r < 26; r++) inv[sigma[r]] = r;      // letter -> trigraph rank

    int *s = g_fm_stream;
    for (int i = 0; i < clen; i++) {
        int r = inv[cipher[i]];
        s[3 * i]     = r / 9;
        s[3 * i + 1] = (r / 3) % 3;
        s[3 * i + 2] = r % 3;
    }
    int slen = 3 * clen;

    int o = 0, ntok = 0, nval = 0;
    int runlen = 0, v = 1;
    for (int i = 0; i <= slen; i++) {
        int sym = (i < slen) ? s[i] : FM_X;              // sentinel X flushes the final run
        if (sym == FM_X) {
            if (runlen > 0) {                            // a completed token
                int letter = (runlen <= 4) ? g_fm_rev[v] : -1;
                if (letter >= 0) { out[o++] = letter; nval++; }
                else             { out[o++] = filler; }
                ntok++;
            }
            runlen = 0; v = 1;
        } else {
            if (runlen < 4) v = v * 2 + (sym == FM_DASH ? 1 : 0);
            runlen++;                                    // runlen may exceed 4 -> invalid
        }
    }
    if (n_tokens) *n_tokens = ntok;
    if (n_valid)  *n_valid  = nval;
    return o;
}

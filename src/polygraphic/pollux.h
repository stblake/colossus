#ifndef POLLUX_H
#define POLLUX_H
#include "colossus.h"

// =====================================================================
//  Pollux cipher primitives (Morse over a digit -> {dot,dash,x} map)
// =====================================================================
//
//  The plaintext is written in International Morse with a single 'x' separator
//  between letters and a double 'xx' between words (a run of three or more x's is
//  impossible in valid Morse). The three stream symbols are DOT, DASH and X. The
//  KEY assigns each of the ten digits 0..9 to one of those three symbols (usually
//  4 dots / 3 dashes / 3 x, but any split with >= 1 of each is legal), and the
//  ciphertext is the resulting DIGIT stream -- one digit per Morse symbol.
//
//  Encryption is POLYPHONIC: for each Morse symbol the encoder picks an arbitrary
//  (here random) digit among those the key assigns to that symbol. Decryption is
//  fully deterministic given the key: map every digit to its symbol, split the
//  symbol stream on x, and decode each maximal DOT/DASH run as a Morse codeword.

// Morse symbols == the three values a key digit maps to.
enum { PX_DOT = 0, PX_DASH = 1, PX_X = 2 };

// Build the A..Z Morse tables (idempotent). Called by encrypt/decrypt; exposed for
// tests that inspect the tables directly.
void pollux_init(void);

// Encipher plaintext pt[0..n-1] under key[0..9] (each entry PX_DOT/PX_DASH/PX_X).
// pt[i] in 0..25 is a letter; pt[i] < 0 is a WORD-DIVIDER sentinel (spaces carried
// through from the input). Letters are joined by a single X, word breaks by XX, with
// no leading/trailing X; each Morse symbol is emitted as a RANDOM digit the key
// assigns to it (uses the global RNG -- seed it for a reproducible cipher). Writes
// the digit stream (values 0..9) into out[] and returns its length C. Requires each
// of PX_DOT/PX_DASH/PX_X to own >= 1 digit; returns -1 otherwise.
int pollux_encrypt(const int pt[], int n, const int key[10], int out[]);

// Decipher the digit stream cipher[0..clen-1] (values 0..9) under key[0..9]. Maps
// each digit to its Morse symbol, splits the symbol stream on X, and decodes each
// maximal DOT/DASH run as a Morse codeword: a legal codeword -> its letter (a valid
// token), anything else (a codeword > 4 symbols, or an unknown pattern) -> `filler`
// (an invalid token). A run of >= 3 consecutive X is illegal Morse and is counted as
// one extra invalid token. Word boundaries are NOT emitted -- out[] receives LETTERS
// ONLY (one per codeword token) so a char-for-char plaintext compare works. Returns
// the number of letters written; *n_tokens / *n_valid receive the total-token /
// legal-codeword-token counts (either pointer may be NULL). out[] must hold up to
// clen letters. On the true key this recovers exactly the plaintext (all valid).
int pollux_decrypt(const int cipher[], int clen, const int key[10], int out[],
                   int filler, int *n_tokens, int *n_valid);
#endif

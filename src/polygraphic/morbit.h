#ifndef MORBIT_H
#define MORBIT_H
#include "colossus.h"

// =====================================================================
//  Morbit cipher primitives (Morse taken in PAIRS over a pair -> digit map)
// =====================================================================
//
//  The plaintext is written in International Morse with a single 'x' separator
//  between letters and a double 'xx' between words -- exactly the Fractionated
//  Morse / Pollux convention (a run of three or more x's is impossible in valid
//  Morse). The three stream symbols are DOT, DASH and X. Morbit is the digraphic
//  sibling of Pollux: instead of one digit per symbol, the symbol stream is taken
//  in PAIRS (units of 2, padding one trailing X if the count is odd), and each of
//  the 9 possible pairs is enciphered as a single digit 1..9.
//
//  The 9 pairs are ordered base-3 (DOT=0, DASH=1, X=2): pair = 3*top + bottom, i.e.
//  (.,.) (.,-) (.,x) (-,.) (-,-) (-,x) (x,.) (x,-) (x,x). A 9-letter keyword assigns
//  the digits 1..9 to these pairs by the keyword letters' alphabetical rank -- so the
//  KEY is a bijection pair <-> digit and the full keyspace is 9! = 362,880.
//
//  Encryption is DETERMINISTIC (a bijection, unlike Pollux's polyphonic symbol map).
//  Decryption is also deterministic: map every digit back to its two Morse symbols,
//  split the symbol stream on X, and decode each maximal DOT/DASH run as a codeword.

// Morse symbols == the three values a pair's top/bottom cell can hold.
enum { MB_DOT = 0, MB_DASH = 1, MB_X = 2 };

// Build the A..Z Morse tables (idempotent). Called by encrypt/decrypt; exposed for
// tests that inspect the tables directly.
void morbit_init(void);

// Encipher plaintext pt[0..n-1] under key[1..9] (each entry a PAIR index 0..8;
// key[0] is unused -- a Morbit ciphertext never contains a 0). pt[i] in 0..25 is a
// letter; pt[i] < 0 is a WORD-DIVIDER sentinel (spaces carried through from the
// input). Letters are joined by a single X, word breaks by XX, with no leading/
// trailing X; the symbol stream is padded with one trailing X if odd, then taken in
// pairs (top, bottom) -> the digit the key assigns to pair 3*top+bottom. Writes the
// digit stream (values 1..9) into out[] and returns its length C = ceil(stream/2).
// Requires key[1..9] to be a bijection onto the 9 pairs; returns -1 otherwise.
int morbit_encrypt(const int pt[], int n, const int key[10], int out[]);

// Decipher the digit stream cipher[0..clen-1] (values 1..9) under key[1..9]. Maps
// each digit to its pair's two Morse symbols, splits the symbol stream on X, and
// decodes each maximal DOT/DASH run as a Morse codeword: a legal codeword -> its
// letter (a valid token), anything else (a codeword > 4 symbols, or an unknown
// pattern) -> `filler` (an invalid token). A run of >= 3 consecutive X is illegal
// Morse and is counted as one extra invalid token. Word boundaries are NOT emitted --
// out[] receives LETTERS ONLY (one per codeword token) so a char-for-char plaintext
// compare works. Returns the number of letters written; *n_tokens / *n_valid receive
// the total-token / legal-codeword-token counts (either pointer may be NULL). out[]
// must hold up to clen letters. On the true key this recovers exactly the plaintext.
int morbit_decrypt(const int cipher[], int clen, const int key[10], int out[],
                   int filler, int *n_tokens, int *n_valid);
#endif

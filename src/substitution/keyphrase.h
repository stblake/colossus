#ifndef KEYPHRASE_H
#define KEYPHRASE_H
#include "colossus.h"

// =====================================================================
//  Key Phrase cipher primitives (TYPE keyphrase, ACA)
// =====================================================================
//
// The cipher alphabet is a 26-letter PHRASE which must be "complete" (exactly 26 letters,
// matched letter-for-letter to a straight plaintext alphabet a..z). Enciphering sends the
// plaintext letter p (0..25) to the phrase letter at position p:  CT[i] = keyphrase[p].
// Because a natural phrase REPEATS letters, several plaintext letters map to the SAME
// ciphertext letter, so DECODE is many-to-one / AMBIGUOUS -- a ciphertext letter stands for
// a SET of plaintext letters, resolved only by context (the solver's inner beam-Viterbi).
// Word divisions and punctuation are retained in the ciphertext.
//
//   Worked example (ACA), key phrase "GIVEMELIBERTYORGIVEMEDEATH" over straight a..z:
//   "a ciphertext letter may stand for more than one plaintext letter"
//     -> "G VBGIMVMMAM TMMMMV YGT EMGOE ERV YRVM MIGO ROM GTGBOMMAM TMMMMV."
//   (e.g. TMMMMV = "letter" = l,e,t,t,e,r: the four M's decode to e,t,t,e -- one ciphertext
//   letter, four plaintext letters.)
//
// The primitive is deterministic ENCODE only (ambiguous decode is the solver's job). The key
// is the phrase as alphabet indices key[0..25] (key[p] = the ciphertext letter for plaintext p).

#define KEYPHRASE_LEN 26   // the phrase spans the 26 plaintext letters a..z

// Encipher a letter-index stream in[0..len-1] (values 0..25) under the 26-letter phrase key[].
// Non-letter sentinels (values < 0) pass through unchanged.
void keyphrase_encrypt(const int in[], int len, const int key[KEYPHRASE_LEN], int out[]);

// Build a phrase key[] from a phrase STRING: the first 26 A..Z letters of `phrase` become
// key[0..25] (key[p] = index of the p-th phrase letter). Returns the number of letters used
// (must be 26 for a valid "complete" phrase); non-letters in the string are skipped.
int keyphrase_build_key(const char *phrase, int key[KEYPHRASE_LEN]);
#endif

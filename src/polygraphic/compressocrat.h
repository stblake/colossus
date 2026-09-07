//
//  Compressocrat cipher primitives (ACA type 92).
//
//  A Compressocrat is the fractionation twin of Fractionated Morse (fracmorse.c): a FIXED
//  prefix-free "irregular" alphabet maps each of the 26 letters to a variable-length code over
//  the digits {1,2,3} (a Huffman-style compression code -- common letters are short: E=31,
//  T=12, A=13; rare letters long: J=321112, Z=321113). The plaintext's codes are concatenated
//  into one {1,2,3} stream, padded with the digit 1 (once or twice) to a multiple of 3, grouped
//  into trigraphs, and each trigraph mapped to one ciphertext LETTER through a keyed 26-letter
//  alphabet (rank -> letter). The 26 trigraphs are all length-3 combinations of {1,2,3} EXCEPT
//  333 (which never occurs -- no code contains "33", and after a leading 3 the next digit is
//  always 1 or 2), in natural order, exactly as fracmorse excludes xxx. So the ciphertext is a
//  length CHANGE from the N plaintext letters (typically SHORTER -- the "compression").
//
//  The keyed alphabet sigma (a permutation of A..Z, built keyword-then-ascending-tail like every
//  ACA mixed alphabet) is the ONLY unknown; the compression table is fixed. Decoding inverts
//  sigma (letter -> rank), expands each ciphertext letter to its 3 digits (a fixed, key-
//  independent 3C-digit stream), then GREEDILY prefix-parses the stream against the compression
//  code (prefix-free => the first complete match is unique): a legal codeword -> its letter (a
//  valid token), a 6-digit run with no match -> a caller-supplied filler (an invalid token).
//  Trailing padding (<=2 ones that complete no code) is dropped. On the true key this recovers
//  exactly the N plaintext letters (all tokens valid); a WRONG key strands unparseable runs, so
//  compressocrat_decrypt reports the token / valid-token counts and the solver folds the valid
//  fraction into its score (the fracmorse Morse-validity analogue).
//
//  The solver needs only compressocrat_decrypt(); compressocrat_encrypt serves the generator +
//  unit tests. (The keyed-alphabet search moves live in compressocrat_solver.h, like fracmorse.)
//

#ifndef COMPRESSOCRAT_H
#define COMPRESSOCRAT_H

// Encipher plaintext (n letters, indices 0..25) under the keyed alphabet sigma (rank -> letter).
// Writes ciphertext letters into out[] and returns the ciphertext length C.
int compressocrat_encrypt(const int plain[], int n, const int sigma[], int out[]);

// Decipher ciphertext (clen letters) under sigma. Writes one plaintext letter per token into
// out[] (filler for an invalid token) and returns the token count; *n_tokens / *n_valid receive
// the total / legal-codeword token counts (either pointer may be NULL).
int compressocrat_decrypt(const int cipher[], int clen, const int sigma[], int out[],
                          int filler, int *n_tokens, int *n_valid);

#endif

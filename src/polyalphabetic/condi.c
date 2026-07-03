//
// Condi Cipher
//

#include "colossus.h"

/*
   Condi Cipher Logic
   ==================

   The ACA Condi ("CONsecutive DIgraphs" feedback) cipher enciphers over a simple KEYED alphabet
   sigma (a keyword with duplicates dropped, then the remaining letters in order -- optionally
   cyclically shifted) using a PLAINTEXT-FEEDBACK running key: the shift applied to each plaintext
   letter is the position of the PRECEDING plaintext letter in the keyed alphabet; the first letter
   uses a supplied starter offset. Per the ACA:

     "With a starter value or off-set of #, substitute the first plaintext letter by the letter
      found # places further along the alphabet. Then the position of that first plaintext letter
      is the new value for #, the off-set for the next plaintext letter. And so on."

   Let sigma[k] be the letter at position k (0..25) of the keyed alphabet and sigma_inv[x] the
   position of letter x. With u_i = sigma_inv[pt_i] (0-indexed position of plaintext letter i):

     off_1 = starter
     ct_i  = sigma[(u_i + off_i) mod 26]
     off_{i+1} = (u_i + 1) mod 26          // the 1-indexed position of plaintext letter i

   i.e. idx(ct_i) = (idx(pt_i) + idx(pt_{i-1}) + 1) mod 26 for i >= 2. Decryption is CAUSAL: each
   offset needs the already-recovered previous plaintext letter, so it walks left-to-right.

   Only alphabetic letters participate (word divisions / punctuation are skipped by the encipherer,
   and Colossus operates on the A..Z-only stream, so the feedback is over consecutive LETTERS).

   Hand-verified against the ACA worked example: keyword STRANGE, keyed alphabet (shifted)
   "VWXYZSTRANGEBCDFHIJKLMOPQU", starter 25, plaintext "OURSISAVERYGREEN...PENCILANDPAPER"
   enciphers to "MORCPPDNBKE...YJQBTTNK" (see tests/test_condi.c). Note off = 26 mod 26 = 0 gives a
   self-encipherment (a plaintext letter at the last keyed-alphabet position maps to itself).

   The primitives operate directly on the integer-index text arrays; the caller supplies the keyed
   alphabet sigma AND its inverse sigma_inv (the solver caches the inverse per decrypt). The whole
   key the solver recovers is sigma (a 26-permutation, in the key lane) plus the starter (0..25).
*/

// off_1 = starter; ct_i = sigma[(u_i + off_i) mod 26]; off_{i+1} = (u_i + 1) mod 26.
void condi_encrypt(const int plaintext_indices[], int plaintext_len,
    const int sigma[], const int sigma_inv[], int starter, int encrypted[]) {

    int off = ((starter % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    for (int i = 0; i < plaintext_len; i++) {
        int u = sigma_inv[plaintext_indices[i]];               // 0-indexed position of the pt letter
        encrypted[i] = sigma[(u + off) % ALPHABET_SIZE];
        off = (u + 1) % ALPHABET_SIZE;                         // 1-indexed position of this pt letter
    }
}

// Causal inverse: pt_i = sigma[(idx(ct_i) - off_i) mod 26]; off_{i+1} = (idx(pt_i) + 1) mod 26.
void condi_decrypt(const int cipher_indices[], int cipher_len,
    const int sigma[], const int sigma_inv[], int starter, int decrypted[]) {

    int off = ((starter % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    for (int i = 0; i < cipher_len; i++) {
        int v = sigma_inv[cipher_indices[i]];                  // 0-indexed position of the ct letter
        int u = ((v - off) % ALPHABET_SIZE + ALPHABET_SIZE) % ALPHABET_SIZE;
        decrypted[i] = sigma[u];
        off = (u + 1) % ALPHABET_SIZE;                         // feedback off the recovered pt letter
    }
}

//
// Interrupted Key Cipher
//

#include "colossus.h"

/*
   Interrupted Key Cipher Logic
   ============================

   The ACA Interrupted Key cipher is a periodic polyalphabetic base cipher (Vigenere / Variant /
   Beaufort) under a P-letter keyword, but the keying is INTERRUPTED: the key index k returns to
   the first key letter (k = 0) at "break" points and otherwise increments (k = (k+1) mod P). Per
   the ACA: "The plaintext is enciphered with 1, 2, 3 or more letters of the keyword which is
   interrupted at random, by plaintext word division, or according to some other scheme. Return to
   the first key letter each time the keyword is interrupted. The entire keyword must be used at
   least once."

   The base shift math is identical to Progressive Key's, so these primitives reuse
   progkey_base_encrypt / progkey_base_decrypt (colossus.h) with base in {VIG, VAR, BEAU}:

     * Vigenere  E(p,k) = (p + k),  D(c,k) = (c - k)
     * Variant   E(p,k) = (p - k),  D(c,k) = (c + k)
     * Beaufort  E(p,k) = (k - p),  D(c,k) = (k - c)   (self-reciprocal)

   How break points arise defines the interruption STRATEGY (all attacked by intkey_solver.c):

     - CIPHERTEXT interruptor: reset AFTER a chosen ciphertext letter. The break positions are
       determined by the ciphertext alone (keyword-independent) -- so the key index at each
       position is known and the keyword recovers like a plain Vigenere. Encryption is causal
       (the reset depends on the ciphertext being produced): intkey_encrypt_ctint. Decryption is
       equivalent to intkey_decrypt_mask over the mask from intkey_build_mask_ct.

     - PLAINTEXT interruptor: reset AFTER a chosen plaintext letter. The break positions depend on
       the plaintext (hence the keyword), so decryption walks causally: intkey_*_ptint.

     - SUPPLIED BREAKS (random / word-division): the encipherer chose the breaks by some scheme;
       given them as an is_break[] mask the cipher is fully determined: intkey_*_mask. The solver
       either takes them from the user (-breaks) or, blind, searches the mask jointly with the
       keyword.

   Convention: "reset AFTER the trigger" -- the interruptor letter is enciphered normally and the
   NEXT letter restarts the key. The mask form uses is_break[i]==1 to mean "position i starts a new
   group (force k=0 here)"; is_break[0] is implicit (k starts at 0).

   Hand-verified against the ACA worked example (Vigenere, keyword ORANGE): the plaintext
   "thiscipher can be used with any of the periodics" broken into groups of lengths
   4,6,2,3,4,3,1,1,5,1,2,3,5 enciphers to "HYIFQZPUKVQRBSEIJEQKZTVOBPOSZVSGSIICUIPY"
   (see tests/test_intkey.c).

   The primitives operate directly on the integer-index text arrays; the keyword shifts (0..25)
   live in the cycleword lane for the solver, so the keyword IS the recovered per-column shift
   sequence.
*/

// --- general break-mask form (breaks known: ciphertext-interruptor / supplied / joint) ---------

void intkey_encrypt_mask(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int base, const int is_break[]) {

    int k = 0;
    for (int i = 0; i < plaintext_len; i++) {
        if (is_break[i]) k = 0;                                     // force key letter 1 at a break
        encrypted[i] = progkey_base_encrypt(plaintext_indices[i], keyword[k], base);
        if (++k >= P) k = 0;                                        // increment, wrap mod P
    }
}

void intkey_decrypt_mask(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int base, const int is_break[]) {

    int k = 0;
    for (int i = 0; i < cipher_len; i++) {
        if (is_break[i]) k = 0;
        decrypted[i] = progkey_base_decrypt(cipher_indices[i], keyword[k], base);
        if (++k >= P) k = 0;
    }
}

// --- plaintext-interruptor form (breaks causal on the plaintext letter) -------------------------

void intkey_encrypt_ptint(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int base, int interruptor) {

    int k = 0;
    for (int i = 0; i < plaintext_len; i++) {
        encrypted[i] = progkey_base_encrypt(plaintext_indices[i], keyword[k], base);
        if (plaintext_indices[i] == interruptor) k = 0;             // reset after the pt interruptor
        else if (++k >= P) k = 0;
    }
}

void intkey_decrypt_ptint(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int base, int interruptor) {

    int k = 0;
    for (int i = 0; i < cipher_len; i++) {
        decrypted[i] = progkey_base_decrypt(cipher_indices[i], keyword[k], base);
        if (decrypted[i] == interruptor) k = 0;                     // same trigger the encryptor used
        else if (++k >= P) k = 0;
    }
}

// --- ciphertext-interruptor form (breaks causal on the ciphertext letter) -----------------------

void intkey_encrypt_ctint(int encrypted[], int plaintext_indices[], int plaintext_len,
    int keyword[], int P, int base, int interruptor) {

    int k = 0;
    for (int i = 0; i < plaintext_len; i++) {
        encrypted[i] = progkey_base_encrypt(plaintext_indices[i], keyword[k], base);
        if (encrypted[i] == interruptor) k = 0;                     // reset after the ct interruptor
        else if (++k >= P) k = 0;
    }
}

void intkey_decrypt_ctint(int decrypted[], int cipher_indices[], int cipher_len,
    int keyword[], int P, int base, int interruptor) {

    int k = 0;
    for (int i = 0; i < cipher_len; i++) {
        decrypted[i] = progkey_base_decrypt(cipher_indices[i], keyword[k], base);
        if (cipher_indices[i] == interruptor) k = 0;                // ciphertext is known: same k as encrypt
        else if (++k >= P) k = 0;
    }
}

// is_break[i]==1 iff cipher[i-1]==interruptor (reset AFTER the interruptor ciphertext letter).
// Keyword-independent -- the solver builds this once per candidate interruptor, then anneals the
// keyword over intkey_decrypt_mask (equivalent to intkey_decrypt_ctint but with the mask cached).
void intkey_build_mask_ct(int is_break[], const int cipher_indices[], int cipher_len, int interruptor) {
    for (int i = 0; i < cipher_len; i++)
        is_break[i] = (i > 0 && cipher_indices[i - 1] == interruptor) ? 1 : 0;
}

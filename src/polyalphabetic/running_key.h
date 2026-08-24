#ifndef RUNNING_KEY_H
#define RUNNING_KEY_H
#include "colossus.h"

//
// Running Key cipher primitive (ACA type running-key).
//
// A Running Key is a Vigenere-family cipher whose key is a long running TEXT (as
// long as the message), not a repeating short cycleword: character i is enciphered
// against key_stream[i] with no modular wrap. In the ACA convention the plaintext is
// one continuous passage of length 2N split in half -- the FIRST half is the running
// key, the SECOND half is the plaintext, and enciphering the second against the first
// gives the length-N ciphertext. The tableau may be any of the "periodics"; the four
// letter-keyed families are supported here.
//
// Sign conventions mirror autokey.c one-for-one (same tableaux), so
// running_key_decrypt(running_key_encrypt(P)) == P for every family.
//

#define RK_VIGENERE 0   // C = (P + K) mod 26 ;  P = (C - K) mod 26
#define RK_BEAUFORT 1   // C = (K - P) mod 26 ;  P = (K - C) mod 26   (reciprocal)
#define RK_VARIANT  2   // C = (P - K) mod 26 ;  P = (C + K) mod 26
#define RK_PORTA    3   // reciprocal half-alphabet, shift = K/2      (P = enc = dec)
#define RK_N_FAMILIES 4

// Single-character encode/decode against one key letter (exposed for the solver's
// beam warm start and incremental fast path).
int rk_encode_char(int p, int k, int family);
int rk_decode_char(int c, int k, int family);

// The key letter that makes rk_decode_char(c, k, family) == p (crib anchoring): given a
// known plaintext letter p at a ciphertext letter c, the running-key letter is forced.
// (For Porta the key is only determined up to its /2 fold; the even representative is
// returned, which decodes identically.)
int rk_key_from_pt(int c, int p, int family);

// Whole-stream encode/decode: key_stream[i] enciphers position i (running, no wrap).
// key_stream must hold at least `len` letters.
void running_key_encrypt(int out[], const int plain[], int len,
                         const int key_stream[], int family);
void running_key_decrypt(int out[], const int cipher[], int len,
                         const int key_stream[], int family);

const char *rk_family_name(int family);

#endif

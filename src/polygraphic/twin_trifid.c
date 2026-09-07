//
//  Twin Trifid cipher primitives (ACA).
//
//  The Trifid analogue of Twin Bifid: TWO ordinary Trifid messages enciphered under the
//  SAME keyed 3x3x3 cube but with DIFFERENT periods, the two plaintexts sharing a common
//  phrase (the ACA con guarantees a >= 16-letter repeat). The shared cube decrypts BOTH
//  messages, so a solver anneals a single cube while scoring both decrypts jointly (~2x the
//  n-gram signal of a lone Trifid -- see twin_trifid_solver.c). The cipher math is plain
//  Trifid applied twice, so these are thin wrappers over trifid_encrypt / trifid_decrypt
//  (which own the keyed-cube convention and the thread-local coordinate-stream scratch),
//  keeping the twin and single-message code from ever drifting.
//
//  twin_trifid_decrypt writes the two plaintexts CONCATENATED into out[] (message 1 in
//  out[0..n1-1], message 2 in out[n1..n1+n2-1]) so the caller can n-gram-score the pair as
//  one buffer. The messages need not be the same length.
//

#include "colossus.h"

// Encipher the two plaintexts under one shared cube at two periods (out1/out2 separate).
void twin_trifid_encrypt(const int plain1[], int n1, const int plain2[], int n2,
                         const int cube[], int side, int period1, int period2,
                         int out1[], int out2[]) {
    trifid_encrypt(plain1, n1, cube, side, period1, out1);
    trifid_encrypt(plain2, n2, cube, side, period2, out2);
}

// Decipher both ciphertexts under one shared cube at two periods into a single concatenated
// plaintext buffer: out[0..n1-1] = decrypt(cipher1), out[n1..n1+n2-1] = decrypt(cipher2).
// trifid_decrypt's thread-local stream scratch is reused sequentially (each call completes
// before the next), so this stays thread-safe.
void twin_trifid_decrypt(const int cipher1[], int n1, const int cipher2[], int n2,
                         const int cube[], int side, int period1, int period2,
                         int out[]) {
    trifid_decrypt(cipher1, n1, cube, side, period1, out);
    trifid_decrypt(cipher2, n2, cube, side, period2, out + n1);
}

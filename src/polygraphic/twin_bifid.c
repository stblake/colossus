//
//  Twin Bifid cipher primitives (ACA).
//
//  A Twin Bifid is TWO ordinary Bifid messages enciphered under the SAME keyed Polybius
//  square but with DIFFERENT periods; the two plaintexts share a common phrase (the ACA
//  con guarantees a >= 18-letter repeat). Cryptanalytically the shared square is the whole
//  point: one square decrypts BOTH messages, so a solver anneals a single square while
//  scoring both decrypts jointly (~2x the n-gram signal of a lone Bifid -- see
//  twin_bifid_solver.c). The cipher math itself is just plain Bifid applied twice, so these
//  primitives are thin wrappers over bifid_encrypt / bifid_decrypt (which own the keyed-
//  square convention and the thread-local coordinate-stream scratch); this keeps the twin
//  and the single-message code from ever drifting.
//
//  twin_bifid_decrypt writes the two plaintexts CONCATENATED into out[] (message 1 in
//  out[0..n1-1], message 2 in out[n1..n1+n2-1]) so the caller can n-gram-score the pair as
//  one buffer. The messages need not be the same length.
//

#include "colossus.h"

// Encipher the two plaintexts under one shared square at two periods (out1/out2 separate).
void twin_bifid_encrypt(const int plain1[], int n1, const int plain2[], int n2,
                        const int grid[], int side, int period1, int period2,
                        int out1[], int out2[]) {
    bifid_encrypt(plain1, n1, grid, side, period1, out1);
    bifid_encrypt(plain2, n2, grid, side, period2, out2);
}

// Decipher both ciphertexts under one shared square at two periods into a single
// concatenated plaintext buffer: out[0..n1-1] = decrypt(cipher1), out[n1..n1+n2-1] =
// decrypt(cipher2). bifid_decrypt's thread-local stream scratch is reused sequentially
// (each call completes before the next), so this stays thread-safe.
void twin_bifid_decrypt(const int cipher1[], int n1, const int cipher2[], int n2,
                        const int grid[], int side, int period1, int period2,
                        int out[]) {
    bifid_decrypt(cipher1, n1, grid, side, period1, out);
    bifid_decrypt(cipher2, n2, grid, side, period2, out + n1);
}

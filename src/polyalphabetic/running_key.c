//
// Running Key Cipher
//
// Vigenere-family cipher with a running-TEXT key (key length == message length). See
// running_key.h for the ACA construction and the family sign conventions. The four
// families reuse the exact tableau math of autokey.c / beaufort.c / porta.c, applied
// per position against key_stream[i] (no cycleword wrap).
//

#include "running_key.h"

// Encipher one plaintext letter p against key letter k.
int rk_encode_char(int p, int k, int family) {
    switch (family) {
    case RK_BEAUFORT:                                   // C = K - P (mod 26)
        return (k - p + g_alpha) % g_alpha;
    case RK_VARIANT:                                    // C = P - K (mod 26)
        return (p - k + g_alpha) % g_alpha;
    case RK_PORTA: {                                    // reciprocal, shift = K/2
        int shift = k / 2;
        if (p < 13) return (p + shift) % 13 + 13;
        return (p - 13 - shift + ALPHABET_SIZE) % 13;
    }
    case RK_VIGENERE:
    default:                                            // C = P + K (mod 26)
        return (p + k) % g_alpha;
    }
}

// Recover one plaintext letter from ciphertext letter c against key letter k.
int rk_decode_char(int c, int k, int family) {
    switch (family) {
    case RK_BEAUFORT:                                   // P = K - C (mod 26)
        return (k - c + g_alpha) % g_alpha;
    case RK_VARIANT:                                    // P = C + K (mod 26)
        return (c + k) % g_alpha;
    case RK_PORTA: {                                    // reciprocal (identical to encode)
        int shift = k / 2;
        if (c < 13) return (c + shift) % 13 + 13;
        return (c - 13 - shift + ALPHABET_SIZE) % 13;
    }
    case RK_VIGENERE:
    default:                                            // P = C - K (mod 26)
        return (c - k + g_alpha) % g_alpha;
    }
}

int rk_key_from_pt(int c, int p, int family) {
    switch (family) {
    case RK_BEAUFORT:                                   // p = k - c  =>  k = p + c
        return (p + c) % g_alpha;
    case RK_VARIANT:                                    // p = c + k  =>  k = p - c
        return (p - c + g_alpha) % g_alpha;
    case RK_PORTA: {                                    // p = porta(c, k), shift = k/2
        int shift = (c < 13) ? (((p - 13 - c) % 13) + 13) % 13
                             : (((c - 13 - p) % 13) + 13) % 13;
        return 2 * shift;                              // even representative of the /2 fold
    }
    case RK_VIGENERE:
    default:                                            // p = c - k  =>  k = c - p
        return (c - p + g_alpha) % g_alpha;
    }
}

void running_key_encrypt(int out[], const int plain[], int len,
                         const int key_stream[], int family) {
    for (int i = 0; i < len; i++)
        out[i] = rk_encode_char(plain[i], key_stream[i], family);
}

void running_key_decrypt(int out[], const int cipher[], int len,
                         const int key_stream[], int family) {
    for (int i = 0; i < len; i++)
        out[i] = rk_decode_char(cipher[i], key_stream[i], family);
}

const char *rk_family_name(int family) {
    switch (family) {
    case RK_BEAUFORT: return "Beaufort";
    case RK_VARIANT:  return "Variant";
    case RK_PORTA:    return "Porta";
    case RK_VIGENERE:
    default:          return "Vigenere";
    }
}

#include "affine.h"

// =====================================================================
//  Affine cipher primitives -- see affine.h
// =====================================================================

// The 12 multipliers a coprime to 26 and their inverses mod 26 (verified: a*inv % 26 == 1).
const int affine_multipliers[AFFINE_N_MULT] = { 1, 3,  5,  7, 9, 11, 15, 17, 19, 21, 23, 25 };
const int affine_inverses[AFFINE_N_MULT]    = { 1, 9, 21, 15, 3, 19,  7, 23, 11,  5, 17, 25 };

int affine_mod_inverse(int a) {
    a = ((a % AFFINE_MOD) + AFFINE_MOD) % AFFINE_MOD;
    for (int i = 1; i < AFFINE_MOD; i++)
        if ((a * i) % AFFINE_MOD == 1) return i;
    return -1;
}

void affine_encrypt_map(int a, int b, int enc[AFFINE_MOD]) {
    a = ((a % AFFINE_MOD) + AFFINE_MOD) % AFFINE_MOD;
    b = ((b % AFFINE_MOD) + AFFINE_MOD) % AFFINE_MOD;
    for (int p = 0; p < AFFINE_MOD; p++) enc[p] = (a * p + b) % AFFINE_MOD;
}

void affine_decrypt_map(int a, int b, int dec[AFFINE_MOD]) {
    int ainv = affine_mod_inverse(a);
    b = ((b % AFFINE_MOD) + AFFINE_MOD) % AFFINE_MOD;
    for (int c = 0; c < AFFINE_MOD; c++)
        dec[c] = (ainv < 0) ? c : (ainv * ((c - b + AFFINE_MOD) % AFFINE_MOD)) % AFFINE_MOD;
}

void affine_encrypt(const int in[], int len, int a, int b, int out[]) {
    int enc[AFFINE_MOD]; affine_encrypt_map(a, b, enc);
    for (int i = 0; i < len; i++)
        out[i] = (in[i] >= 0 && in[i] < AFFINE_MOD) ? enc[in[i]] : in[i];
}

void affine_decrypt(const int in[], int len, int a, int b, int out[]) {
    int dec[AFFINE_MOD]; affine_decrypt_map(a, b, dec);
    for (int i = 0; i < len; i++)
        out[i] = (in[i] >= 0 && in[i] < AFFINE_MOD) ? dec[in[i]] : in[i];
}

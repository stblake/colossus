#include "chaocipher.h"

// =====================================================================
//  Chaocipher primitive (see chaocipher.h for the algorithm description)
// =====================================================================

// Cyclically rotate a[0..25] LEFT by s positions (element at index (i+s) mod 26 -> i).
static void chao_rotate_left(int a[CHAO_N], int s) {
    int tmp[CHAO_N];
    s %= CHAO_N;
    if (s < 0) s += CHAO_N;
    for (int i = 0; i < CHAO_N; i++) tmp[i] = a[(i + s) % CHAO_N];
    for (int i = 0; i < CHAO_N; i++) a[i] = tmp[i];
}

// Permute the LEFT (ciphertext) alphabet. PRECONDITION: the just-used ciphertext letter is
// already at the zenith (index 0). Remove index 1, shift [2..13] down one, reinsert at 13.
static void chao_permute_left(int a[CHAO_N]) {
    int t = a[CHAO_ZENITH + 1];
    for (int i = CHAO_ZENITH + 1; i < CHAO_NADIR; i++) a[i] = a[i + 1];
    a[CHAO_NADIR] = t;
}

// Permute the RIGHT (plaintext) alphabet. PRECONDITION: the just-used plaintext letter is
// already at the zenith (index 0). Rotate one MORE (its successor to the zenith), then
// remove index 2, shift [3..13] down one, reinsert at 13.
static void chao_permute_right(int a[CHAO_N]) {
    chao_rotate_left(a, 1);
    int t = a[CHAO_ZENITH + 2];
    for (int i = CHAO_ZENITH + 2; i < CHAO_NADIR; i++) a[i] = a[i + 1];
    a[CHAO_NADIR] = t;
}

// --- public disk steppers (single source of truth for the permutation math) ---

void chaocipher_step_left(int a[CHAO_N], int pos) {
    chao_rotate_left(a, pos);      // bring the used ciphertext letter to the zenith
    chao_permute_left(a);
}

void chaocipher_step_right(int a[CHAO_N], int pos) {
    chao_rotate_left(a, pos);      // bring the used plaintext letter to the zenith
    chao_permute_right(a);
}

// Encrypt and decrypt share the same per-step stepping (find the known letter in one disk,
// read the answer from the other at the SAME position, then step BOTH disks by that position).
// They differ only in which disk the lookup happens in.

void chaocipher_encrypt(const int plain[], int n,
                        const int left0[CHAO_N], const int right0[CHAO_N], int out[]) {
    int left[CHAO_N], right[CHAO_N];
    for (int k = 0; k < CHAO_N; k++) { left[k] = left0[k]; right[k] = right0[k]; }

    for (int s = 0; s < n; s++) {
        int p = plain[s];
        int i = 0;
        while (i < CHAO_N && right[i] != p) i++;   // plaintext letter position in RIGHT disk
        out[s] = left[i];                          // ciphertext letter at same position in LEFT
        chaocipher_step_left(left, i);
        chaocipher_step_right(right, i);
    }
}

void chaocipher_decrypt(const int cipher[], int n,
                        const int left0[CHAO_N], const int right0[CHAO_N], int out[]) {
    int left[CHAO_N], right[CHAO_N];
    for (int k = 0; k < CHAO_N; k++) { left[k] = left0[k]; right[k] = right0[k]; }

    for (int s = 0; s < n; s++) {
        int c = cipher[s];
        int i = 0;
        while (i < CHAO_N && left[i] != c) i++;    // ciphertext letter position in LEFT disk
        out[s] = right[i];                         // plaintext letter at same position in RIGHT
        chaocipher_step_left(left, i);
        chaocipher_step_right(right, i);
    }
}

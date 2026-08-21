//
//  Primitive known-answer tests for the Affine cipher (type 87). Framework-free;
//  built by `make test` (links utils.c + affine.c). Checks the modular-inverse table,
//  encrypt/decrypt round-trip over every valid key, and a Wikipedia KAT vector.
//

#include "colossus.h"
#include "affine.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

int main(void) {
    init_alphabet(NULL);                            // 26-letter default; populates g_char_to_idx

    // 1. The 12 units mod 26 and their inverses.
    for (int i = 0; i < AFFINE_N_MULT; i++) {
        int a = affine_multipliers[i];
        CHECK(affine_mod_inverse(a) == affine_inverses[i], "affine inverse of %d wrong", a);
        CHECK((a * affine_inverses[i]) % AFFINE_MOD == 1, "affine %d * inv != 1 mod 26", a);
    }
    CHECK(affine_mod_inverse(2) == -1, "affine: 2 must be non-invertible mod 26");
    CHECK(affine_mod_inverse(13) == -1, "affine: 13 must be non-invertible mod 26");

    // 2. encrypt then decrypt is the identity for every coprime a and every b.
    int pt[AFFINE_MOD]; for (int i = 0; i < AFFINE_MOD; i++) pt[i] = i;
    for (int mi = 0; mi < AFFINE_N_MULT; mi++) {
        int a = affine_multipliers[mi];
        for (int b = 0; b < AFFINE_MOD; b++) {
            int ct[AFFINE_MOD], back[AFFINE_MOD];
            affine_encrypt(pt, AFFINE_MOD, a, b, ct);
            affine_decrypt(ct, AFFINE_MOD, a, b, back);
            int ok = 1; for (int i = 0; i < AFFINE_MOD; i++) if (back[i] != pt[i]) ok = 0;
            CHECK(ok, "affine round-trip failed for a=%d b=%d", a, b);
        }
    }

    // 3. KAT (Wikipedia): a=5, b=8 encrypts AFFINECIPHER -> IHHWVCSWFRCP.
    const char *P = "AFFINECIPHER", *C = "IHHWVCSWFRCP";
    int pin[12], cout[12];
    for (int i = 0; i < 12; i++) pin[i] = g_char_to_idx[(unsigned char) P[i]];
    affine_encrypt(pin, 12, 5, 8, cout);
    int ok = 1; for (int i = 0; i < 12; i++) if (cout[i] != g_char_to_idx[(unsigned char) C[i]]) ok = 0;
    CHECK(ok, "affine KAT a=5,b=8 AFFINECIPHER->IHHWVCSWFRCP failed");

    printf("test_affine: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}

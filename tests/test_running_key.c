// Running Key primitive tests (framework-free).
//
// Pins the ACA worked example (RunningKey.pdf) as a known-answer vector and round-trips
// encrypt/decrypt over random key streams x lengths x all four families, plus the crib
// anchoring inverse (rk_key_from_pt).

#include "running_key.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

static int arrays_equal(const int a[], const int b[], int len) {
    for (int i = 0; i < len; i++) if (a[i] != b[i]) return 0;
    return 1;
}

// A..Z string -> alphabet indices; returns length.
static int str_to_pt(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c >= 'A' && c <= 'Z') out[n++] = g_char_to_idx[c];
    }
    return n;
}

// ---- ACA worked example (Vigenere, CT = key + pt) ------------------------------------
static const char *ACA_KEY = "THISCIPHERCANBEUSEDW";
static const char *ACA_PT  = "ITHANYOFTHEPERIODICS";
static const char *ACA_CT  = "BAPSPGDMXYGPRSMIVMFO";

static void test_kat_vigenere(void) {
    int key[64], pt[64], expect_ct[64];
    int nk = str_to_pt(ACA_KEY, key);
    int np = str_to_pt(ACA_PT, pt);
    int nc = str_to_pt(ACA_CT, expect_ct);
    CHECK(nk == 20 && np == 20 && nc == 20, "KAT lengths %d/%d/%d, expected 20", nk, np, nc);

    int ct[64];
    running_key_encrypt(ct, pt, np, key, RK_VIGENERE);
    CHECK(arrays_equal(ct, expect_ct, np), "encode KAT != published ciphertext BAPSPGDMXYGPRSMIVMFO");

    int back[64];
    running_key_decrypt(back, ct, np, key, RK_VIGENERE);
    CHECK(arrays_equal(back, pt, np), "decode KAT != plaintext ITHANYOFTHEPERIODICS");
}

// ---- Round-trip over random key streams x lengths x families -------------------------
static void test_roundtrip(void) {
    int fails = 0;
    for (int fam = 0; fam < RK_N_FAMILIES; fam++) {
        for (int t = 0; t < 4000; t++) {
            int n = rand_int(1, 300);
            int pt[320], key[320], ct[320], back[320];
            for (int i = 0; i < n; i++) { pt[i] = rand_int(0, 26); key[i] = rand_int(0, 26); }
            running_key_encrypt(ct, pt, n, key, fam);
            running_key_decrypt(back, ct, n, key, fam);
            if (!arrays_equal(back, pt, n)) fails++;
        }
    }
    CHECK(fails == 0, "round-trip failed in %d trials across families", fails);
}

// ---- Crib anchoring inverse: rk_key_from_pt makes decode(c,k)==p ---------------------
// Iterate over KEYS (not arbitrary (c,p) pairs): p = decode(c,k) is an ACHIEVABLE pair,
// and rk_key_from_pt(c,p) must return a key that decodes c to p. (Porta maps between
// alphabet halves, so ~half of all (c,p) pairs are unachievable -- a crib is only
// consistent with the ciphertext in the opposite half; those are not valid inputs.)
static void test_key_from_pt(void) {
    int fails = 0;
    for (int fam = 0; fam < RK_N_FAMILIES; fam++) {
        for (int c = 0; c < 26; c++) {
            for (int k = 0; k < 26; k++) {
                int p = rk_decode_char(c, k, fam);
                int k2 = rk_key_from_pt(c, p, fam);
                CHECK(k2 >= 0 && k2 < 26, "rk_key_from_pt out of range (fam %d)", fam);
                if (rk_decode_char(c, k2, fam) != p) fails++;
            }
        }
    }
    CHECK(fails == 0, "rk_key_from_pt inverse failed in %d achievable (c,p) cases", fails);
}

// ---- Family sign conventions match the primitive comments (spot checks) --------------
static void test_family_conventions(void) {
    // Vigenere C=P+K, Variant C=P-K, Beaufort C=K-P; Porta reciprocal (enc==dec).
    CHECK(rk_encode_char(8, 19, RK_VIGENERE) == (8 + 19) % 26, "Vigenere encode");
    CHECK(rk_decode_char(1, 19, RK_VIGENERE) == 8, "Vigenere decode B,T -> I");
    CHECK(rk_encode_char(8, 19, RK_VARIANT) == (8 - 19 + 26) % 26, "Variant encode");
    CHECK(rk_encode_char(8, 19, RK_BEAUFORT) == (19 - 8) % 26, "Beaufort encode");
    // Porta is an involution: encode then encode returns the input.
    for (int p = 0; p < 26; p++)
        for (int k = 0; k < 26; k++)
            CHECK(rk_encode_char(rk_encode_char(p, k, RK_PORTA), k, RK_PORTA) == p,
                  "Porta not reciprocal at p=%d k=%d", p, k);
}

// ---- Edge cases ---------------------------------------------------------------------
static void test_edges(void) {
    int pt[1] = {5}, key[1] = {7}, ct[1], back[1];
    for (int fam = 0; fam < RK_N_FAMILIES; fam++) {
        running_key_encrypt(ct, pt, 1, key, fam);
        running_key_decrypt(back, ct, 1, key, fam);
        CHECK(back[0] == pt[0], "length-1 round-trip failed (fam %d)", fam);
    }
    // Length 0 must be a no-op (no OOB).
    running_key_encrypt(ct, pt, 0, key, RK_VIGENERE);
    CHECK(1, "length-0 no-op");
}

int main(void) {
    seed_rand(20260824u);
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_kat_vigenere();
    test_roundtrip();
    test_key_from_pt();
    test_family_conventions();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

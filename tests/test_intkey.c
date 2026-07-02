//
//  Unit tests for the Interrupted Key cipher primitives (intkey.c).
//
//  Framework-free: build with `make test`, which links this against intkey.c +
//  vigenere.c + beaufort.c + quagmire.c + utils.c. Exits non-zero if any check fails.
//
//  The Interrupted Key cipher is a periodic base cipher (Vigenere / Variant / Beaufort) under a
//  letter keyword whose key index RESETS to the first key letter at break points, then increments
//  (mod P). Where breaks arise is the STRATEGY: after a chosen ciphertext letter (ct), after a
//  chosen plaintext letter (pt), or at supplied positions (mask). The invariants tested:
//
//    1. The ACA worked-example known-answer vector (Vigenere, keyword ORANGE, the word-division
//       break lengths 4,6,2,3,4,3,1,1,5,1,2,3,5) pinning the whole convention via the mask form.
//    2. Encrypt/decrypt round-trips over random keyword x length x period x interruptor for the
//       ct / pt / mask forms and all three bases, incl. ragged final group, P=1, len=1.
//    3. Agreement with an INDEPENDENT reference (raw shift formulas + a hand-rolled reset walk,
//       sharing no code with intkey.c) for all three forms and bases.
//    4. Edge cases: an interruptor that never occurs (pt form -> plain periodic base cipher); an
//       interruptor that occurs at every position (-> monoalphabetic on the first key letter).
//    5. Consistency: intkey_decrypt_ctint == intkey_decrypt_mask over intkey_build_mask_ct, and
//       encrypt_ctint/ptint/mask each invert under the matching decrypt.
//

#include "colossus.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

static int arrays_equal(int a[], int b[], int len) {
    for (int i = 0; i < len; i++) if (a[i] != b[i]) return 0;
    return 1;
}

static void random_text(int P[], int len) {
    for (int i = 0; i < len; i++) P[i] = rand_int(0, ALPHABET_SIZE);
}
static void random_keyword(int kw[], int len) {
    for (int i = 0; i < len; i++) kw[i] = rand_int(0, ALPHABET_SIZE);
}
static void str_to_indices(const char *s, int out[]) {
    for (int i = 0; s[i]; i++) out[i] = s[i] - 'A';
}

static const int BASES[3] = { IK_BASE_VIG, IK_BASE_VAR, IK_BASE_BEAU };
static const char *BASE_NAME[3] = { "vig", "var", "beau" };

// --- Independent reference (raw formulas + hand-rolled reset walk) --------------------
// Deliberately re-derives the base shift math and the key-index walk WITHOUT calling any
// intkey.c / progkey_base_* code, so it is a genuine cross-check.

static int ref_enc(int p, int k, int base) {
    if (base == IK_BASE_VAR)  return ((p - k) % 26 + 26) % 26;
    if (base == IK_BASE_BEAU) return ((k - p) % 26 + 26) % 26;
    return (p + k) % 26;
}
static void ref_encrypt_mask(int out[], int in[], int n, int kw[], int P, int base, const int is_break[]) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (is_break[i]) k = 0;
        out[i] = ref_enc(in[i], kw[k], base);
        k = (k + 1) % P;
    }
}
static void ref_encrypt_ptint(int out[], int in[], int n, int kw[], int P, int base, int intr) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        out[i] = ref_enc(in[i], kw[k], base);
        if (in[i] == intr) k = 0; else k = (k + 1) % P;
    }
}
static void ref_encrypt_ctint(int out[], int in[], int n, int kw[], int P, int base, int intr) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        out[i] = ref_enc(in[i], kw[k], base);
        if (out[i] == intr) k = 0; else k = (k + 1) % P;
    }
}

// A random break mask: groups of length 1..P; first group = P so the whole keyword is used.
static void random_break_mask(int is_break[], int n, int P) {
    for (int i = 0; i < n; i++) is_break[i] = 0;
    int pos = 0, first = 1;
    while (pos < n) {
        is_break[pos] = 1;
        int L = first ? P : rand_int(1, P + 1);
        first = 0;
        pos += L;
    }
}

// --- Known-answer vector (ACA worked example) ---------------------------------------

static void test_intkey_known_answer(void) {
    // Vigenere, keyword ORANGE (P=6), word-division break lengths 4,6,2,3,4,3,1,1,5,1,2,3,5:
    //   "THISCIPHERCANBEUSEDWITHANYOFTHEPERIODICS" -> "HYIFQZPUKVQRBSEIJEQKZTVOBPOSZVSGSIICUIPY".
    const char *pts = "THISCIPHERCANBEUSEDWITHANYOFTHEPERIODICS";
    const char *cts = "HYIFQZPUKVQRBSEIJEQKZTVOBPOSZVSGSIICUIPY";
    int n = 40, P = 6;
    int pt[64], expected[64], kw[8], C[64], back[64];
    str_to_indices(pts, pt);
    str_to_indices(cts, expected);
    str_to_indices("ORANGE", kw);

    int lens[] = {4, 6, 2, 3, 4, 3, 1, 1, 5, 1, 2, 3, 5};
    int is_break[64] = {0};
    int pos = 0;
    for (int g = 0; g < 13; g++) { is_break[pos] = 1; pos += lens[g]; }
    CHECK(pos == n, "intkey ACA KAT break lengths sum to %d, expected %d", pos, n);

    intkey_encrypt_mask(C, pt, n, kw, P, IK_BASE_VIG, is_break);
    CHECK(arrays_equal(C, expected, n), "intkey ACA KAT encrypt mismatch");

    intkey_decrypt_mask(back, C, n, kw, P, IK_BASE_VIG, is_break);
    CHECK(arrays_equal(back, pt, n), "intkey ACA KAT decrypt mismatch");
}

// --- Round-trip stress + reference agreement (all forms, all bases) -----------------

static void test_intkey_roundtrip(void) {
    int lens[]    = {1, 2, 17, 50, 97, 336, 600};
    int periods[] = {1, 3, 5, 7, 10};
    for (int bi = 0; bi < 3; bi++) {
        int base = BASES[bi];
        for (int li = 0; li < 7; li++) {
            int len = lens[li];
            for (int pi = 0; pi < 5; pi++) {
                int P = periods[pi];
                int pt[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH], back[MAX_CIPHER_LENGTH];
                int ref[MAX_CIPHER_LENGTH], kw[MAX_CYCLEWORD_LEN], is_break[MAX_CIPHER_LENGTH];
                random_text(pt, len);
                random_keyword(kw, P);
                int intr = rand_int(0, ALPHABET_SIZE);

                // mask form
                random_break_mask(is_break, len, P);
                intkey_encrypt_mask(C, pt, len, kw, P, base, is_break);
                ref_encrypt_mask(ref, pt, len, kw, P, base, is_break);
                CHECK(arrays_equal(C, ref, len), "intkey mask vs ref base=%s len=%d P=%d", BASE_NAME[bi], len, P);
                intkey_decrypt_mask(back, C, len, kw, P, base, is_break);
                CHECK(arrays_equal(back, pt, len), "intkey mask round-trip base=%s len=%d P=%d", BASE_NAME[bi], len, P);

                // ct form
                intkey_encrypt_ctint(C, pt, len, kw, P, base, intr);
                ref_encrypt_ctint(ref, pt, len, kw, P, base, intr);
                CHECK(arrays_equal(C, ref, len), "intkey ct vs ref base=%s len=%d P=%d", BASE_NAME[bi], len, P);
                intkey_decrypt_ctint(back, C, len, kw, P, base, intr);
                CHECK(arrays_equal(back, pt, len), "intkey ct round-trip base=%s len=%d P=%d", BASE_NAME[bi], len, P);

                // pt form
                intkey_encrypt_ptint(C, pt, len, kw, P, base, intr);
                ref_encrypt_ptint(ref, pt, len, kw, P, base, intr);
                CHECK(arrays_equal(C, ref, len), "intkey pt vs ref base=%s len=%d P=%d", BASE_NAME[bi], len, P);
                intkey_decrypt_ptint(back, C, len, kw, P, base, intr);
                CHECK(arrays_equal(back, pt, len), "intkey pt round-trip base=%s len=%d P=%d", BASE_NAME[bi], len, P);
            }
        }
    }
}

// --- ct mask consistency: decrypt_ctint == decrypt_mask(build_mask_ct) --------------

static void test_intkey_ct_mask_consistency(void) {
    int len = 300, P = 7, base = IK_BASE_VIG;
    int pt[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH], a[MAX_CIPHER_LENGTH], b[MAX_CIPHER_LENGTH];
    int kw[MAX_CYCLEWORD_LEN], mask[MAX_CIPHER_LENGTH];
    random_text(pt, len);
    random_keyword(kw, P);
    int intr = 4;
    intkey_encrypt_ctint(C, pt, len, kw, P, base, intr);
    intkey_decrypt_ctint(a, C, len, kw, P, base, intr);
    intkey_build_mask_ct(mask, C, len, intr);
    intkey_decrypt_mask(b, C, len, kw, P, base, mask);
    CHECK(arrays_equal(a, b, len), "intkey decrypt_ctint != decrypt_mask(build_mask_ct)");
    CHECK(arrays_equal(a, pt, len), "intkey ct consistency round-trip mismatch");
}

// --- Edge cases -------------------------------------------------------------------

static void test_intkey_edge_cases(void) {
    int len = 200, P = 6, base = IK_BASE_VIG;
    int pt[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH], ref[MAX_CIPHER_LENGTH], back[MAX_CIPHER_LENGTH];
    int kw[MAX_CYCLEWORD_LEN];
    random_keyword(kw, P);

    // (a) pt-interruptor that NEVER occurs -> plain periodic base cipher (k = i % P).
    for (int i = 0; i < len; i++) pt[i] = rand_int(0, 25);   // avoid letter 25 ('Z')
    intkey_encrypt_ptint(C, pt, len, kw, P, base, 25);       // interruptor 'Z' never appears
    for (int i = 0; i < len; i++) ref[i] = ref_enc(pt[i], kw[i % P], base);
    CHECK(arrays_equal(C, ref, len), "intkey pt absent-interruptor != plain periodic");
    intkey_decrypt_ptint(back, C, len, kw, P, base, 25);
    CHECK(arrays_equal(back, pt, len), "intkey pt absent-interruptor round-trip");

    // (b) pt-interruptor at EVERY position (constant plaintext == interruptor) -> monoalphabetic on
    //     the first key letter (k resets after every letter).
    for (int i = 0; i < len; i++) pt[i] = 0;                 // all 'A'
    intkey_encrypt_ptint(C, pt, len, kw, P, base, 0);        // interruptor 'A' == every letter
    for (int i = 0; i < len; i++) CHECK(C[i] == ref_enc(0, kw[0], base),
        "intkey pt every-position not monoalphabetic at %d", i);

    // (c) P == 1: key index is always 0 -> monoalphabetic on kw[0], regardless of breaks.
    int is_break[MAX_CIPHER_LENGTH];
    random_text(pt, len);
    random_break_mask(is_break, len, 3);
    intkey_encrypt_mask(C, pt, len, kw, 1, base, is_break);
    for (int i = 0; i < len; i++) CHECK(C[i] == ref_enc(pt[i], kw[0], base),
        "intkey P=1 not monoalphabetic at %d", i);
}

int main(void) {
    seed_rand(20240617u);
    init_alphabet(NULL);            // full 26-letter alphabet, as the real binary does

    test_intkey_known_answer();
    test_intkey_roundtrip();
    test_intkey_ct_mask_consistency();
    test_intkey_edge_cases();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) {
        printf("TESTS FAILED\n");
        return 1;
    }
    printf("ALL TESTS PASSED\n");
    return 0;
}

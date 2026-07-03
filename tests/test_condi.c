//
//  Unit tests for the Condi cipher primitives (condi.c).
//
//  Framework-free: build with `make test`, which links this against condi.c + utils.c. Exits
//  non-zero if any check fails.
//
//  Condi is a plaintext-feedback substitution over a keyed alphabet sigma: the shift for each
//  letter is the 1-indexed position of the PRECEDING plaintext letter in sigma; the first letter
//  uses a starter offset. idx(ct_i) = (idx(pt_i) + idx(pt_{i-1}) + 1) mod 26 for i >= 2; decryption
//  is causal. Invariants tested:
//
//    1. The ACA worked-example known-answer vector (keyword STRANGE, shifted keyed alphabet
//       VWXYZSTRANGEBCDFHIJKLMOPQU, starter 25) pinning the whole convention -- encrypt asserted
//       cell-for-cell and decrypt round-trip.
//    2. Agreement with an INDEPENDENT reference that follows the ACA VERBAL description via a linear
//       search of the alphabet (no shared code with condi.c's inverse-table implementation), over
//       random keyword x starter x length.
//    3. Encrypt/decrypt round-trips over random sigma x starter x length, incl. len=1.
//    4. Edge cases: len=1 (starter-only); starter=0; the offset-wrap self-encipherment (a previous
//       plaintext letter at the last sigma position -> shift 0 -> the next letter maps to itself);
//       a cyclically shifted keyed alphabet.
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
static void str_to_indices(const char *s, int out[]) {
    for (int i = 0; s[i]; i++) out[i] = s[i] - 'A';
}
static void sigma_inverse(const int sigma[], int inv[]) {
    for (int k = 0; k < ALPHABET_SIZE; k++) inv[sigma[k]] = k;
}
// A random keyed alphabet from a random keyword (via the real make_keyed_alphabet), optionally
// cyclically left-shifted (the ACA "the keyed alphabet may or may not be shifted").
static void random_sigma(int sigma[], int shift) {
    char kw[16];
    int L = rand_int(1, 9);
    for (int i = 0; i < L; i++) kw[i] = 'A' + rand_int(0, ALPHABET_SIZE);
    kw[L] = '\0';
    int keyed[ALPHABET_SIZE];
    make_keyed_alphabet(kw, keyed);
    int s = ((shift % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    for (int k = 0; k < ALPHABET_SIZE; k++) sigma[k] = keyed[(k + s) % ALPHABET_SIZE];
}

// --- Independent reference: the ACA verbal description via a LINEAR search of the alphabet -------
// "substitute the plaintext letter by the one found # places further along the alphabet; the new #
// is the position of that plaintext letter." No inverse table, no shared code with condi.c.
static void ref_encrypt(const int *pt, int n, const int *sigma, int starter, int *ct) {
    int off = ((starter % ALPHABET_SIZE) + ALPHABET_SIZE) % ALPHABET_SIZE;
    for (int i = 0; i < n; i++) {
        int pos = -1;
        for (int k = 0; k < ALPHABET_SIZE; k++) if (sigma[k] == pt[i]) { pos = k; break; }
        ct[i] = sigma[(pos + off) % ALPHABET_SIZE];
        off = (pos + 1) % ALPHABET_SIZE;
    }
}

// --- Tests --------------------------------------------------------------------------------------

static void test_known_answer(void) {
    int sigma[ALPHABET_SIZE], inv[ALPHABET_SIZE];
    str_to_indices("VWXYZSTRANGEBCDFHIJKLMOPQU", sigma);
    sigma_inverse(sigma, inv);
    int starter = 25;
    const char *PT = "OURSISAVERYGREENPASTIMETHEWIDEVARIETYOFCIPHERSWEUSECANALLBESOLVEDWITHPENCILANDPAPER";
    const char *CT = "MORCPPDNBKEDJKPMRTDBQCRJPXCKTVBNHUYJGVBYSFDXKCRCESIUOJJYFRQIXIMBVHKQPDNMPSBYJQBTTNK";
    int n = (int) strlen(PT);
    int pt[128], expected[128], out[128];
    str_to_indices(PT, pt);
    str_to_indices(CT, expected);

    condi_encrypt(pt, n, sigma, inv, starter, out);
    CHECK(arrays_equal(out, expected, n), "ACA known-answer encrypt mismatch");

    int dec[128];
    condi_decrypt(expected, n, sigma, inv, starter, dec);
    CHECK(arrays_equal(dec, pt, n), "ACA known-answer decrypt round-trip mismatch");
}

static void test_reference_and_roundtrip(void) {
    int lens[] = { 1, 2, 17, 50, 97, 240, 600 };
    for (int li = 0; li < (int) (sizeof(lens) / sizeof(lens[0])); li++) {
        int n = lens[li];
        for (int trial = 0; trial < 40; trial++) {
            int sigma[ALPHABET_SIZE], inv[ALPHABET_SIZE];
            random_sigma(sigma, rand_int(0, ALPHABET_SIZE));
            sigma_inverse(sigma, inv);
            int starter = rand_int(0, ALPHABET_SIZE);
            int pt[700], ct[700], ref[700], dec[700];
            random_text(pt, n);

            condi_encrypt(pt, n, sigma, inv, starter, ct);
            ref_encrypt(pt, n, sigma, starter, ref);
            CHECK(arrays_equal(ct, ref, n),
                  "encrypt != reference (len=%d starter=%d trial=%d)", n, starter, trial);

            condi_decrypt(ct, n, sigma, inv, starter, dec);
            CHECK(arrays_equal(dec, pt, n),
                  "decrypt(encrypt) != plaintext (len=%d starter=%d trial=%d)", n, starter, trial);
        }
    }
}

static void test_edge_cases(void) {
    int sigma[ALPHABET_SIZE], inv[ALPHABET_SIZE];
    random_sigma(sigma, 0);
    sigma_inverse(sigma, inv);

    // len=1: only the starter is used; ct[0] = sigma[(idx(pt0) + starter) mod 26].
    for (int starter = 0; starter < ALPHABET_SIZE; starter++) {
        int pt = rand_int(0, ALPHABET_SIZE), ct, dec;
        condi_encrypt(&pt, 1, sigma, inv, starter, &ct);
        CHECK(ct == sigma[(inv[pt] + starter) % ALPHABET_SIZE], "len=1 starter=%d wrong", starter);
        condi_decrypt(&ct, 1, sigma, inv, starter, &dec);
        CHECK(dec == pt, "len=1 round-trip starter=%d", starter);
    }

    // Offset-wrap self-encipherment: after a plaintext letter at the LAST sigma position (index 25,
    // 1-indexed position 26 == 0 mod 26), the NEXT plaintext letter maps to itself.
    int last = sigma[ALPHABET_SIZE - 1];             // letter at position 25
    int mid  = sigma[7];
    int pt[2] = { last, mid }, ct[2];
    condi_encrypt(pt, 2, sigma, inv, /*starter*/3, ct);
    CHECK(ct[1] == mid, "offset-wrap self-encipherment failed (ct[1]=%d mid=%d)", ct[1], mid);

    // starter=0: the first letter maps to itself too (off_1 = 0).
    int p0 = rand_int(0, ALPHABET_SIZE), c0;
    condi_encrypt(&p0, 1, sigma, inv, 0, &c0);
    CHECK(c0 == p0, "starter=0 first letter should be self-enciphered");

    // A shifted keyed alphabet still round-trips.
    int sig2[ALPHABET_SIZE], inv2[ALPHABET_SIZE];
    random_sigma(sig2, 13);
    sigma_inverse(sig2, inv2);
    int q[64], qc[64], qd[64];
    random_text(q, 64);
    condi_encrypt(q, 64, sig2, inv2, 9, qc);
    condi_decrypt(qc, 64, sig2, inv2, 9, qd);
    CHECK(arrays_equal(qd, q, 64), "shifted-alphabet round-trip failed");
}

int main(void) {
    seed_rand(20240703u);
    init_alphabet(NULL);
    if (g_alpha != ALPHABET_SIZE) { printf("FAIL: alphabet is %d, need 26\n", g_alpha); return 1; }

    test_known_answer();
    test_reference_and_roundtrip();
    test_edge_cases();

    printf("%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

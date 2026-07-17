//
//  Unit tests for the Sequence Transposition primitives (pi_from_keyword / encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against sequence_transposition.c +
//  gromark.c + utils.c. Exits non-zero if any check fails.
//
//  Strategy: the ACA worked example anchors two KATs -- a KEYWORD KAT (GUMMYBEARS -> the read
//  order 4 9 5 6 0 2 3 1 7 8) and an ENCODE KAT (THE EARLY BIRD GETS THE WORM under keyword
//  GUMMYBEARS, primer 69315 -> ciphertext YHOMARTBDETHIGWLRESEERT, plus its exact decrypt).
//  Then encode->decode round-trips == identity over random read orders x primers x lengths, a
//  structural check that pi_from_keyword always yields a permutation of 0..9, and the edge
//  cases (single letter, an empty column, a 2-digit primer).
//

#include "colossus.h"

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

// A..Z string -> alphabet indices (non-letters skipped).
static int str_to_pt(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c >= 'A' && c <= 'Z') out[n++] = g_char_to_idx[c];
    }
    return n;
}

// The ACA worked example.
static const char *ACA_PT_STR = "THEEARLYBIRDGETSTHEWORM";
static const char *ACA_KEYWORD = "GUMMYBEARS";
static const int   ACA_PRIMER[5] = {6, 9, 3, 1, 5};
static const char *ACA_CT_STR = "YHOMARTBDETHIGWLRESEERT";
static const int   ACA_PI[10] = {4, 9, 5, 6, 0, 2, 3, 1, 7, 8};

static void test_kat_keyword(void) {
    int pi[10];
    sequence_transposition_pi_from_keyword(ACA_KEYWORD, pi);
    CHECK(arrays_equal(pi, ACA_PI, 10),
          "keyword KAT: GUMMYBEARS -> read order mismatch");
}

static void test_kat_encode_decode(void) {
    int pt[64], expect_ct[64];
    int n = str_to_pt(ACA_PT_STR, pt);
    int m = str_to_pt(ACA_CT_STR, expect_ct);
    CHECK(n == 23 && m == 23, "KAT lengths %d/%d, expected 23/23", n, m);

    int pi[10];
    sequence_transposition_pi_from_keyword(ACA_KEYWORD, pi);

    int ct[64];
    sequence_transposition_encrypt(pt, n, ACA_PRIMER, 5, pi, ct);
    CHECK(arrays_equal(ct, expect_ct, n), "encode KAT != published ciphertext");

    int back[64];
    sequence_transposition_decrypt(ct, n, ACA_PRIMER, 5, pi, back);
    CHECK(arrays_equal(back, pt, n), "decode KAT != plaintext");
}

// pi_from_keyword must always produce a permutation of 0..9 for a 10-letter keyword
// (including repeated letters, resolved by stable position order).
static void test_pi_is_permutation(void) {
    int bad = 0;
    for (int t = 0; t < 2000; t++) {
        char kw[11];
        for (int i = 0; i < 10; i++) kw[i] = 'A' + rand_int(0, 26);
        kw[10] = '\0';
        int pi[10], seen[10] = {0};
        sequence_transposition_pi_from_keyword(kw, pi);
        for (int k = 0; k < 10; k++) if (pi[k] >= 0 && pi[k] < 10) seen[pi[k]]++;
        for (int d = 0; d < 10; d++) if (seen[d] != 1) bad++;
    }
    CHECK(bad == 0, "pi_from_keyword produced a non-permutation in %d slot(s)", bad);
}

static void random_permutation(int pi[], int k) {
    for (int i = 0; i < k; i++) pi[i] = i;
    for (int i = k - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = pi[i]; pi[i] = pi[j]; pi[j] = t; }
}

static void test_roundtrip(void) {
    int fails = 0;
    for (int t = 0; t < 5000; t++) {
        int pi[10]; random_permutation(pi, 10);
        int P = rand_int(2, 6);                     // primer length 2..5
        int primer[8]; for (int i = 0; i < P; i++) primer[i] = rand_int(0, 10);
        int n = rand_int(1, 400);
        int pt[420]; for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int ct[420], back[420];
        sequence_transposition_encrypt(pt, n, primer, P, pi, ct);
        sequence_transposition_decrypt(ct, n, primer, P, pi, back);
        if (!arrays_equal(back, pt, n)) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/5000 trials", fails);
}

// The transposition is a pure permutation, so the ciphertext is always an anagram of the
// plaintext (a structural invariant worth pinning).
static void test_anagram_invariant(void) {
    int fails = 0;
    for (int t = 0; t < 1000; t++) {
        int pi[10]; random_permutation(pi, 10);
        int primer[5]; for (int i = 0; i < 5; i++) primer[i] = rand_int(0, 10);
        int n = rand_int(1, 300);
        int pt[320], ct[320];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        sequence_transposition_encrypt(pt, n, primer, 5, pi, ct);
        int hp[26] = {0}, hc[26] = {0};
        for (int i = 0; i < n; i++) { hp[pt[i]]++; hc[ct[i]]++; }
        for (int c = 0; c < 26; c++) if (hp[c] != hc[c]) { fails++; break; }
    }
    CHECK(fails == 0, "ciphertext was not an anagram of plaintext in %d/1000 trials", fails);
}

static void test_edges(void) {
    int pi[10]; random_permutation(pi, 10);
    int primer[5] = {1, 2, 3, 4, 5};

    // single letter round-trips.
    int one[1] = {7}, c1[1], b1[1];
    sequence_transposition_encrypt(one, 1, primer, 5, pi, c1);
    sequence_transposition_decrypt(c1, 1, primer, 5, pi, b1);
    CHECK(b1[0] == 7, "single-letter round-trip failed (%d)", b1[0]);

    // a constant primer (all zeros) puts every letter in column 0 -> the identity transposition.
    int zprimer[5] = {0, 0, 0, 0, 0};
    int pt[10] = {0,1,2,3,4,5,6,7,8,9}, ct[10], bk[10];
    sequence_transposition_encrypt(pt, 10, zprimer, 5, pi, ct);
    CHECK(arrays_equal(ct, pt, 10), "all-zero primer should be the identity");
    sequence_transposition_decrypt(ct, 10, zprimer, 5, pi, bk);
    CHECK(arrays_equal(bk, pt, 10), "all-zero primer round-trip failed");
}

int main(void) {
    seed_rand(20260717u);
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_kat_keyword();
    test_kat_encode_decode();
    test_pi_is_permutation();
    test_roundtrip();
    test_anagram_invariant();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

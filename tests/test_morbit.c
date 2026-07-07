//
//  Unit tests for the Morbit primitives (encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against morbit.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: the ACA worked example (keyword WISECRACK -> pair->digit 9 5 8 4 2 7 1 3 6;
//  `Once upon a time.` -> CT `27435 88151 28274 65679 378`) anchors two KATs -- a DECODE
//  KAT (the published CT decodes to ONCEUPONATIME, deterministic) and an ENCODE KAT
//  (encoding `ONCE UPON A TIME` with the word spaces reproduces the exact 23-digit CT,
//  exercising the xx word-divider AND the odd-length trailing-X pad). Then encode->decode
//  round-trips == identity over random bijective keys x plaintexts x lengths, and the edge
//  paths (illegal xxx run, a codeword > 4 symbols, a non-bijection key, single letter).
//

#include "morbit.h"

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

// A..Z string -> alphabet indices; a space becomes a WORD-DIVIDER sentinel (-1).
static int str_to_pt(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c == ' ') { out[n++] = -1; continue; }
        if (c < 'A' || c > 'Z') continue;
        out[n++] = g_char_to_idx[c];
    }
    return n;
}

// The ACA key for keyword WISECRACK, as digit(1..9) -> pair index 0..8 (index 0 unused).
// pair->digit is [9,5,8,4,2,7,1,3,6]; inverted, digit->pair is {_,6,4,7,3,1,8,5,2,0}.
static const int ACA_KEY[10] = { 0, 6, 4, 7, 3, 1, 8, 5, 2, 0 };
// The published ciphertext (23 digits, spaces/period stripped).
static const int ACA_CT[23] = {
    2,7,4,3,5, 8,8,1,5,1, 2,8,2,7,4, 6,5,6,7,9, 3,7,8
};

static void test_kat_decode(void) {
    int out[64], nt = 0, nv = 0;
    int n = morbit_decrypt(ACA_CT, 23, ACA_KEY, out, 23, &nt, &nv);
    int expect[16]; str_to_pt("ONCEUPONATIME", expect);
    CHECK(n == 13, "KAT decode length %d, expected 13", n);
    CHECK(n == 13 && arrays_equal(out, expect, 13), "KAT decode != ONCEUPONATIME");
    CHECK(nt == 13 && nv == 13, "KAT decode tokens %d/%d, expected 13/13", nv, nt);
}

static void test_kat_encode(void) {
    // Encode `ONCE UPON A TIME` (with the word spaces) -> must equal the published CT.
    // The stream is 45 symbols; morbit pads one trailing X (-> the final digit 8), so the
    // padding and the xx word divider are both exercised.
    int pt[24]; int n = str_to_pt("ONCE UPON A TIME", pt);
    int cipher[128];
    int clen = morbit_encrypt(pt, n, ACA_KEY, cipher);
    CHECK(clen == 23, "encode length %d, expected 23 (odd-pad + xx divider)", clen);
    CHECK(clen == 23 && arrays_equal(cipher, ACA_CT, 23), "encoded CT != published CT");

    // And it must round-trip to the bare letters.
    int out[64], nt = 0, nv = 0;
    int m = morbit_decrypt(cipher, clen, ACA_KEY, out, 23, &nt, &nv);
    int expect[16]; str_to_pt("ONCEUPONATIME", expect);
    CHECK(m == 13 && arrays_equal(out, expect, 13), "encode->decode != ONCEUPONATIME");
}

// A random valid key: a uniform permutation of the 9 pairs over digits 1..9.
static void random_bijection(int key[10]) {
    int perm[9]; for (int i = 0; i < 9; i++) perm[i] = i;
    for (int i = 8; i > 0; i--) { int j = rand_int(0, i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
    key[0] = 0;
    for (int d = 1; d <= 9; d++) key[d] = perm[d - 1];
}

static void test_roundtrip(void) {
    int fails = 0;
    for (int t = 0; t < 4000; t++) {
        int key[10]; random_bijection(key);
        int n = rand_int(1, 120);
        int pt[130];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);   // letters only (no dividers)
        int cipher[1024];
        int clen = morbit_encrypt(pt, n, key, cipher);
        if (clen < 0) { fails++; continue; }
        int out[1024], nt = 0, nv = 0;
        int m = morbit_decrypt(cipher, clen, key, out, 23, &nt, &nv);
        if (m != n || !arrays_equal(out, pt, n) || nt != n || nv != n) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/4000 trials", fails);
}

static void test_odd_pad(void) {
    // A single letter E (a single dot) is an odd-length stream -> padded to (DOT, X) =
    // one pair -> one digit, and round-trips to E.
    int key[10]; random_bijection(key);
    int e[1] = {4}, ce[16];
    int cl = morbit_encrypt(e, 1, key, ce);
    int oe[16], nt = 0, nv = 0;
    int re = morbit_decrypt(ce, cl, key, oe, 23, &nt, &nv);
    CHECK(cl == 1 && re == 1 && oe[0] == 4 && nt == 1 && nv == 1,
          "odd-pad single-letter E: clen=%d n=%d out=%d", cl, re, re ? oe[0] : -1);
}

static void test_edges(void) {
    int nt = 0, nv = 0;

    // Illegal xxx run: a digit mapping to the (X,X) pair repeated -> a run of X's decodes
    // to no letters and at least one invalid token.
    int kxx[10]; kxx[0] = 0;
    kxx[1] = 8;                                   // digit 1 -> pair (X,X)
    for (int d = 2; d <= 9; d++) kxx[d] = d - 2;  // 0..7 to the rest (a valid bijection)
    int c1[3] = {1, 1, 1}, o1[8];
    int r1 = morbit_decrypt(c1, 3, kxx, o1, 23, &nt, &nv);
    CHECK(r1 == 0 && nt >= 1 && nv == 0, "xxx run: got %d letters, %d/%d valid", r1, nv, nt);

    // A codeword of > 4 symbols is invalid -> filler, not valid. A digit mapping to the
    // (DOT,DOT) pair, three times, is six dots (no X) -> the sentinel flushes one long,
    // invalid codeword.
    int kdd[10]; kdd[0] = 0;
    kdd[1] = 0;                                   // digit 1 -> pair (DOT,DOT)
    for (int d = 2; d <= 9; d++) kdd[d] = d - 1;  // 1..8 to the rest
    int c2[3] = {1, 1, 1}, o2[8];
    int r2 = morbit_decrypt(c2, 3, kdd, o2, 23, &nt, &nv);
    CHECK(r2 == 1 && nv == 0 && o2[0] == 23, "long codeword: got %d letters, %d valid, out[0]=%d",
          r2, nv, r2 ? o2[0] : -1);

    // A non-bijection key (a duplicate pair) cannot encrypt.
    int dup[10]; for (int d = 0; d <= 9; d++) dup[d] = 0;   // every digit -> pair 0
    int pt[4] = {4, 19, 4, 19}, cbuf[64];
    CHECK(morbit_encrypt(pt, 4, dup, cbuf) < 0, "encrypt with a duplicate pair should fail");
}

int main(void) {
    seed_rand(20240707u);
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_kat_decode();
    test_kat_encode();
    test_roundtrip();
    test_odd_pad();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

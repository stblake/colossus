//
//  Unit tests for the Pollux primitives (encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against pollux.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: the ACA worked example (key 1->x 2->- 3->. 4->. 5->x 6->. 7->- 8->- 9->x 0->.;
//  `LUCK HELPS` -> CT `08639 34257 02417 68596 30414 56234 90874 5360`) anchors two KATs --
//  a DECODE KAT (the published CT decodes to LUCKHELPS, deterministic) and a MORSE-STREAM KAT
//  (encoding `LUCK HELPS` with the word space, then mapping the polyphonic digits back through
//  the key, reproduces the exact 39-symbol Morse stream -- exercising the xx word-divider). Then
//  encode->decode round-trips == identity over random keys x plaintexts x lengths, and the edge
//  paths (illegal xxx run, a codeword > 4 symbols, a key missing an element, single letter).
//

#include "pollux.h"

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

// The ACA key, indexed by digit 0..9 (0->. 1->x 2->- 3->. 4->. 5->x 6->. 7->- 8->- 9->x).
static const int ACA_KEY[10] = {
    PX_DOT, PX_X, PX_DASH, PX_DOT, PX_DOT, PX_X, PX_DOT, PX_DASH, PX_DASH, PX_X
};
// The published ciphertext (39 digits, spaces/period stripped).
static const int ACA_CT[39] = {
    0,8,6,3,9, 3,4,2,5,7, 0,2,4,1,7, 6,8,5,9,6, 3,0,4,1,4, 5,6,2,3,4, 9,0,8,7,4, 5,3,6,0
};
// The expected Morse symbol stream for `LUCK HELPS` (DOT=0, DASH=1, X=2).
static const int ACA_STREAM[39] = {
    0,1,0,0, 2, 0,0,1, 2, 1,0,1,0, 2, 1,0,1, 2, 2, 0,0,0,0, 2, 0, 2, 0,1,0,0, 2, 0,1,1,0, 2, 0,0,0
};

static void test_kat_decode(void) {
    int out[64], nt = 0, nv = 0;
    int n = pollux_decrypt(ACA_CT, 39, ACA_KEY, out, 23, &nt, &nv);
    int expect[9]; str_to_pt("LUCKHELPS", expect);
    CHECK(n == 9, "KAT decode length %d, expected 9", n);
    CHECK(n == 9 && arrays_equal(out, expect, 9), "KAT decode != LUCKHELPS");
    CHECK(nt == 9 && nv == 9, "KAT decode tokens %d/%d, expected 9/9", nv, nt);
}

static void test_kat_stream(void) {
    // Map the published CT digits through the key -> must equal the expected Morse stream.
    for (int i = 0; i < 39; i++)
        CHECK(ACA_KEY[ACA_CT[i]] == ACA_STREAM[i], "CT[%d] maps to wrong Morse symbol", i);

    // Encode `LUCK HELPS` (with the word space) and map its polyphonic digits back through the
    // key: the element stream is deterministic (only the digit chosen is random) and must match.
    seed_rand(12345u);
    int pt[16]; int n = str_to_pt("LUCK HELPS", pt);
    int cipher[128];
    int clen = pollux_encrypt(pt, n, ACA_KEY, cipher);
    CHECK(clen == 39, "encode length %d, expected 39 (xx word divider)", clen);
    int ok = (clen == 39);
    for (int i = 0; i < clen && ok; i++) ok = (ACA_KEY[cipher[i]] == ACA_STREAM[i]);
    CHECK(ok, "encoded Morse stream != expected");
    // And it must round-trip to the bare letters.
    int out[64], nt = 0, nv = 0;
    int m = pollux_decrypt(cipher, clen, ACA_KEY, out, 23, &nt, &nv);
    int expect[9]; str_to_pt("LUCKHELPS", expect);
    CHECK(m == 9 && arrays_equal(out, expect, 9), "encode->decode != LUCKHELPS");
}

// A random valid key: forces digit 0/1/2 -> DOT/DASH/X (>= 1 of each), randomises the rest.
static void random_key(int key[10]) {
    key[0] = PX_DOT; key[1] = PX_DASH; key[2] = PX_X;
    for (int d = 3; d < 10; d++) key[d] = rand_int(0, 3);
}

static void test_roundtrip(void) {
    int fails = 0;
    for (int t = 0; t < 4000; t++) {
        int key[10]; random_key(key);
        int n = rand_int(1, 120);
        int pt[130];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int cipher[1024];
        int clen = pollux_encrypt(pt, n, key, cipher);
        if (clen < 0) { fails++; continue; }
        int out[1024], nt = 0, nv = 0;
        int m = pollux_decrypt(cipher, clen, key, out, 23, &nt, &nv);
        if (m != n || !arrays_equal(out, pt, n) || nt != n || nv != n) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/4000 trials", fails);
}

static void test_edges(void) {
    // Illegal xxx run: three X's decode to no letters, one invalid token.
    int key[10]; for (int d = 0; d < 10; d++) key[d] = PX_DOT;
    key[0] = PX_X;
    int c1[3] = {0, 0, 0}, o1[8], nt = 0, nv = 0;
    int r1 = pollux_decrypt(c1, 3, key, o1, 23, &nt, &nv);
    CHECK(r1 == 0 && nt >= 1 && nv == 0, "xxx run: got %d letters, %d/%d valid", r1, nv, nt);

    // A codeword of > 4 symbols is invalid -> filler, not valid.
    for (int d = 0; d < 10; d++) key[d] = PX_DOT;
    key[1] = PX_X;
    int c2[6] = {0, 0, 0, 0, 0, 1}, o2[8];
    int r2 = pollux_decrypt(c2, 6, key, o2, 23, &nt, &nv);
    CHECK(r2 == 1 && nv == 0 && o2[0] == 23, "long codeword: got %d letters, %d valid, out[0]=%d",
          r2, nv, r2 ? o2[0] : -1);

    // A key missing an element cannot encrypt.
    int allx[10]; for (int d = 0; d < 10; d++) allx[d] = PX_X;
    int pt[4] = {4, 19, 4, 19}, cbuf[64];
    CHECK(pollux_encrypt(pt, 4, allx, cbuf) < 0, "encrypt with no dot/dash should fail");

    // Single letter E (a single dot) round-trips.
    int kok[10]; random_key(kok);
    int e[1] = {4}, ce[16];
    int cl = pollux_encrypt(e, 1, kok, ce);
    int oe[16];
    int re = pollux_decrypt(ce, cl, kok, oe, 23, &nt, &nv);
    CHECK(cl == 1 && re == 1 && oe[0] == 4, "single-letter E round-trip");
}

int main(void) {
    seed_rand(20240707u);
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_kat_decode();
    test_kat_stream();
    test_roundtrip();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

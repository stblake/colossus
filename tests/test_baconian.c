//
//  Unit tests for the Baconian primitives (encrypt / decode).
//
//  Framework-free: build with `make test`, which links this against baconian.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: the TWO ACA worked examples (Baconian.pdf) anchor the KATs. The concealment
//  classifier (A-M = 'a', N-Z = 'b') lives in the solver, but it is a couple of lines, so we
//  replicate it here to decode the PUBLISHED cover text end-to-end:
//    Ex1 (per-word): "Now is a good time to attend college. ..."  -> SUCCESS
//    Ex2 (per-letter): "BOWED ASTER PINED JOKED THEIR ..."        -> NOWISAGOODT
//  Plus: the fixed biliteral table both directions, encode->decode round-trips == fold(pt)
//  (J->I, U->V) over random plaintexts x lengths, and the edge paths (impossible >=24 code,
//  trailing partial group, single letter).
//

#include "baconian.h"

#define BAC_FILLER_TEST 23      // arbitrary filler letter index ('X') for invalid groups

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

// A..Z string -> alphabet indices (non-letters skipped). Returns the count.
static int str_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c >= 'A' && c <= 'Z') out[n++] = g_char_to_idx[c];
    }
    return n;
}

// The canonical ACA classifier used by both PDF examples: A-M -> 'a' (0), N-Z -> 'b' (1).
static int classify_amnz(int l) { return (l >= 13) ? 1 : 0; }

// Per-letter concealment: every cover letter is one a/b bit.
static int bits_per_letter(const char *cover, int bits[]) {
    int nb = 0;
    for (int i = 0; cover[i]; i++) {
        int c = toupper((unsigned char) cover[i]);
        if (c >= 'A' && c <= 'Z') bits[nb++] = classify_amnz(g_char_to_idx[c]);
    }
    return nb;
}

// Per-word concealment: the first letter of each whitespace-delimited word is one a/b bit.
static int bits_per_word(const char *cover, int bits[]) {
    int nb = 0, in_word = 0;
    for (int i = 0; cover[i]; i++) {
        int c = toupper((unsigned char) cover[i]);
        int alpha = (c >= 'A' && c <= 'Z');
        if (alpha && !in_word) bits[nb++] = classify_amnz(g_char_to_idx[c]);
        in_word = alpha;
    }
    return nb;
}

// ---- Ex1: per-word cover text -> SUCCESS -------------------------------------
static const char *EX1_COVER =
    "Now is a good time to attend college. School work is a good teacher "
    "and a good builder of character. Every man should be a student and "
    "learn all that there is about a subject.";

// ---- Ex2: per-letter cover text -> NOWISAGOODT -------------------------------
static const char *EX2_COVER =
    "BOWED ASTER PINED JOKED THEIR BLACK HASTE ARRAY INSET CHEST SLING";

static void test_kat_examples(void) {
    int bits[512], out[128], nt = 0, nv = 0;

    int nb1 = bits_per_word(EX1_COVER, bits);
    CHECK(nb1 == 35, "Ex1 per-word bit count %d, expected 35", nb1);
    int m1 = baconian_decode(bits, nb1, out, BAC_FILLER_TEST, &nt, &nv);
    int exp1[16]; int e1 = str_to_idx("SUCCESS", exp1);
    CHECK(m1 == e1 && arrays_equal(out, exp1, e1), "Ex1 decode != SUCCESS");
    CHECK(nt == 7 && nv == 7, "Ex1 tokens %d/%d, expected 7/7", nv, nt);

    int nb2 = bits_per_letter(EX2_COVER, bits);
    CHECK(nb2 == 55, "Ex2 per-letter bit count %d, expected 55", nb2);
    int m2 = baconian_decode(bits, nb2, out, BAC_FILLER_TEST, &nt, &nv);
    int exp2[16]; int e2 = str_to_idx("NOWISAGOODT", exp2);
    CHECK(m2 == e2 && arrays_equal(out, exp2, e2), "Ex2 decode != NOWISAGOODT");
    CHECK(nt == 11 && nv == 11, "Ex2 tokens %d/%d, expected 11/11", nv, nt);
}

// The exact bit pattern for SUCCESS (a=0/b=1), most-significant symbol first, from the table.
static void test_kat_encode(void) {
    int pt[16]; int n = str_to_idx("SUCCESS", pt);
    int bits[128];
    int nb = baconian_encrypt(pt, n, bits);
    CHECK(nb == 35, "encode length %d, expected 35", nb);
    // S=10001 U=10011 C=00010 C=00010 E=00100 S=10001 S=10001
    static const int expect[35] = {
        1,0,0,0,1, 1,0,0,1,1, 0,0,0,1,0, 0,0,0,1,0, 0,0,1,0,0, 1,0,0,0,1, 1,0,0,0,1
    };
    CHECK(arrays_equal(bits, expect, 35), "SUCCESS bit pattern mismatch");
}

// fold J->I (9->8) and V->U (21->20): what the 24-letter decode recovers.
static int fold(int l) { return (l == 9) ? 8 : (l == 21) ? 20 : l; }

static void test_roundtrip(void) {
    int fails = 0;
    for (int t = 0; t < 5000; t++) {
        int n = rand_int(1, 60);
        int pt[80];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int bits[512];
        int nb = baconian_encrypt(pt, n, bits);
        if (nb != 5 * n) { fails++; continue; }
        int out[80], nt = 0, nv = 0;
        int m = baconian_decode(bits, nb, out, BAC_FILLER_TEST, &nt, &nv);
        if (m != n || nt != n || nv != n) { fails++; continue; }
        for (int i = 0; i < n; i++) if (out[i] != fold(pt[i])) { fails++; break; }
    }
    CHECK(fails == 0, "round-trip failed in %d/5000 trials", fails);
}

static void test_edges(void) {
    int out[8], nt = 0, nv = 0;
    // An impossible code (v = 24..31): bits 11000 -> filler, invalid.
    int inval[5] = {1, 1, 0, 0, 0};
    int r1 = baconian_decode(inval, 5, out, BAC_FILLER_TEST, &nt, &nv);
    CHECK(r1 == 1 && nt == 1 && nv == 0 && out[0] == BAC_FILLER_TEST,
          "impossible code: got %d letters, %d/%d valid, out=%d", r1, nv, nt, out[0]);

    // A trailing partial group (< 5 bits) is dropped: 7 bits -> 1 letter.
    int seven[7] = {0,0,0,0,1, 1,0};    // first group = B, then 2 leftover bits
    int r2 = baconian_decode(seven, 7, out, BAC_FILLER_TEST, &nt, &nv);
    CHECK(r2 == 1 && out[0] == 1 && nt == 1 && nv == 1, "partial group: got %d letters, out=%d", r2, out[0]);

    // Single letter round-trips: Z = babbb = 10111 (value 23).
    int z[1] = {25}; int zb[8];
    int nb = baconian_encrypt(z, 1, zb);
    int r3 = baconian_decode(zb, nb, out, BAC_FILLER_TEST, &nt, &nv);
    CHECK(nb == 5 && r3 == 1 && out[0] == 25, "single-letter Z round-trip");
}

int main(void) {
    seed_rand(20260824u);
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_kat_examples();
    test_kat_encode();
    test_roundtrip();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

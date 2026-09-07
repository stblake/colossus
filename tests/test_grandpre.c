//
//  Unit / stress tests for the Grandpre primitive (grandpre.c). Framework-free: a CHECK
//  macro, deterministic seed. Covers the ACA worked-example known-answer vectors (the 8x8
//  LACQUERS square), the polyphonic homophone structure, build validation (bad side / missing
//  letter), and heavy encrypt->decrypt round-trip stress over random squares x plaintexts x
//  sides (6..10). Grandpre uses the full 26-letter alphabet ('A'+i <-> i); the ciphertext is a
//  stream of 2-digit (row,col) codes and decoding is unique while encoding is many-to-one.
//

#include "grandpre.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

// The ACA 8x8 worked example (first column spells LACQUERS).
static const char *const ACA_WORDS[8] = {
    "LADYBUGS", "AZIMUTHS", "CALFSKIN", "QUACKISH",
    "UNJOVIAL", "EVULSION", "ROWDYISM", "SEXTUPLY"
};

static int text_of(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        char c = s[i];
        if (c >= 'A' && c <= 'Z') out[n++] = c - 'A';
        else if (c >= 'a' && c <= 'z') out[n++] = c - 'a';
    }
    return n;
}

// ---- KAT: decode the two ACA example ciphertexts back to their plaintexts ----
static void test_kat_aca(void) {
    GrandpreSquare g;
    CHECK(grandpre_build_from_words(&g, 8, ACA_WORDS) == 0, "ACA 8x8 square build failed");

    // "thefirstcolum" -> 84 27 82 34 56 71 77 26 44 54 64 63 78 (line 1 of "the first column
    // is the keyword"; the trailing N begins line 2 below).
    int codes1[] = {84, 27, 82, 34, 56, 71, 77, 26, 44, 54, 64, 63, 78};
    int pt1[32];
    int n1 = grandpre_decrypt(codes1, 13, &g, pt1);
    int exp1[32]; int en1 = text_of("THEFIRSTCOLUM", exp1);
    int ok1 = (n1 == en1);
    for (int i = 0; i < n1 && ok1; i++) ok1 = (pt1[i] == exp1[i]);
    CHECK(ok1, "ACA decode 1 mismatch (n=%d)", n1);

    // "nisthekeyword" -> 52 66 65 84 27 82 36 61 88 73 54 71 13.
    int codes2[] = {52, 66, 65, 84, 27, 82, 36, 61, 88, 73, 54, 71, 13};
    int pt2[32];
    int n2 = grandpre_decrypt(codes2, 13, &g, pt2);
    int exp2[32]; int en2 = text_of("NISTHEKEYWORD", exp2);
    int ok2 = (n2 == en2);
    for (int i = 0; i < n2 && ok2; i++) ok2 = (pt2[i] == exp2[i]);
    CHECK(ok2, "ACA decode 2 mismatch (n=%d)", n2);
}

// ---- Homophone structure: T occupies two cells (26 and 84); all 26 letters present ----
static void test_homophones(void) {
    GrandpreSquare g;
    grandpre_build_from_words(&g, 8, ACA_WORDS);

    // Every letter has >= 1 code and each code decodes back to that letter.
    int allok = 1;
    for (int L = 0; L < ALPHABET_SIZE; L++) {
        if (g.ncodes[L] < 1) { allok = 0; break; }
        for (int k = 0; k < g.ncodes[L]; k++) {
            int code = g.codes[L][k], dec[1];
            grandpre_decrypt(&code, 1, &g, dec);
            if (dec[0] != L) { allok = 0; break; }
        }
    }
    CHECK(allok, "some letter's codes do not resolve back to it");

    // T is at (2,6)->26 and (8,4)->84 (two homophones).
    int T = 'T' - 'A';
    CHECK(g.ncodes[T] == 2, "T should have 2 codes, got %d", g.ncodes[T]);
    int has26 = 0, has84 = 0;
    for (int k = 0; k < g.ncodes[T]; k++) { if (g.codes[T][k] == 26) has26 = 1; if (g.codes[T][k] == 84) has84 = 1; }
    CHECK(has26 && has84, "T codes should be {26,84}");
}

// ---- Build validation: bad side, missing letter, short word rejected ----
static void test_build_validation(void) {
    GrandpreSquare g;
    int letters[GRANDPRE_MAX_GRID];
    for (int i = 0; i < 36; i++) letters[i] = i % 26;              // 6x6, all 26 present (0..25,0..9)
    CHECK(grandpre_build(&g, 6, letters) == 0, "valid 6x6 build failed");

    CHECK(grandpre_build(&g, 5, letters) == -1, "side 5 must be rejected");
    CHECK(grandpre_build(&g, 11, letters) == -1, "side 11 must be rejected");

    int missing[GRANDPRE_MAX_GRID];
    for (int i = 0; i < 36; i++) missing[i] = i % 25;             // never places letter 25 (Z)
    CHECK(grandpre_build(&g, 6, missing) == -1, "missing-letter square must be rejected");

    const char *const shortw[6] = {"ABCDE", "FGHIJK", "LMNOPQ", "RSTUVW", "XYZABC", "DEFGHI"};
    CHECK(grandpre_build_from_words(&g, 6, shortw) == -1, "row shorter than N must be rejected");
}

// Random N x N square with all 26 letters: seed cells 0..25 with 0..25, rest random, then shuffle.
static void rand_square(int letters[], int N) {
    int M = N * N;
    for (int i = 0; i < 26; i++) letters[i] = i;
    for (int i = 26; i < M; i++) letters[i] = rand_int(0, 26);
    for (int i = M - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = letters[i]; letters[i] = letters[j]; letters[j] = t; }
}

// ---- Round-trip stress: random squares (sides 6..10) x random plaintexts ----
static void test_roundtrip(void) {
    int fails = 0, trials = 0;
    for (int t = 0; t < 6000; t++) {
        int N = rand_int(6, 11);                                  // 6..10
        int letters[GRANDPRE_MAX_GRID];
        rand_square(letters, N);
        GrandpreSquare g;
        CHECK(grandpre_build(&g, N, letters) == 0, "stress build failed at trial %d (N=%d)", t, N);

        int n = rand_int(1, 220), plain[256];
        for (int i = 0; i < n; i++) plain[i] = rand_int(0, 26);

        int codes[256], dec[256];
        int clen = grandpre_encrypt(plain, n, &g, codes);
        int dlen = grandpre_decrypt(codes, clen, &g, dec);
        trials++;

        int ok = (clen == n && dlen == n);
        for (int i = 0; i < n && ok; i++) ok = (dec[i] == plain[i]);
        if (!ok) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/%d trials", fails, trials);
    CHECK(trials == 6000, "not all round-trip trials ran (%d)", trials);
}

int main(void) {
    init_alphabet(NULL);
    seed_rand(20260719u);
    test_kat_aca();
    test_homophones();
    test_build_validation();
    test_roundtrip();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

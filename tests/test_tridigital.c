//
//  Unit / stress tests for the Tridigital primitive (tridigital.c). Framework-free: a CHECK
//  macro, deterministic seed. Covers the ACA worked-example known-answer vector (keywords
//  NOVELCRAFT / DRAGONFLY), the column-membership structure, grid-build validation, heavy
//  encrypt round-trip stress over random grids x random spaced plaintexts, and the space /
//  separator edges. Tridigital uses the full 26-letter alphabet, so letters map directly
//  ('A'+i <-> i) and a word break is the TRI_SPACE sentinel.
//

#include "tridigital.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

// Map a string to an index stream: 'A'..'Z' -> 0..25, ' ' -> TRI_SPACE (word break).
static int text_of(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        char c = s[i];
        if (c >= 'A' && c <= 'Z') out[n++] = c - 'A';
        else if (c >= 'a' && c <= 'z') out[n++] = c - 'a';
        else out[n++] = TRI_SPACE;
    }
    return n;
}

static void rand_perm(int p[], int n) {
    for (int i = 0; i < n; i++) p[i] = i;
    for (int i = n - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = p[i]; p[i] = p[j]; p[j] = t; }
}

// ---- KAT: the ACA worked example (NOVELCRAFT / DRAGONFLY) ----
//        6 7 0 3 5 2 8 1 4 9      column digit-labels (col 9 = separator)
//        D R A G O N F L Y -      row 0
//        B C E H I J K M P -      row 1
//        Q S T U V W X Z - -      row 2
static void test_kat_aca(void) {
    TridigitalGrid g;
    CHECK(tridigital_build_from_keywords(&g, "NOVELCRAFT", "DRAGONFLY") == 0, "ACA build failed");
    CHECK(g.sep_digit == 9, "ACA separator digit %d, expected 9", g.sep_digit);

    // Column digit-labels from NOVELCRAFT: 6 7 0 3 5 2 8 1 4 9.
    int exp_digits[TRI_COLS] = {6, 7, 0, 3, 5, 2, 8, 1, 4, 9};
    int dok = 1;
    for (int c = 0; c < TRI_COLS; c++) if (g.digit_of_col[c] != exp_digits[c]) dok = 0;
    CHECK(dok, "ACA column digit-labels mismatch");

    // Encipher "the ides of march" -> 0 3 0 9 5 6 0 7 9 5 8 9 1 0 7 7 3.
    int plain[64], n = text_of("THE IDES OF MARCH", plain);
    int out[64];
    int clen = tridigital_encrypt(plain, n, &g, out);
    const char *expect = "03095607958910773";
    int elen = (int) strlen(expect);
    CHECK(clen == n, "ACA cipher length %d, expected %d (1:1)", clen, n);
    CHECK(clen == elen, "ACA cipher length %d, expected %d", clen, elen);
    int ok = (clen == elen);
    for (int i = 0; i < clen && ok; i++) ok = (out[i] == expect[i] - '0');
    CHECK(ok, "ACA digit stream mismatch");
}

// ---- Column membership: the (up to 3) letters sharing each digit's column ----
static void test_membership(void) {
    TridigitalGrid g;
    tridigital_build_from_keywords(&g, "NOVELCRAFT", "DRAGONFLY");
    int c[TRI_MAXCOLLET];

    // digit 0 -> column 2 -> {A, E, T}.
    int n0 = tridigital_candidates(0, &g, c);
    CHECK(n0 == 3 && c[0] == 'A' - 'A' && c[1] == 'E' - 'A' && c[2] == 'T' - 'A',
          "digit 0 candidates n=%d {%d,%d,%d}", n0, c[0], c[1], c[2]);

    // digit 6 -> column 0 -> {D, B, Q}.
    int n6 = tridigital_candidates(6, &g, c);
    CHECK(n6 == 3 && c[0] == 'D' - 'A' && c[1] == 'B' - 'A' && c[2] == 'Q' - 'A',
          "digit 6 candidates n=%d {%d,%d,%d}", n6, c[0], c[1], c[2]);

    // digit 4 -> column 8 -> {Y, P} (the 2-letter column).
    int n4 = tridigital_candidates(4, &g, c);
    CHECK(n4 == 2 && c[0] == 'Y' - 'A' && c[1] == 'P' - 'A',
          "digit 4 candidates n=%d {%d,%d}", n4, c[0], c[1]);

    // separator digit 9 -> no letters.
    CHECK(tridigital_candidates(9, &g, c) == 0, "separator digit 9 must have 0 candidates");

    // Every letter enciphers to a digit whose candidate list contains it.
    int allok = 1;
    for (int L = 0; L < TRI_NALPHA; L++) {
        int d = g.digit_of_letter[L];
        int cand[TRI_MAXCOLLET], nc = tridigital_candidates(d, &g, cand), in = 0;
        for (int i = 0; i < nc; i++) if (cand[i] == L) in = 1;
        if (!in || d == g.sep_digit) allok = 0;
    }
    CHECK(allok, "some letter's digit does not resolve back to it");
}

// ---- Build validation: non-permutation inputs and a too-short column keyword rejected ----
static void test_build_validation(void) {
    TridigitalGrid g;
    int good_alpha[TRI_NALPHA]; for (int i = 0; i < TRI_NALPHA; i++) good_alpha[i] = i;
    int good_cols[TRI_COLS] = {6, 7, 0, 3, 5, 2, 8, 1, 4, 9};
    CHECK(tridigital_build_grid(&g, good_alpha, good_cols) == 0, "valid grid build failed");

    int bad_alpha[TRI_NALPHA]; for (int i = 0; i < TRI_NALPHA; i++) bad_alpha[i] = i;
    bad_alpha[0] = 1;                                        // not a permutation (1 twice)
    CHECK(tridigital_build_grid(&g, bad_alpha, good_cols) == -1, "non-perm alphabet must be rejected");

    int bad_cols[TRI_COLS] = {6, 7, 0, 3, 5, 2, 8, 1, 4, 4};  // repeated digit 4
    CHECK(tridigital_build_grid(&g, good_alpha, bad_cols) == -1, "non-perm col digits must be rejected");

    // kw_cols with fewer than 10 distinct letters is rejected.
    CHECK(tridigital_build_from_keywords(&g, "AABBCCDDEE", "SECRET") == -1,
          "column keyword without 10 distinct letters must be rejected");
    CHECK(tridigital_build_from_keywords(&g, "SHORT", "SECRET") == -1,
          "column keyword shorter than 10 letters must be rejected");
}

// ---- Round-trip stress: random grids x random spaced plaintexts ----
static void test_roundtrip(void) {
    int fails = 0, trials = 0;
    for (int t = 0; t < 8000; t++) {
        int keyed[TRI_NALPHA], col_digit[TRI_COLS];
        rand_perm(keyed, TRI_NALPHA);
        rand_perm(col_digit, TRI_COLS);
        TridigitalGrid g;
        CHECK(tridigital_build_grid(&g, keyed, col_digit) == 0, "stress build failed at trial %d", t);

        int n = rand_int(1, 160), plain[200];
        for (int i = 0; i < n; i++) plain[i] = (frand() < 0.18) ? TRI_SPACE : rand_int(0, TRI_NALPHA);

        int cipher[200];
        int clen = tridigital_encrypt(plain, n, &g, cipher);
        trials++;

        int ok = (clen == n);
        for (int i = 0; i < n && ok; i++) {
            if (plain[i] == TRI_SPACE) {
                ok = (cipher[i] == g.sep_digit);
            } else {
                if (cipher[i] != g.digit_of_letter[plain[i]]) { ok = 0; break; }
                int cand[TRI_MAXCOLLET], nc = tridigital_candidates(cipher[i], &g, cand), in = 0;
                for (int j = 0; j < nc; j++) if (cand[j] == plain[i]) in = 1;
                ok = in;
            }
        }
        if (!ok) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/%d trials", fails, trials);
    CHECK(trials == 8000, "not all round-trip trials ran (%d)", trials);
}

// ---- Edges: leading / trailing separators, a lone letter, punctuation -> separator ----
static void test_edges(void) {
    TridigitalGrid g;
    tridigital_build_from_keywords(&g, "NOVELCRAFT", "DRAGONFLY");

    // Leading + trailing spaces both encipher to the separator digit.
    int plain[16], n = text_of(" AB ", plain);          // ' ' 'A' 'B' ' '
    int out[16], clen = tridigital_encrypt(plain, n, &g, out);
    CHECK(clen == 4 && out[0] == g.sep_digit && out[3] == g.sep_digit,
          "leading/trailing separators (clen=%d %d..%d)", clen, out[0], out[3]);
    CHECK(out[1] == g.digit_of_letter['A' - 'A'] && out[2] == g.digit_of_letter['B' - 'A'],
          "AB digits mismatch");

    // A lone letter -> a single digit.
    int one[2], n1 = text_of("A", one), o1[2];
    int c1 = tridigital_encrypt(one, n1, &g, o1);
    CHECK(c1 == 1 && o1[0] == g.digit_of_letter[0], "single letter A (c=%d d=%d)", c1, o1[0]);

    // Punctuation is treated as a word break (separator digit).
    int p2[8], n2 = text_of("A.B", p2), o2[8];           // '.' -> TRI_SPACE via text_of
    int c2 = tridigital_encrypt(p2, n2, &g, o2);
    CHECK(c2 == 3 && o2[1] == g.sep_digit, "punctuation -> separator (c=%d mid=%d)", c2, o2[1]);
}

int main(void) {
    seed_rand(20260715u);
    test_kat_aca();
    test_membership();
    test_build_validation();
    test_roundtrip();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

//
//  Unit / stress tests for the Monome-Dinome primitive (monome_dinome.c). Framework-free:
//  a CHECK macro, deterministic seed. Covers the ACA worked-example known-answer vector
//  (keyword RMASTERTON), a keyed-column KAT, heavy encrypt->decrypt round-trip stress over
//  random boards x random plaintexts, and the invalid-token edges (truncated dinome, a
//  dinome whose second digit is an indicator, a lone monome).
//

#include "monome_dinome.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

// The 24-letter alphabet (A..Z with J->I, Z->Y). Set up like init_alphabet_monome_dinome().
static void init_md_alphabet(void) {
    init_alphabet("JZ");                                 // 24 letters: A..Z minus J and Z
    g_char_to_idx[(int) 'J'] = g_char_to_idx[(int) 'I'];
    g_char_to_idx[(int) 'Z'] = g_char_to_idx[(int) 'Y'];
}

// Map a letter string (A..Y, J->I, Z->Y) to 24-alphabet indices via g_char_to_idx.
static void letters_of(const char *s, int out[]) {
    for (int i = 0; s[i]; i++) out[i] = g_char_to_idx[(int) s[i]];
}

// ---- KAT: the ACA worked example (keyword RMASTERTON -> 6318927054) ----
//        1 8 9 2 7 0 5 4      column labels
//        N O T A R I E S      top row
//     6  B C D F G H K L      indicator 6
//     3  M P Q U V W X Y      indicator 3
static void test_kat_aca(void) {
    int letters[MD_NALPHA];
    letters_of("NOTARIES" "BCDFGHKL" "MPQUVWXY", letters);
    int col_label[MD_COLS] = {1, 8, 9, 2, 7, 0, 5, 4};

    MonomeDinomeBoard b;
    CHECK(monome_dinome_build_board(&b, letters, col_label, 6, 3) == 0, "ACA build failed");

    int plain[64], n = 0;
    letters_of("HIGHFREQUENCYKEYSSHORTENCIPHERTEXT", plain);
    n = 34;

    int out[128];
    int clen = monome_dinome_encrypt(plain, n, &b, out);
    const char *expect = "6006760627539325168346553444608795168038605795359";
    int elen = (int) strlen(expect);
    CHECK(clen == elen, "ACA cipher length %d, expected %d", clen, elen);
    int ok = (clen == elen);
    for (int i = 0; i < clen && ok; i++) ok = (out[i] == expect[i] - '0');
    CHECK(ok, "ACA digit stream mismatch");

    int dec[64], nt = 0, nv = 0;
    int m = monome_dinome_decrypt(out, clen, &b, dec, 23, &nt, &nv);
    CHECK(m == n, "ACA decrypt length %d, expected %d", m, n);
    int dok = (m == n);
    for (int i = 0; i < m && dok; i++) dok = (dec[i] == plain[i]);
    CHECK(dok, "ACA decrypt != plaintext");
    CHECK(nt == n && nv == n, "ACA tokens %d/%d, expected %d/%d valid", nv, nt, n, n);
}

// ---- KAT: a bad column-label set (overlaps an indicator) must be rejected; a genuine
//           keyed set still round-trips. ----
static void test_build_validation(void) {
    int letters[MD_NALPHA];
    letters_of("NOTARIES" "BCDFGHKL" "MPQUVWXY", letters);
    MonomeDinomeBoard b;

    int bad_overlap[MD_COLS] = {1, 8, 9, 2, 7, 0, 5, 6};   // label 6 == indicator 6 -> invalid
    CHECK(monome_dinome_build_board(&b, letters, bad_overlap, 6, 3) == -1,
          "label overlapping an indicator must be rejected");

    int bad_dup[MD_COLS] = {1, 8, 9, 2, 7, 0, 5, 5};       // repeated label 5 -> invalid
    CHECK(monome_dinome_build_board(&b, letters, bad_dup, 6, 3) == -1,
          "repeated column label must be rejected");

    int bad_letters[MD_NALPHA];
    for (int i = 0; i < MD_NALPHA; i++) bad_letters[i] = i;
    bad_letters[0] = 1;                                    // not a permutation
    int good_labels[MD_COLS] = {1, 8, 9, 2, 7, 0, 5, 4};
    CHECK(monome_dinome_build_board(&b, bad_letters, good_labels, 6, 3) == -1,
          "non-permutation letters must be rejected");

    // A different keyed column order + indicator pair still round-trips.
    int labels2[MD_COLS] = {4, 5, 0, 7, 2, 9, 8, 1};       // 8 distinct digits, none == 3 or 9? (9 present)
    CHECK(monome_dinome_build_board(&b, letters, labels2, 3, 6) == 0, "keyed build failed");
    int plain[64], n = 40;
    letters_of("MEETATDAWNBYTHEOLDMILLRACEATNOONTOMORROW", plain);   // no J/Z
    int out[128];
    int clen = monome_dinome_encrypt(plain, n, &b, out);
    CHECK(clen > 0, "keyed encrypt failed");
    int dec[128], nt = 0, nv = 0;
    int m = monome_dinome_decrypt(out, clen, &b, dec, 23, &nt, &nv);
    int ok = (m == n && nv == nt);
    for (int i = 0; i < m && ok; i++) ok = (dec[i] == plain[i]);
    CHECK(ok, "keyed round-trip != plaintext (m=%d n=%d valid=%d/%d)", m, n, nv, nt);
}

// ---- Round-trip stress: random boards x random plaintexts (letters only, all in 0..23) ----
static void rand_perm(int p[], int n) {
    for (int i = 0; i < n; i++) p[i] = i;
    for (int i = n - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = p[i]; p[i] = p[j]; p[j] = t; }
}

static void test_roundtrip(void) {
    int fails = 0, trials = 0;
    for (int t = 0; t < 8000; t++) {
        int letters[MD_NALPHA];
        rand_perm(letters, MD_NALPHA);
        int dperm[10];
        rand_perm(dperm, 10);
        int r0 = dperm[0], r1 = dperm[1];
        int col_label[MD_COLS];
        for (int c = 0; c < MD_COLS; c++) col_label[c] = dperm[2 + c];   // remaining 8 digits

        MonomeDinomeBoard b;
        CHECK(monome_dinome_build_board(&b, letters, col_label, r0, r1) == 0,
              "stress build failed at trial %d", t);

        int n = rand_int(1, 160), plain[200];
        for (int i = 0; i < n; i++) plain[i] = rand_int(0, MD_NALPHA);

        int cipher[512];
        int clen = monome_dinome_encrypt(plain, n, &b, cipher);
        trials++;

        int dec[400], nt = 0, nv = 0;
        int m = monome_dinome_decrypt(cipher, clen, &b, dec, 23, &nt, &nv);
        int ok = (m == n && nv == nt && clen <= 2 * n && clen >= n);
        for (int i = 0; i < m && ok; i++) ok = (dec[i] == plain[i]);
        if (!ok) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/%d trials", fails, trials);
    CHECK(trials == 8000, "not all round-trip trials ran (%d)", trials);
}

// ---- Edges: truncated final dinome, indicator-as-second-digit, a lone monome ----
static void test_edges(void) {
    int letters[MD_NALPHA];
    letters_of("NOTARIES" "BCDFGHKL" "MPQUVWXY", letters);
    int col_label[MD_COLS] = {1, 8, 9, 2, 7, 0, 5, 4};
    MonomeDinomeBoard b;
    monome_dinome_build_board(&b, letters, col_label, 6, 3);

    // A lone indicator digit (6) at end -> a truncated, invalid dinome.
    int trunc[1] = {6}, o1[4], nt = 0, nv = 0;
    int r1 = monome_dinome_decrypt(trunc, 1, &b, o1, 23, &nt, &nv);
    CHECK(r1 == 1 && nt == 1 && nv == 0, "truncated dinome: got %d out, %d/%d valid", r1, nv, nt);

    // A dinome whose second digit is also an indicator (6 then 3) -> invalid on a real board.
    int badpair[2] = {6, 3}, o2[4];
    int r2 = monome_dinome_decrypt(badpair, 2, &b, o2, 23, &nt, &nv);
    CHECK(r2 == 1 && nt == 1 && nv == 0 && o2[0] == 23,
          "indicator-2nd-digit token: got %d out, %d/%d valid, out0=%d", r2, nv, nt, o2[0]);

    // A single top-row letter (N: column label 1) -> one digit, round-trips.
    int one[1]; letters_of("N", one);
    int cbuf[8];
    int cl = monome_dinome_encrypt(one, 1, &b, cbuf);
    int od[8];
    int re = monome_dinome_decrypt(cbuf, cl, &b, od, 23, &nt, &nv);
    CHECK(cl == 1 && cbuf[0] == 1 && re == 1 && od[0] == one[0],
          "single-letter N monome (cl=%d d0=%d re=%d o=%d)", cl, cbuf[0], re, od[0]);

    // A single 3rd-row letter (X: indicator 3, column label 5) -> two digits "35".
    int one2[1]; letters_of("X", one2);
    int cl2 = monome_dinome_encrypt(one2, 1, &b, cbuf);
    CHECK(cl2 == 2 && cbuf[0] == 3 && cbuf[1] == 5, "single-letter X dinome (cl=%d %d%d)",
          cl2, cbuf[0], cbuf[1]);
}

int main(void) {
    seed_rand(20260713u);
    init_md_alphabet();
    CHECK(g_alpha == MD_NALPHA, "alphabet %d, expected %d", g_alpha, MD_NALPHA);

    test_kat_aca();
    test_build_validation();
    test_roundtrip();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

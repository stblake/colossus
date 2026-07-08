//
//  Unit / stress tests for the Straddling Checkerboard primitive (straddling_checkerboard.c).
//  Framework-free: a CHECK macro, deterministic seed. Covers the Wikipedia known-answer
//  vector (fixed 0..9 labels), a keyed-label KAT, a figure-shift round-trip, heavy
//  encrypt->decrypt round-trip stress (letters + digit runs over random boards), and edges.
//

#include "straddling_checkerboard.h"
#include "scoring.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

// Build a 27-symbol board arrangement from a string (A..Z letters, '/' = FIG marker).
static void build_seq(const char *s, int seq[SC_NSYM]) {
    for (int i = 0; i < SC_NSYM; i++)
        seq[i] = (s[i] == '/') ? SC_FIG : (s[i] - 'A');
}

// ---- KAT: the Wikipedia worked example (fixed 0..9 labels, blanks at 2 & 6) ----
static void test_kat_wikipedia(void) {
    int labels[10]; for (int c = 0; c < 10; c++) labels[c] = c;
    // Reading order: top cols 0,1,3,4,5,7,8,9 = E T A O N R I S; row 2 = B C D F G H J K L M;
    // row 6 cols 0..8 = P Q / U V W X Y Z (col 9 = '.' = the NULL cell, not in seq).
    int seq[SC_NSYM];
    build_seq("ETAONRIS" "BCDFGHJKLM" "PQ/UVWXYZ", seq);

    StraddlingBoard b;
    CHECK(straddling_build_board(&b, seq, labels, 2, 6) == 0, "build failed");

    const char *pt = "ATTACKATDAWN";
    int n = 0, plain[32];
    for (int i = 0; pt[i]; i++) plain[n++] = pt[i] - 'A';

    int out[64];
    int clen = straddling_encrypt(plain, n, &b, out);
    const char *expect = "3113212731223655";
    int elen = (int) strlen(expect);
    CHECK(clen == elen, "KAT cipher length %d, expected %d", clen, elen);
    int ok = (clen == elen);
    for (int i = 0; i < clen && ok; i++) ok = (out[i] == expect[i] - '0');
    CHECK(ok, "KAT digit stream mismatch");

    int dec[64], nt = 0, nv = 0;
    int m = straddling_decrypt(out, clen, &b, dec, 23, &nt, &nv);
    CHECK(m == n, "KAT decrypt length %d, expected %d", m, n);
    int dok = (m == n);
    for (int i = 0; i < m && dok; i++) dok = (dec[i] == plain[i]);
    CHECK(dok, "KAT decrypt != plaintext");
    CHECK(nt == n && nv == n, "KAT tokens %d/%d, expected %d/%d", nv, nt, n, n);
}

// ---- KAT: a keyed (scrambled) label permutation still round-trips ----
static void test_kat_keyed_labels(void) {
    int labels[10] = {3, 1, 4, 1, 5, 9, 2, 6, 8, 7};    // NOT a perm on purpose -> must fail
    int seq[SC_NSYM];
    build_seq("ETAONRIS" "BCDFGHJKLM" "PQ/UVWXYZ", seq);
    StraddlingBoard b;
    CHECK(straddling_build_board(&b, seq, labels, 2, 6) == -1, "non-perm labels must be rejected");

    int lab2[10] = {7, 2, 9, 0, 5, 1, 8, 3, 6, 4};      // a genuine permutation of 0..9
    CHECK(straddling_build_board(&b, seq, lab2, lab2[3], lab2[6]) == 0, "keyed-label build failed");
    // indicators are the labels of two chosen physical columns (3 and 6): lab2[3]=0, lab2[6]=8.
    const char *pt = "MEETATDAWNBYTHERIVERX";
    int n = 0, plain[64];
    for (int i = 0; pt[i]; i++) plain[n++] = pt[i] - 'A';
    int out[128];
    int clen = straddling_encrypt(plain, n, &b, out);
    CHECK(clen > 0, "keyed-label encrypt failed");
    int dec[128], nt = 0, nv = 0;
    int m = straddling_decrypt(out, clen, &b, dec, 23, &nt, &nv);
    int ok = (m == n && nv == nt);
    for (int i = 0; i < m && ok; i++) ok = (dec[i] == plain[i]);
    CHECK(ok, "keyed-label round-trip != plaintext (m=%d n=%d valid=%d/%d)", m, n, nv, nt);
}

// ---- Figure-shift: numeric plaintext round-trips ----
static void test_figure_shift(void) {
    int labels[10]; for (int c = 0; c < 10; c++) labels[c] = c;
    int seq[SC_NSYM];
    build_seq("ETAONRIS" "BCDFGHJKLM" "PQ/UVWXYZ", seq);
    StraddlingBoard b;
    straddling_build_board(&b, seq, labels, 2, 6);

    // "MEET AT 0700 HOURS" -> letters + a digit run (0,7,0,0). Digits are 26..35.
    int plain[32], n = 0;
    const char *letters = "MEETAT";
    for (int i = 0; letters[i]; i++) plain[n++] = letters[i] - 'A';
    int digits[] = {0, 7, 0, 0};
    for (int i = 0; i < 4; i++) plain[n++] = 26 + digits[i];
    const char *tail = "HOURS";
    for (int i = 0; tail[i]; i++) plain[n++] = tail[i] - 'A';

    int out[128];
    int clen = straddling_encrypt(plain, n, &b, out);
    CHECK(clen > 0, "figure-shift encrypt failed");
    int dec[128], nt = 0, nv = 0;
    int m = straddling_decrypt(out, clen, &b, dec, 23, &nt, &nv);
    int ok = (m == n && nv == nt);
    for (int i = 0; i < m && ok; i++) ok = (dec[i] == plain[i]);
    CHECK(ok, "figure-shift round-trip != plaintext (m=%d n=%d valid=%d/%d)", m, n, nv, nt);
}

// ---- Round-trip stress: random boards x random plaintexts (letters + digit runs) ----
static void rand_perm(int p[], int n) {
    for (int i = 0; i < n; i++) p[i] = i;
    for (int i = n - 1; i > 0; i--) { int j = rand_int(0, i + 1); int t = p[i]; p[i] = p[j]; p[j] = t; }
}

static void test_roundtrip(void) {
    int fails = 0, trials = 0, skipped = 0;
    for (int t = 0; t < 6000; t++) {
        int seq[SC_NSYM], labels[10];
        rand_perm(seq, SC_NSYM);
        rand_perm(labels, 10);
        int r0 = rand_int(0, 10), r1 = rand_int(0, 10);
        if (r0 == r1) { skipped++; continue; }

        StraddlingBoard b;
        if (straddling_build_board(&b, seq, labels, r0, r1) != 0) { skipped++; continue; }

        int n = rand_int(1, 160), plain[200];
        for (int i = 0; i < n; i++) {
            // ~15% of symbols are digits (exercise figure-shift); rest letters.
            if (rand_int(0, 100) < 15) plain[i] = 26 + rand_int(0, 10);
            else                       plain[i] = rand_int(0, 26);
        }

        int cipher[512];
        int clen = straddling_encrypt(plain, n, &b, cipher);
        if (clen < 0) { skipped++; continue; }             // board can't encipher a needed digit
        trials++;

        int dec[400], nt = 0, nv = 0;
        int m = straddling_decrypt(cipher, clen, &b, dec, 23, &nt, &nv);
        int ok = (m == n && nv == nt);
        for (int i = 0; i < m && ok; i++) ok = (dec[i] == plain[i]);
        if (!ok) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/%d trials (%d skipped)", fails, trials, skipped);
    CHECK(trials > 4000, "too few round-trip trials ran (%d) -- board generation too fragile", trials);
}

// ---- Edges: truncated final token, NULL-cell hit, single symbol ----
static void test_edges(void) {
    int labels[10]; for (int c = 0; c < 10; c++) labels[c] = c;
    int seq[SC_NSYM];
    build_seq("ETAONRIS" "BCDFGHJKLM" "PQ/UVWXYZ", seq);
    StraddlingBoard b;
    straddling_build_board(&b, seq, labels, 2, 6);

    // A lone indicator digit (2) at end -> a truncated, invalid token.
    int trunc[1] = {2}, o1[4], nt = 0, nv = 0;
    int r1 = straddling_decrypt(trunc, 1, &b, o1, 23, &nt, &nv);
    CHECK(r1 == 1 && nt == 1 && nv == 0, "truncated token: got %d out, %d/%d valid", r1, nv, nt);

    // The NULL cell has code (ind[1]=6, labels[9]=9) = "69" -> an invalid token.
    int nullc[2] = {6, 9}, o2[4];
    int r2 = straddling_decrypt(nullc, 2, &b, o2, 23, &nt, &nv);
    CHECK(r2 == 1 && nt == 1 && nv == 0 && o2[0] == 23,
          "NULL-cell token: got %d out, %d/%d valid, out0=%d", r2, nv, nt, o2[0]);

    // A single letter round-trips.
    int one[1] = {4 /*E, a top-row single*/}, cbuf[8];
    int cl = straddling_encrypt(one, 1, &b, cbuf);
    int od[8];
    int re = straddling_decrypt(cbuf, cl, &b, od, 23, &nt, &nv);
    CHECK(cl == 1 && re == 1 && od[0] == 4, "single-letter E round-trip (cl=%d re=%d o=%d)", cl, re, od[0]);
}

int main(void) {
    seed_rand(20240708u);
    init_alphabet(NULL);

    test_kat_wikipedia();
    test_kat_keyed_labels();
    test_figure_shift();
    test_roundtrip();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

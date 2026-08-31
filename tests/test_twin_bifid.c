//
//  Unit tests for the Twin Bifid primitives (twin_bifid_encrypt / twin_bifid_decrypt).
//
//  Framework-free: build with `make test`, which links this against twin_bifid.c + bifid.c +
//  utils.c. Exits non-zero if any check fails.
//
//  A Twin Bifid is just plain Bifid applied to two messages under ONE shared square at two
//  periods, so the primitives are thin wrappers over bifid_encrypt/bifid_decrypt. The tests
//  pin exactly that: (1) an "anchor" check that twin_bifid_encrypt is byte-identical to two
//  independent bifid_encrypt calls (ties the twin to the proven single-message primitive and
//  catches a swapped-period / swapped-message wiring bug); (2) the known Wikipedia Bifid
//  vector as message 1 so a convention change is caught; (3) the CONCATENATION layout of the
//  decrypt (msg1 in out[0..n1-1], msg2 in out[n1..n1+n2-1]); (4) decrypt(encrypt(.)) == (P1,P2)
//  over random squares, lengths and periods (incl. incomplete final blocks, period 1, and
//  period > len); and (5) a 6x6 (36-cell) side-generic round-trip.
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

// A..Z string -> alphabet indices, merging J into I (the 25-letter convention).
static int str_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c == 'J') c = 'I';
        if (c < 'A' || c > 'Z') continue;
        out[n++] = g_char_to_idx[c];
    }
    return n;
}

static void idx_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = index_to_char(a[i]);
    out[len] = '\0';
}

// --- Known-answer vector: message 1 is the Wikipedia Bifid example -------------
//
// square BGWKZQPNDSIOAXEFCLUMTHYVR, FLEEATONCE @ period 10 -> UAEOLWRINS. Twin Bifid must
// reproduce that for message 1 and (anchor) match two separate bifid_encrypt calls.

static void test_twin_bifid_known_answer(void) {
    int grid[PLAYFAIR_GRID];
    int gn = str_to_idx("BGWKZQPNDSIOAXEFCLUMTHYVR", grid);
    CHECK(gn == PLAYFAIR_GRID, "KAT square is not 25 letters (%d)", gn);

    int p1[64], p2[64];
    int n1 = str_to_idx("FLEEATONCE", p1);          // Wikipedia message
    int n2 = str_to_idx("ATTACKATDAWN", p2);        // a second, different message
    int period1 = 10, period2 = 5;

    int o1[64], o2[64];
    twin_bifid_encrypt(p1, n1, p2, n2, grid, 5, period1, period2, o1, o2);

    char cbuf[65]; idx_to_str(o1, n1, cbuf);
    CHECK(strcmp(cbuf, "UAEOLWRINS") == 0,
        "twin bifid msg1 KAT mismatch: got '%s', want 'UAEOLWRINS'", cbuf);

    // Anchor: twin encrypt == two independent bifid encrypts under the same square.
    int ref1[64], ref2[64];
    bifid_encrypt(p1, n1, grid, 5, period1, ref1);
    bifid_encrypt(p2, n2, grid, 5, period2, ref2);
    CHECK(arrays_equal(o1, ref1, n1) && arrays_equal(o2, ref2, n2),
        "twin_bifid_encrypt is not two independent bifid_encrypt calls");

    // Decrypt: concatenated layout, msg1 then msg2, both recovered.
    int dec[128];
    twin_bifid_decrypt(o1, n1, o2, n2, grid, 5, period1, period2, dec);
    CHECK(arrays_equal(dec, p1, n1), "twin decrypt msg1 mismatch");
    CHECK(arrays_equal(dec + n1, p2, n2), "twin decrypt msg2 (concatenation offset n1) mismatch");
}

// --- Round-trip over random squares, lengths and periods -----------------------

static void test_twin_bifid_roundtrip(void) {
    for (int t = 0; t < 4000; t++) {
        int grid[PLAYFAIR_GRID];
        for (int i = 0; i < PLAYFAIR_GRID; i++) grid[i] = i;
        shuffle(grid, PLAYFAIR_GRID);
        int n1 = 1 + rand_int(0, 400), n2 = 1 + rand_int(0, 400);
        int period1 = 1 + rand_int(0, 40), period2 = 1 + rand_int(0, 40);
        int p1[440], p2[440], o1[440], o2[440], dec[880];
        for (int i = 0; i < n1; i++) p1[i] = rand_int(0, PLAYFAIR_GRID);
        for (int i = 0; i < n2; i++) p2[i] = rand_int(0, PLAYFAIR_GRID);
        twin_bifid_encrypt(p1, n1, p2, n2, grid, 5, period1, period2, o1, o2);
        twin_bifid_decrypt(o1, n1, o2, n2, grid, 5, period1, period2, dec);
        CHECK(arrays_equal(dec, p1, n1) && arrays_equal(dec + n1, p2, n2),
            "twin round-trip mismatch (n1=%d p1=%d, n2=%d p2=%d)", n1, period1, n2, period2);
    }
}

// --- Period-1 identity: each message is a monoalphabetic map (letter -> grid[pos]) ---

static void test_twin_bifid_period_one(void) {
    int grid[PLAYFAIR_GRID];
    for (int i = 0; i < PLAYFAIR_GRID; i++) grid[i] = i;
    shuffle(grid, PLAYFAIR_GRID);
    int p1[50], p2[30], o1[50], o2[30], dec[80];
    int n1 = 50, n2 = 30;
    for (int i = 0; i < n1; i++) p1[i] = rand_int(0, PLAYFAIR_GRID);
    for (int i = 0; i < n2; i++) p2[i] = rand_int(0, PLAYFAIR_GRID);
    twin_bifid_encrypt(p1, n1, p2, n2, grid, 5, 1, 1, o1, o2);
    twin_bifid_decrypt(o1, n1, o2, n2, grid, 5, 1, 1, dec);
    CHECK(arrays_equal(dec, p1, n1) && arrays_equal(dec + n1, p2, n2),
        "twin period-1 round-trip mismatch");
}

// --- 6x6 (36-cell) side-generic round-trip -------------------------------------

static void test_twin_bifid_6x6(void) {
    int grid[36];
    for (int i = 0; i < 36; i++) grid[i] = i;
    shuffle(grid, 36);
    int p1[200], p2[200], o1[200], o2[200], dec[400];
    int n1 = 137, n2 = 90;
    for (int i = 0; i < n1; i++) p1[i] = rand_int(0, 36);
    for (int i = 0; i < n2; i++) p2[i] = rand_int(0, 36);
    twin_bifid_encrypt(p1, n1, p2, n2, grid, 6, 7, 9, o1, o2);
    twin_bifid_decrypt(o1, n1, o2, n2, grid, 6, 7, 9, dec);
    CHECK(arrays_equal(dec, p1, n1) && arrays_equal(dec + n1, p2, n2),
        "twin 6x6 round-trip mismatch");
}

int main(void) {
    seed_rand(20240831u);
    init_alphabet("J");                  // 25-letter alphabet (J merged into I)
    CHECK(g_alpha == PLAYFAIR_GRID, "alphabet size %d, expected %d", g_alpha, PLAYFAIR_GRID);
    test_twin_bifid_known_answer();
    test_twin_bifid_roundtrip();
    test_twin_bifid_period_one();
    test_twin_bifid_6x6();
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

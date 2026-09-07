//
//  Unit tests for the Twin Trifid primitives (twin_trifid_encrypt / twin_trifid_decrypt).
//
//  Framework-free: build with `make test`, which links this against twin_trifid.c + trifid.c +
//  utils.c. Exits non-zero if any check fails.
//
//  A Twin Trifid is plain Trifid applied to two messages under ONE shared 3x3x3 cube at two
//  periods, so the primitives are thin wrappers over trifid_encrypt/trifid_decrypt. The tests
//  pin: (1) an anchor check that twin_trifid_encrypt is byte-identical to two independent
//  trifid_encrypt calls; (2) the known Wikipedia / Practical-Cryptography Trifid vectors as the
//  two messages (cube FELIXMARDSTBCGHJKNOPQUVWYZ+, AIDET@5 -> FMJFV, OILEC@5 -> OISSU); (3) the
//  CONCATENATION layout of the decrypt; (4) decrypt(encrypt(.)) == (P1,P2) over random cubes,
//  lengths and periods; and (5) a 2x2x2 side-generic round-trip.
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

// A..Z + '+' string -> alphabet indices (the 27-symbol Trifid alphabet; no J merge).
static int str_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c >= 128) continue;
        int idx = g_char_to_idx[c];
        if (idx < 0) continue;
        out[n++] = idx;
    }
    return n;
}

static void idx_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = index_to_char(a[i]);
    out[len] = '\0';
}

// --- Known-answer vectors: the Wikipedia Trifid worked example as both messages ---

static void test_twin_trifid_known_answer(void) {
    int cube[TRIFID_CELLS];
    int cn = str_to_idx("FELIXMARDSTBCGHJKNOPQUVWYZ+", cube);
    CHECK(cn == TRIFID_CELLS, "KAT cube is not 27 symbols (%d)", cn);

    int p1[64], p2[64];
    int n1 = str_to_idx("AIDET", p1);        // AIDET @ period 5 -> FMJFV
    int n2 = str_to_idx("OILEC", p2);        // OILEC @ period 5 -> OISSU
    int period1 = 5, period2 = 5;

    int o1[64], o2[64];
    twin_trifid_encrypt(p1, n1, p2, n2, cube, TRIFID_SIDE, period1, period2, o1, o2);

    char b1[65], b2[65]; idx_to_str(o1, n1, b1); idx_to_str(o2, n2, b2);
    CHECK(strcmp(b1, "FMJFV") == 0, "twin trifid msg1 KAT mismatch: got '%s', want 'FMJFV'", b1);
    CHECK(strcmp(b2, "OISSU") == 0, "twin trifid msg2 KAT mismatch: got '%s', want 'OISSU'", b2);

    // Anchor: twin encrypt == two independent trifid encrypts under the same cube.
    int ref1[64], ref2[64];
    trifid_encrypt(p1, n1, cube, TRIFID_SIDE, period1, ref1);
    trifid_encrypt(p2, n2, cube, TRIFID_SIDE, period2, ref2);
    CHECK(arrays_equal(o1, ref1, n1) && arrays_equal(o2, ref2, n2),
        "twin_trifid_encrypt is not two independent trifid_encrypt calls");

    // Decrypt: concatenated layout, msg1 then msg2, both recovered.
    int dec[128];
    twin_trifid_decrypt(o1, n1, o2, n2, cube, TRIFID_SIDE, period1, period2, dec);
    CHECK(arrays_equal(dec, p1, n1), "twin trifid decrypt msg1 mismatch");
    CHECK(arrays_equal(dec + n1, p2, n2), "twin trifid decrypt msg2 (offset n1) mismatch");
}

// --- Round-trip over random cubes, lengths and periods -------------------------

static void test_twin_trifid_roundtrip(void) {
    for (int t = 0; t < 4000; t++) {
        int cube[TRIFID_CELLS];
        for (int i = 0; i < TRIFID_CELLS; i++) cube[i] = i;
        shuffle(cube, TRIFID_CELLS);
        int n1 = 1 + rand_int(0, 400), n2 = 1 + rand_int(0, 400);
        int period1 = 1 + rand_int(0, 40), period2 = 1 + rand_int(0, 40);
        int p1[440], p2[440], o1[440], o2[440], dec[880];
        for (int i = 0; i < n1; i++) p1[i] = rand_int(0, TRIFID_CELLS);
        for (int i = 0; i < n2; i++) p2[i] = rand_int(0, TRIFID_CELLS);
        twin_trifid_encrypt(p1, n1, p2, n2, cube, TRIFID_SIDE, period1, period2, o1, o2);
        twin_trifid_decrypt(o1, n1, o2, n2, cube, TRIFID_SIDE, period1, period2, dec);
        CHECK(arrays_equal(dec, p1, n1) && arrays_equal(dec + n1, p2, n2),
            "twin trifid round-trip mismatch (n1=%d p1=%d, n2=%d p2=%d)", n1, period1, n2, period2);
    }
}

// --- 2x2x2 (8-cell) side-generic round-trip ------------------------------------

static void test_twin_trifid_2x2x2(void) {
    int cube[8];
    for (int i = 0; i < 8; i++) cube[i] = i;
    shuffle(cube, 8);
    int p1[120], p2[120], o1[120], o2[120], dec[240];
    int n1 = 83, n2 = 47;
    for (int i = 0; i < n1; i++) p1[i] = rand_int(0, 8);
    for (int i = 0; i < n2; i++) p2[i] = rand_int(0, 8);
    twin_trifid_encrypt(p1, n1, p2, n2, cube, 2, 5, 6, o1, o2);
    twin_trifid_decrypt(o1, n1, o2, n2, cube, 2, 5, 6, dec);
    CHECK(arrays_equal(dec, p1, n1) && arrays_equal(dec + n1, p2, n2),
        "twin trifid 2x2x2 round-trip mismatch");
}

int main(void) {
    seed_rand(20240831u);
    init_alphabet_trifid();              // 27-symbol alphabet (A..Z + '+')
    CHECK(g_alpha == TRIFID_CELLS, "alphabet size %d, expected %d", g_alpha, TRIFID_CELLS);
    test_twin_trifid_known_answer();
    test_twin_trifid_roundtrip();
    test_twin_trifid_2x2x2();
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

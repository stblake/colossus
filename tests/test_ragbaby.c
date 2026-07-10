//
//  Unit tests for the Ragbaby primitives (build / number / encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against ragbaby.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: the ACA worked-example known-answer vector (keyword GROSBEAK -> KA
//  GROSBEAKCDFHILMNPQTUVWYZ; plaintext "word divisions are kept" -> CT "YBBL HNGQDUFGL DEF HFYR")
//  pins the keyed-alphabet build, the word-position numbering, and the mod-24 forward shift cell
//  for cell. Then: decrypt inverts encrypt on the KAT; encrypt<->decrypt round-trips == identity
//  over random keyed alphabets x random word-divided plaintexts; the mod-24 wrap (number 24 == 0,
//  25 == 1) is exercised directly; and the I/J, W/X pairing (plaintext J -> I, X -> W) is checked.
//

#include "colossus.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define RAG_ALPHA 24

static int arrays_equal(const int a[], const int b[], int len) {
    for (int i = 0; i < len; i++) if (a[i] != b[i]) return 0;
    return 1;
}

static void invert_ka(const int ka[], int ka_inv[]) {
    for (int i = 0; i < RAG_ALPHA; i++) ka_inv[ka[i]] = i;
}

static void idx_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = index_to_char(a[i]);
    out[len] = '\0';
}

// A..Z string (spaces/other stripped) -> alphabet indices under the live Ragbaby alphabet.
static int letters_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        unsigned char c = (unsigned char) toupper((unsigned char) s[i]);
        if (c < 'A' || c > 'Z') continue;
        int idx = g_char_to_idx[c];
        if (idx >= 0) out[n++] = idx;
    }
    return n;
}

int main(void) {
    init_alphabet_ragbaby();
    CHECK(g_alpha == RAG_ALPHA, "alphabet size %d (want 24)", g_alpha);
    CHECK(index_to_char(g_char_to_idx[(int) 'J']) == 'I', "J should pair to I");
    CHECK(index_to_char(g_char_to_idx[(int) 'X']) == 'W', "X should pair to W");

    // --- KAT 1: keyed-alphabet build ---------------------------------------------------
    int ka[RAG_ALPHA], ka_inv[RAG_ALPHA];
    ragbaby_build_keyed_alphabet(ka, "GROSBEAK", RAG_ALPHA);
    char ka_str[RAG_ALPHA + 1];
    idx_to_str(ka, RAG_ALPHA, ka_str);
    CHECK(strcmp(ka_str, "GROSBEAKCDFHILMNPQTUVWYZ") == 0,
          "KA build: got %s want GROSBEAKCDFHILMNPQTUVWYZ", ka_str);
    invert_ka(ka, ka_inv);

    // --- KAT 2: numbering + encrypt (the ACA worked example) ---------------------------
    int letters[64], num[64], L = 0;
    ragbaby_number_stream("word divisions are kept", RAG_ALPHA, letters, num, &L);
    CHECK(L == 20, "letter count %d (want 20)", L);
    int expect_num[20] = {1,2,3,4, 2,3,4,5,6,7,8,9,10, 3,4,5, 4,5,6,7};
    CHECK(arrays_equal(num, expect_num, 20), "numbering mismatch");

    int ct[64];
    ragbaby_encrypt(letters, num, L, ka, ka_inv, RAG_ALPHA, ct);
    char ct_str[64];
    idx_to_str(ct, L, ct_str);
    CHECK(strcmp(ct_str, "YBBLHNGQDUFGLDEFHFYR") == 0,
          "encrypt KAT: got %s want YBBLHNGQDUFGLDEFHFYR", ct_str);

    // --- KAT 3: decrypt inverts encrypt on the KAT -------------------------------------
    int back[64];
    ragbaby_decrypt(ct, num, L, ka, ka_inv, RAG_ALPHA, back);
    CHECK(arrays_equal(back, letters, L), "decrypt(encrypt(KAT)) != plaintext");

    // --- I/J, W/X pairing: a J and an X fold to I and W --------------------------------
    {
        int fold[8], fn = letters_to_idx("JX", fold);
        CHECK(fn == 2 && index_to_char(fold[0]) == 'I' && index_to_char(fold[1]) == 'W',
              "J/X should fold to I/W");
        // ragbaby_number_stream folds too (parses "JX" as one word "IW", numbered 1,2)
        int fl[8], fnum[8], flen = 0;
        ragbaby_number_stream("JX", RAG_ALPHA, fl, fnum, &flen);
        CHECK(flen == 2 && index_to_char(fl[0]) == 'I' && index_to_char(fl[1]) == 'W'
              && fnum[0] == 1 && fnum[1] == 2, "number_stream J/X fold + numbering");
    }

    // --- mod-24 wrap: number sequence wraps (24 == 0, 25 == 1); shift 24 == no move ----
    {
        // A single 30-letter word: numbers 1..24 then wrap to 0(=24), 1, ... Check positions.
        int wl[40], wnum[40], wlen = 0;
        ragbaby_number_stream("abcdefghiklmnopqrstuvwyzabcdef", RAG_ALPHA, wl, wnum, &wlen);
        // within-word pos 0..: num = (1 + pos) % 24. pos 23 -> 24 % 24 = 0; pos 24 -> 25 % 24 = 1.
        CHECK(wlen >= 26, "wrap word length %d", wlen);
        CHECK(wnum[0] == 1 && wnum[22] == 23 && wnum[23] == 0 && wnum[24] == 1,
              "wrap numbering: %d %d %d %d", wnum[0], wnum[22], wnum[23], wnum[24]);
        // shift of 24 (== 0 mod 24) leaves a letter unchanged; a fresh KA, one letter.
        int one_in = ka[5], one_out;
        int n24 = 24, p = ka_inv[one_in];
        ragbaby_encrypt(&one_in, &n24, 1, ka, ka_inv, RAG_ALPHA, &one_out);
        CHECK(one_out == ka[(p + 0) % RAG_ALPHA], "shift 24 should equal shift 0");
    }

    // --- round-trips: random KA x random word-divided plaintext ------------------------
    rng_state = 0x9e3779b9u;
    for (int t = 0; t < 4000; t++) {
        // random keyed alphabet (permutation of 0..23)
        int rka[RAG_ALPHA], rinv[RAG_ALPHA];
        for (int i = 0; i < RAG_ALPHA; i++) rka[i] = i;
        for (int i = RAG_ALPHA - 1; i > 0; i--) {
            int j = rand_int(0, i + 1);
            int tmp = rka[i]; rka[i] = rka[j]; rka[j] = tmp;
        }
        invert_ka(rka, rinv);
        // random word-divided plaintext (letters from the live 24-letter alphabet)
        char buf[300]; int bn = 0, nwords = rand_int(1, 12);
        for (int w = 0; w < nwords; w++) {
            if (w) buf[bn++] = ' ';
            int wl = rand_int(1, 14);
            for (int k = 0; k < wl; k++) buf[bn++] = g_idx_to_char_arr[rand_int(0, RAG_ALPHA)];
        }
        buf[bn] = '\0';
        int pl[300], pnum[300], plen = 0;
        ragbaby_number_stream(buf, RAG_ALPHA, pl, pnum, &plen);
        if (plen == 0) continue;
        int c[300], b[300];
        ragbaby_encrypt(pl, pnum, plen, rka, rinv, RAG_ALPHA, c);
        ragbaby_decrypt(c, pnum, plen, rka, rinv, RAG_ALPHA, b);
        CHECK(arrays_equal(b, pl, plen), "round-trip t=%d", t);
    }

    printf("test_ragbaby: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}

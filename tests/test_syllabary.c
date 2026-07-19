//
//  Unit / stress tests for the Syllabary primitive (syllabary.c). Framework-free: a CHECK
//  macro, deterministic seed. Covers the ACA worked-example square (the 100-token "Unknown
//  coordinates, Unknown keysquare" example), decode of all four isolog ciphertexts back to
//  ORDERSRECEIVED, the token-table integrity (100 distinct, all 26 letters + 9 digits + 64
//  syllables + null), and heavy encrypt->decrypt round-trip stress over random squares x
//  plaintexts x tokenization modes (single / greedy / random isolog).
//

#include "syllabary.h"

static int failures = 0, checks = 0;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { \
    failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

// The ACA worked-example square in physical reading order (row 0 top). "" is the null element.
static const char *const PDF_CELLS[SYL_NTOKENS] = {
    "C","3","H","8","AR","M","ING","P","RI","N",
    "CE","A","1","AL","AN","AND","ARE","AS","AT","ATE",
    "ATI","B","2","BE","CA","CO","COM","D","4","DA",
    "DE","E","5","EA","ED","EN","ENT","ER","ERE","ERS",
    "ES","EST","F","6","G","7","HAS","HE","I","9",
    "IN","ION","IS","IT","IVE","J","","K","L","LA",
    "LE","ME","ND","NE","NT","O","OF","ON","OR","OU",
    "Q","R","RA","RE","RED","RES","RO","S","SE","SH",
    "ST","STO","T","TE","TED","TER","TH","THE","THI","THR",
    "TI","TO","U","V","VE","W","WE","X","Y","Z"
};
static const int PDF_ROW_LABELS[SYL_SIDE] = {8, 5, 0, 2, 3, 4, 1, 6, 7, 9};
static const int PDF_COL_LABELS[SYL_SIDE] = {6, 7, 1, 9, 4, 3, 2, 5, 0, 8};

static int text_of(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++)
        if (s[i] >= 'A' && s[i] <= 'Z') out[n++] = s[i] - 'A';
    return n;
}

static int decode_eq(const SyllabarySquare *sq, const int *codes, int nc, const char *expect) {
    int pt[64], n = syllabary_decrypt(codes, nc, sq, pt);
    int exp[64], en = text_of(expect, exp);
    if (n != en) return 0;
    for (int i = 0; i < n; i++) if (pt[i] != exp[i]) return 0;
    return 1;
}

// ---- Token table integrity ----
static void test_token_table(void) {
    // 100 distinct tokens; indices 0..25 are A..Z; the null token round-trips by index.
    int distinct = 1;
    for (int i = 0; i < SYL_NTOKENS; i++)
        for (int j = i + 1; j < SYL_NTOKENS; j++)
            if (strcmp(syllabary_tokens[i], syllabary_tokens[j]) == 0) distinct = 0;
    CHECK(distinct, "syllabary token table has duplicates");
    int letters_ok = 1;
    for (int L = 0; L < 26; L++) { char s[2] = {(char)('A' + L), 0}; if (syllabary_token_index(s) != L) letters_ok = 0; }
    CHECK(letters_ok, "single-letter tokens not at canonical indices 0..25");
    CHECK(syllabary_token_index("") == SYL_NULL_IDX, "null token index wrong");
    CHECK(syllabary_token_index("THE") >= 0 && syllabary_token_index("ING") >= 0, "expected syllables missing");
    CHECK(syllabary_token_index("QQ") == -1, "unknown token should return -1");
}

// ---- KAT: decode the four ACA isolog ciphertexts back to ORDERSRECEIVED ----
static void test_kat_aca(void) {
    SyllabarySquare sq;
    CHECK(syllabary_build_from_strings(&sq, PDF_CELLS, PDF_ROW_LABELS, PDF_COL_LABELS) == 0,
          "ACA square build failed");

    int c1[] = {13,67,5,27,67,65,67,27,86,27,30,99,27,5};     // o r d e r s r e c e i v e d
    int c2[] = {10,26,67,65,69,56,44,5};                      // or de r s re ce ive d
    int c3[] = {10,5,25,65,69,56,30,94,5};                    // or d ers re ce i ve d
    int c4[] = {13,67,26,67,65,67,27,86,27,30,99,24};         // o r de r s r e c e i v ed
    CHECK(decode_eq(&sq, c1, 14, "ORDERSRECEIVED"), "ACA isolog 1 decode mismatch");
    CHECK(decode_eq(&sq, c2,  8, "ORDERSRECEIVED"), "ACA isolog 2 decode mismatch");
    CHECK(decode_eq(&sq, c3,  9, "ORDERSRECEIVED"), "ACA isolog 3 decode mismatch");
    CHECK(decode_eq(&sq, c4, 12, "ORDERSRECEIVED"), "ACA isolog 4 decode mismatch");

    // The single-letter spelling of ORDERSRECEIVED reproduces isolog 1 exactly.
    int pt[32], n = text_of("ORDERSRECEIVED", pt);
    int out[32], m = syllabary_encrypt(pt, n, &sq, SYL_TOK_SINGLE, out);
    int ok = (m == 14);
    for (int i = 0; i < m && ok; i++) ok = (out[i] == c1[i]);
    CHECK(ok, "single-letter encrypt of ORDERSRECEIVED != isolog 1 (m=%d)", m);
}

// ---- Round-trip stress: random squares x random plaintexts x tokenization modes ----
static void test_roundtrip(void) {
    int fails = 0, trials = 0;
    for (int t = 0; t < 4000; t++) {
        SyllabarySquare sq;
        syllabary_build_random(&sq, 100u + t);

        int n = rand_int(1, 200), plain[256];
        for (int i = 0; i < n; i++) plain[i] = rand_int(0, 26);

        int mode = t % 3;                                     // single / greedy / random
        int codes[256], dec[256];
        int clen = syllabary_encrypt(plain, n, &sq, mode, codes);
        int dlen = syllabary_decrypt(codes, clen, &sq, dec);
        trials++;

        int ok = (dlen == n);
        for (int i = 0; i < n && ok; i++) ok = (dec[i] == plain[i]);
        if (!ok) fails++;
    }
    CHECK(fails == 0, "round-trip failed in %d/%d trials", fails, trials);
    CHECK(trials == 4000, "not all round-trip trials ran (%d)", trials);
}

int main(void) {
    init_alphabet(NULL);
    seed_rand(20260719u);
    test_token_table();
    test_kat_aca();
    test_roundtrip();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

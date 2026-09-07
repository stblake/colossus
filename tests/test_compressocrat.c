//
//  Unit tests for the Compressocrat primitives (encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against compressocrat.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: the ACA worked-example known-answer vector (keyed alphabet YZFRACTIONBDEGHJKLMPQSUVWX,
//  a 60-letter plaintext -> a 62-letter ciphertext) pins the fixed {1,2,3} Huffman table, the
//  trigraph rank (r = 9(a-1)+3(b-1)+(c-1)), the padding, and the keyed-alphabet mapping cell for
//  cell; a small HAND-COMPUTED vector (TA -> DS under identity sigma) pins the rank/padding
//  arithmetic directly. Then a structural check of the fixed code (prefix-free, and the trigraph
//  333 unreachable), an INDEPENDENT in-test reference encrypt (fresh code sharing nothing with
//  compressocrat.c) checked against compressocrat_encrypt over random inputs, encrypt->decrypt
//  round-trips == identity over random permutations sigma and lengths (incl. the three padding
//  residues and rare-letter-heavy inputs), and the invalid-token / trailing-partial paths.
//

#include "colossus.h"
#include "compressocrat.h"

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

// A..Z string -> alphabet indices (default 26-letter alphabet, no J->I).
static int str_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c < 'A' || c > 'Z') continue;
        out[n++] = g_char_to_idx[c];
    }
    return n;
}

static void idx_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = index_to_char(a[i]);
    out[len] = '\0';
}

// Keyed alphabet sigma (rank -> letter): keyword letters (dedup) then the rest ascending.
static void build_sigma(const char *kw, int sigma[]) {
    char used[26];
    for (int i = 0; i < 26; i++) used[i] = 0;
    int m = 0;
    for (int i = 0; kw[i]; i++) {
        int c = toupper((unsigned char) kw[i]);
        if (c < 'A' || c > 'Z') continue;
        int l = g_char_to_idx[c];
        if (used[l]) continue;
        used[l] = 1; sigma[m++] = l;
    }
    for (int l = 0; l < 26 && m < 26; l++) if (!used[l]) { used[l] = 1; sigma[m++] = l; }
}

// --- INDEPENDENT reference (fresh code, own copy of the compression table) -------

static const char *const REF_CODE[26] = {
    "13",     "32112",  "1112",   "213",    "31",     "3213",   "32113",  "113",
    "322",    "321112", "11112",  "212",    "2111",   "23",     "22",     "3212",
    "11113",  "323",    "112",    "12",     "1113",   "11111",  "2112",   "321111",
    "2113",   "321113"
};

// Concatenate codes -> {1,2,3} stream (digits 1..3), pad with 1, group into trigraph ranks,
// emit sigma[rank]. Returns ciphertext length; also (optionally) exposes the padded stream.
static int ref_encrypt(const int *pt, int n, const int *sigma, int *out, int *stream, int *slen_out) {
    static int s[6 * 4000 + 8];
    int *sp = stream ? stream : s;
    int slen = 0;
    for (int i = 0; i < n; i++)
        for (const char *c = REF_CODE[pt[i]]; *c; c++) sp[slen++] = *c - '0';
    while (slen % 3 != 0) sp[slen++] = 1;                          // pad with the digit 1
    if (slen_out) *slen_out = slen;
    int clen = slen / 3;
    for (int g = 0; g < clen; g++)
        out[g] = sigma[9 * (sp[3 * g] - 1) + 3 * (sp[3 * g + 1] - 1) + (sp[3 * g + 2] - 1)];
    return clen;
}

// --- ACA worked-example KAT (verified: 60-letter pt <-> 62-letter ct) ------------

static void test_kat_aca(void) {
    int sigma[26];
    // Keyed alphabet read off the example's trigraph digit rows: 111->Y 112->Z ... 332->X.
    int ns = str_to_idx("YZFRACTIONBDEGHJKLMPQSUVWX", sigma);
    CHECK(ns == 26, "ACA KAT: keyed alphabet is %d letters, expected 26", ns);

    const char *PT = "THEKEYWORDMAYBESHIFTEDTOAVOIDWXYZALWAYSENCODINGTHESAMEDIGITS";
    const char *CT = "ROYZPFNGVDNFNXZMROHDCRPEMYHGIZSYNXYQBNDNMJJZGOGXFROYDNTDUSOENN";

    int pt[128], ct[128], enc[128], dec[128];
    int n = str_to_idx(PT, pt);
    int clen = str_to_idx(CT, ct);
    CHECK(n == 60 && clen == 62, "ACA KAT: pt=%d (want 60) ct=%d (want 62)", n, clen);

    int ce = compressocrat_encrypt(pt, n, sigma, enc);
    char ebuf[128]; idx_to_str(enc, ce, ebuf);
    CHECK(ce == clen && arrays_equal(enc, ct, clen), "ACA KAT: encrypt mismatch, got '%s'", ebuf);

    int nt = 0, nv = 0;
    int m = compressocrat_decrypt(ct, clen, sigma, dec, 23, &nt, &nv);
    CHECK(m == n && nt == n && nv == n && arrays_equal(dec, pt, n),
        "ACA KAT: decrypt mismatch (m=%d nt=%d nv=%d, want %d valid)", m, nt, nv, n);
}

// --- small hand-computed KAT (identity sigma): TA -> DS, EN -> TS ------------------

static void test_kat_small(void) {
    int sigma[26];
    for (int i = 0; i < 26; i++) sigma[i] = i;               // identity: ct letter i == rank i

    // T=12, A=13 -> stream 1,2,1,3 -> pad two 1s -> 1,2,1,3,1,1 -> ranks (1,2,1)=3, (3,1,1)=18
    //   -> letters D(3), S(18).
    // E=31, N=23 -> stream 3,1,2,3 -> pad two 1s -> 3,1,2,3,1,1 -> ranks (3,1,2)=19, (3,1,1)=18
    //   -> letters T(19), S(18).
    struct { const char *pt; const char *ct; } kat[] = { { "TA", "DS" }, { "EN", "TS" } };
    for (int t = 0; t < 2; t++) {
        int pt[8], ct[8], back[8]; char cbuf[9];
        int n = str_to_idx(kat[t].pt, pt);
        int clen = compressocrat_encrypt(pt, n, sigma, ct);
        idx_to_str(ct, clen, cbuf);
        CHECK(strcmp(cbuf, kat[t].ct) == 0, "KAT %s: got '%s', want '%s'", kat[t].pt, cbuf, kat[t].ct);
        int nt = 0, nv = 0;
        int m = compressocrat_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
        CHECK(m == n && nt == n && nv == n && arrays_equal(back, pt, n),
            "KAT %s: decrypt round-trip mismatch (m=%d nt=%d nv=%d)", kat[t].pt, m, nt, nv);
    }
}

// --- structural: the fixed code is prefix-free, and trigraph 333 is unreachable ---

static void test_code_structure(void) {
    // Prefix-free: no code is a prefix of another distinct code.
    for (int i = 0; i < 26; i++)
        for (int j = 0; j < 26; j++) {
            if (i == j) continue;
            size_t li = strlen(REF_CODE[i]);
            CHECK(!(li <= strlen(REF_CODE[j]) && strncmp(REF_CODE[i], REF_CODE[j], li) == 0),
                "code %c=%s is a prefix of %c=%s", 'A' + i, REF_CODE[i], 'A' + j, REF_CODE[j]);
        }

    // 333 never occurs in a real {1,2,3} stream (no code contains "33"; after a leading 3 the
    // next digit is always 1 or 2), so encrypt never yields rank 26. Verify on random plaintexts.
    int sigma[26]; for (int i = 0; i < 26; i++) sigma[i] = i;
    for (int t = 0; t < 500; t++) {
        int n = 1 + rand_int(0, 300);
        int pt[320], out[2 * 320 + 8], stream[6 * 320 + 8], slen = 0;
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        ref_encrypt(pt, n, sigma, out, stream, &slen);
        int bad = 0;
        for (int i = 0; i + 2 < slen; i++)
            if (stream[i] == 3 && stream[i + 1] == 3 && stream[i + 2] == 3) { bad = 1; break; }
        CHECK(!bad, "trial %d: the substring 333 appeared in the {1,2,3} stream", t);
    }
}

// --- compressocrat_encrypt vs the independent reference, random inputs ------------

static void test_reference_agreement(void) {
    const char *kws[] = { "SHADOW", "KRYPTOS", "CIPHERMACHINE", "ZEBRA", "" };
    for (int t = 0; t < 3000; t++) {
        int sigma[26];
        build_sigma(kws[t % 5], sigma);
        if ((t % 5) == 4) { for (int i = 0; i < 26; i++) sigma[i] = i; shuffle(sigma, 26); }
        int n = 1 + rand_int(0, 400);
        int pt[420], a[2 * 420 + 8], b[2 * 420 + 8];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int ca = compressocrat_encrypt(pt, n, sigma, a);
        int cb = ref_encrypt(pt, n, sigma, b, NULL, NULL);
        CHECK(ca == cb && arrays_equal(a, b, ca),
            "reference disagreement (trial %d, n=%d): ca=%d cb=%d", t, n, ca, cb);
    }
}

// --- encrypt -> decrypt round-trip == identity over random sigma and lengths ------

static void test_roundtrip(void) {
    for (int t = 0; t < 6000; t++) {
        int sigma[26];
        for (int i = 0; i < 26; i++) sigma[i] = i;
        shuffle(sigma, 26);                                  // any permutation, not just keyed
        int n = 1 + rand_int(0, 400);
        int pt[420], ct[2 * 420 + 8], back[2 * 420 + 8];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int clen = compressocrat_encrypt(pt, n, sigma, ct);
        int nt = 0, nv = 0;
        int m = compressocrat_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
        CHECK(m == n && nt == n && nv == n && arrays_equal(back, pt, n),
            "round-trip mismatch (trial %d, n=%d): m=%d nt=%d nv=%d", t, n, m, nt, nv);
    }
}

// --- edges: every single letter, rare-letter-heavy (longest codes), padding residues ---

static void test_edges(void) {
    int sigma[26];
    for (int i = 0; i < 26; i++) sigma[i] = i;               // identity sigma

    // Every single letter round-trips (exercises each code + the trailing 1-pad drop).
    for (int l = 0; l < 26; l++) {
        int pt[1] = { l }, ct[8], back[8], nt = 0, nv = 0;
        int clen = compressocrat_encrypt(pt, 1, sigma, ct);
        int m = compressocrat_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
        CHECK(m == 1 && back[0] == l && nv == 1, "single-letter %c round-trip failed (m=%d)", 'A' + l, m);
    }

    // Rare-letter-heavy (J/X/Z, the 6-digit codes) and E/T (2-digit) at all padding residues.
    const int hard[] = { 9, 23, 25, 4, 19 };                 // J, X, Z, E, T
    for (int h = 0; h < 5; h++)
        for (int L = 1; L <= 20; L++) {
            int pt[20], ct[128], back[128], nt = 0, nv = 0;
            for (int i = 0; i < L; i++) pt[i] = hard[h];
            int clen = compressocrat_encrypt(pt, L, sigma, ct);
            int m = compressocrat_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
            CHECK(m == L && nv == L && arrays_equal(back, pt, L),
                "all-%c round-trip failed (L=%d, m=%d nv=%d)", 'A' + hard[h], L, m, nv);
        }
}

// --- invalid-token + trailing-partial paths --------------------------------------

static void test_invalid_paths(void) {
    int sigma[26];
    for (int i = 0; i < 26; i++) sigma[i] = i;               // identity: ct letter r -> rank r

    // Rank 24 = trigraph (3,3,1) = digits "331". "331331" starts with the dead prefix "33"
    // (no code contains 33), so 6 digits pass with no match -> a single invalid token (filler).
    int cipher[2] = { 24, 24 }, out[8], nt = 0, nv = 0;
    int m = compressocrat_decrypt(cipher, 2, sigma, out, 23, &nt, &nv);
    CHECK(m == 1 && nt == 1 && nv == 0 && out[0] == 23,
        "invalid-token path wrong: m=%d nt=%d nv=%d out0=%d", m, nt, nv, out[0]);

    // Rank 5 = (1,2,3) = "123": parses "12"->T then strands a trailing "3" (a non-1 partial),
    // which is NOT padding -> one extra invalid token.
    int c2[1] = { 5 }, o2[8]; nt = 0; nv = 0;
    m = compressocrat_decrypt(c2, 1, sigma, o2, 23, &nt, &nv);
    CHECK(m == 2 && nt == 2 && nv == 1 && o2[0] == g_char_to_idx['T'] && o2[1] == 23,
        "trailing-partial path wrong: m=%d nt=%d nv=%d o0=%d o1=%d", m, nt, nv, o2[0], o2[1]);
}

int main(void) {
    seed_rand(20250825u);
    init_alphabet(NULL);                                     // full 26-letter A..Z alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_kat_aca();
    test_kat_small();
    test_code_structure();
    test_reference_agreement();
    test_roundtrip();
    test_edges();
    test_invalid_paths();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

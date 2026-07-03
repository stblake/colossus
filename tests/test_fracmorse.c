//
//  Unit tests for the Fractionated Morse primitives (encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against fracmorse.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: three HAND-COMPUTED known-answer vectors under our single-'x'-separator convention
//  (keyword MORSE; SOS -> MWLR, EE -> B, ET -> C) pin the Morse table, the trigraph rank
//  (r = 9a+3b+c), the padding, and the keyed-alphabet mapping cell for cell. Then an INDEPENDENT
//  in-test reference encrypt (a fresh implementation sharing no code with fracmorse.c) is checked
//  against fracmorse_encrypt over random plaintexts / keyed alphabets, encrypt->decrypt round-trips
//  == identity over random permutations sigma and lengths (incl. the three padding residues and
//  single-letter / extreme all-E / all-T inputs), and the invalid-token path (a run of >4 dots ->
//  a filler, counted as not valid) is exercised directly.
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

// --- INDEPENDENT reference encrypt (fresh code, own Morse table) ----------------

static const char *const REF_MORSE[26] = {
    ".-",   "-...", "-.-.", "-..",  ".",    "..-.", "--.",  "....", "..",   ".---",
    "-.-",  ".-..", "--",   "-.",   "---",  ".--.", "--.-", ".-.",  "...",  "-",
    "..-",  "...-", ".--",  "-..-", "-.--", "--.."
};

static int ref_encrypt(const int *pt, int n, const int *sigma, int *out) {
    static int s[5 * 4000 + 8];
    int slen = 0;
    for (int i = 0; i < n; i++) {
        const char *m = REF_MORSE[pt[i]];
        for (int k = 0; m[k]; k++) s[slen++] = (m[k] == '-') ? 1 : 0;   // 0=dot, 1=dash
        if (i + 1 < n) s[slen++] = 2;                                   // 2 = x separator
    }
    while (slen % 3 != 0) s[slen++] = 2;
    int clen = slen / 3;
    for (int g = 0; g < clen; g++)
        out[g] = sigma[9 * s[3 * g] + 3 * s[3 * g + 1] + s[3 * g + 2]];
    return clen;
}

// --- known-answer vectors (hand-computed, single-x convention, keyword MORSE) ---

static void test_known_answer(void) {
    int sigma[26];
    build_sigma("MORSE", sigma);
    char sbuf[27]; idx_to_str(sigma, 26, sbuf);
    CHECK(strcmp(sbuf, "MORSEABCDFGHIJKLNPQTUVWXYZ") == 0,
        "MORSE keyed alphabet mismatch: got '%s'", sbuf);

    struct { const char *pt; const char *ct; } kat[] = {
        { "SOS", "MWLR" },   // ... x --- x ... , pad 1 -> ...x---x...x  ranks 0,22,15,2
        { "EE",  "B"    },   // . x .            -> .x.               rank 6
        { "ET",  "C"    },   // . x -            -> .x-               rank 7
    };
    for (int t = 0; t < 3; t++) {
        int pt[8], ct[16], back[16]; char cbuf[17];
        int n = str_to_idx(kat[t].pt, pt);
        int clen = fracmorse_encrypt(pt, n, sigma, ct);
        idx_to_str(ct, clen, cbuf);
        CHECK(strcmp(cbuf, kat[t].ct) == 0,
            "KAT %s: got '%s', want '%s'", kat[t].pt, cbuf, kat[t].ct);
        int nt = 0, nv = 0;
        int m = fracmorse_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
        CHECK(m == n && nt == n && nv == n && arrays_equal(back, pt, n),
            "KAT %s: decrypt round-trip mismatch (m=%d nt=%d nv=%d)", kat[t].pt, m, nt, nv);
    }
}

// --- fracmorse_encrypt vs the independent reference, random inputs ---------------

static void test_reference_agreement(void) {
    const char *kws[] = { "SHADOW", "KRYPTOS", "CIPHERMACHINE", "ZEBRA", "" };
    for (int t = 0; t < 2000; t++) {
        int sigma[26];
        build_sigma(kws[t % 5], sigma);
        if ((t % 5) == 4) { for (int i = 0; i < 26; i++) sigma[i] = i; shuffle(sigma, 26); }
        int n = 1 + rand_int(0, 400);
        int pt[420], a[2 * 420 + 8], b[2 * 420 + 8];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int ca = fracmorse_encrypt(pt, n, sigma, a);
        int cb = ref_encrypt(pt, n, sigma, b);
        CHECK(ca == cb && arrays_equal(a, b, ca),
            "reference disagreement (trial %d, n=%d): ca=%d cb=%d", t, n, ca, cb);
    }
}

// --- encrypt -> decrypt round-trip == identity over random sigma and lengths ----

static void test_roundtrip(void) {
    for (int t = 0; t < 5000; t++) {
        int sigma[26];
        for (int i = 0; i < 26; i++) sigma[i] = i;
        shuffle(sigma, 26);                              // any permutation, not just keyed
        int n = 1 + rand_int(0, 400);
        int pt[420], ct[2 * 420 + 8], back[2 * 420 + 8];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int clen = fracmorse_encrypt(pt, n, sigma, ct);
        int nt = 0, nv = 0;
        int m = fracmorse_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
        CHECK(m == n && nt == n && nv == n && arrays_equal(back, pt, n),
            "round-trip mismatch (trial %d, n=%d): m=%d nt=%d nv=%d", t, n, m, nt, nv);
    }
}

// --- edge cases: single letter, all-E (.), all-T (-), padding residues ----------

static void test_edges(void) {
    int sigma[26];
    for (int i = 0; i < 26; i++) sigma[i] = i;           // identity sigma
    // Single letter E ('.'): S=[.], pad xx -> [.,x,x] (rank 8) -> one cipher letter.
    int e[1] = { g_char_to_idx['E'] }, ct[8], back[8];
    int clen = fracmorse_encrypt(e, 1, sigma, ct);
    CHECK(clen == 1, "single-E cipher length %d, expected 1", clen);
    int nt = 0, nv = 0, m = fracmorse_decrypt(ct, clen, sigma, back, 23, &nt, &nv);
    CHECK(m == 1 && back[0] == e[0] && nv == 1, "single-E round-trip failed (m=%d)", m);

    // All-E and all-T (the shortest codes) at several lengths hit every padding residue.
    for (int L = 1; L <= 20; L++) {
        int pE[20], pT[20], c[64], bk[64];
        for (int i = 0; i < L; i++) { pE[i] = g_char_to_idx['E']; pT[i] = g_char_to_idx['T']; }
        int cl = fracmorse_encrypt(pE, L, sigma, c);
        int mm = fracmorse_decrypt(c, cl, sigma, bk, 23, &nt, &nv);
        CHECK(mm == L && arrays_equal(bk, pE, L), "all-E round-trip failed (L=%d)", L);
        cl = fracmorse_encrypt(pT, L, sigma, c);
        mm = fracmorse_decrypt(c, cl, sigma, bk, 23, &nt, &nv);
        CHECK(mm == L && arrays_equal(bk, pT, L), "all-T round-trip failed (L=%d)", L);
    }
}

// --- invalid-token path: a run of >4 dots is not a legal codeword -> filler ------

static void test_invalid_token(void) {
    int sigma[26];
    for (int i = 0; i < 26; i++) sigma[i] = i;           // identity: cipher letter i -> rank i
    // Cipher "AA" = ranks 0,0 -> (0,0,0)(0,0,0) = six dots, no x -> one run of length 6 (invalid).
    int cipher[2] = { 0, 0 }, out[8], nt = 0, nv = 0;
    int m = fracmorse_decrypt(cipher, 2, sigma, out, 23, &nt, &nv);
    CHECK(m == 1 && nt == 1 && nv == 0 && out[0] == 23,
        "invalid-token path wrong: m=%d nt=%d nv=%d out0=%d", m, nt, nv, out[0]);
}

int main(void) {
    seed_rand(20240703u);
    init_alphabet(NULL);                                 // full 26-letter A..Z alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_known_answer();
    test_reference_agreement();
    test_roundtrip();
    test_edges();
    test_invalid_token();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

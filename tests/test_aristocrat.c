//
//  Unit tests for the Aristocrat / Patristocrat primitives (build map / apply / invert).
//
//  Framework-free: build with `make test`, which links this against aristocrat.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy: a worked-example known-answer vector (K2, keyword KRYPTOS -> keyed ciphertext
//  alphabet KRYPTOSABCDEFGHIJLMNQUVWXZ; "THE QUICK BROWN FOX" -> "NATJQBYDRLHVGOHW") pins the
//  keyed-alphabet build and the enciphering direction. Then: the map is always a BIJECTION;
//  invert()+apply() inverts encipher for every keying mode (K1..K4); and encrypt<->decrypt
//  round-trips == identity over random maps x random plaintext (the Aristocrat and Patristocrat
//  share the same substitution math -- only the presentation differs -- so one primitive suffices).
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

static void idx_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = index_to_char(a[i]);
    out[len] = '\0';
}

// A..Z string (spaces/other stripped) -> alphabet indices.
static int letters_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        unsigned char c = (unsigned char) toupper((unsigned char) s[i]);
        if (c < 'A' || c > 'Z') continue;
        out[n++] = g_char_to_idx[c];
    }
    return n;
}

static int is_bijection(const int map[], int alpha) {
    int seen[ALPHABET_SIZE] = {0};
    for (int i = 0; i < alpha; i++) {
        if (map[i] < 0 || map[i] >= alpha || seen[map[i]]) return 0;
        seen[map[i]] = 1;
    }
    return 1;
}

int main(void) {
    init_alphabet(NULL);
    CHECK(g_alpha == 26, "alphabet size %d (want 26)", g_alpha);

    // --- KAT: K2 keyed ciphertext alphabet (keyword KRYPTOS) ---------------------------
    int cmap[26], cinv[26];
    aristocrat_build_map(ARIST_K2, "KRYPTOS", "KRYPTOS", 0, cmap);
    CHECK(is_bijection(cmap, 26), "K2 map is not a bijection");
    char alpha_str[27];
    idx_to_str(cmap, 26, alpha_str);          // K2: cmap[pt] read over pt=A..Z IS the keyed CT alphabet
    CHECK(strcmp(alpha_str, "KRYPTOSABCDEFGHIJLMNQUVWXZ") == 0,
          "K2 keyed alphabet: got %s want KRYPTOSABCDEFGHIJLMNQUVWXZ", alpha_str);

    int pt[64], ct[64], back[64];
    int L = letters_to_idx("THEQUICKBROWNFOX", pt);
    aristocrat_apply(pt, L, cmap, ct);
    char ct_str[64]; idx_to_str(ct, L, ct_str);
    CHECK(strcmp(ct_str, "NATJQBYDRLHVGOHW") == 0,
          "encrypt KAT: got %s want NATJQBYDRLHVGOHW", ct_str);

    // decrypt inverts encrypt on the KAT.
    aristocrat_invert(cmap, cinv, 26);
    aristocrat_apply(ct, L, cinv, back);
    CHECK(arrays_equal(back, pt, L), "decrypt(encrypt(KAT)) != plaintext");

    // --- every keying mode yields a bijection that inverts cleanly -----------------------
    {
        int keyings[4] = { ARIST_K1, ARIST_K2, ARIST_K3, ARIST_K4 };
        for (int ki = 0; ki < 4; ki++) {
            int m[26], mi[26];
            aristocrat_build_map(keyings[ki], "CRYPTOGRAM", "PORTABLE", 5, m);
            CHECK(is_bijection(m, 26), "keying %d not a bijection", keyings[ki]);
            aristocrat_invert(m, mi, 26);
            int p2[64], c2[64], b2[64];
            int n = letters_to_idx("PACKMYBOXWITHFIVEDOZENLIQUORJUGS", p2);
            aristocrat_apply(p2, n, m, c2);
            aristocrat_apply(c2, n, mi, b2);
            CHECK(arrays_equal(b2, p2, n), "keying %d round-trip failed", keyings[ki]);
        }
    }

    // --- round-trips: random bijection x random plaintext --------------------------------
    rng_state = 0x9e3779b9u;
    for (int t = 0; t < 4000; t++) {
        int rmap[26], rinv[26];
        for (int i = 0; i < 26; i++) rmap[i] = i;
        for (int i = 25; i > 0; i--) {
            int j = rand_int(0, i + 1);
            int tmp = rmap[i]; rmap[i] = rmap[j]; rmap[j] = tmp;
        }
        CHECK(is_bijection(rmap, 26), "random map t=%d not a bijection", t);
        aristocrat_invert(rmap, rinv, 26);
        int plen = rand_int(1, 260);
        int p[300], c[300], b[300];
        for (int i = 0; i < plen; i++) p[i] = rand_int(0, 26);
        aristocrat_apply(p, plen, rmap, c);
        aristocrat_apply(c, plen, rinv, b);
        if (!arrays_equal(b, p, plen)) { CHECK(0, "round-trip t=%d", t); break; }
    }

    printf("test_aristocrat: %d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}

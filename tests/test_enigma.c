// Framework-free unit tests for the Enigma machine primitive (src/machine/enigma.c).
// Build+run via `make test`. Pins the published KAT vectors (Ostwald-Weierud B432 and
// Gillogly's 647-letter message), the reference wiring strings, and the structural
// invariants (self-reciprocity, no fixed point, double-step, M4->M3 reduction).

#include <stdio.h>
#include <string.h>
#include "colossus.h"   // seed_rand, rand_int (inline)
#include "enigma.h"

static int failures = 0;
static int checks = 0;
#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

static int letters_to_idx(const char *s, int *out) {
    int n = 0;
    for (const char *p = s; *p; p++) {
        if (*p >= 'A' && *p <= 'Z') out[n++] = *p - 'A';
        else if (*p >= 'a' && *p <= 'z') out[n++] = *p - 'a';
    }
    return n;
}

static int arrays_equal(const int *a, const int *b, int n) {
    for (int i = 0; i < n; i++) if (a[i] != b[i]) return 0;
    return 1;
}

// Set the plugboard from a "AB CD EF" pair string.
static void set_plugs(int plug[26], const char *pairs) {
    enigma_plug_identity(plug);
    int idx[64], n = letters_to_idx(pairs, idx);
    for (int i = 0; i + 1 < n; i += 2) enigma_plug_set_pair(plug, idx[i], idx[i + 1]);
}

// ---- KAT-1: Ostwald-Weierud "Modern Breaking" p.6 -----------------------------------
// PT EINSXEINSXVIERXNULLXNULL, key B432 rit VOR, steckers AH BO CG DP FL JQ KS MU TZ WY
// -> CT FVKFCDWRIICYFHVSKQOWQTTH.
static void test_kat_b432(void) {
    EnigmaKey k;
    k.reflector = ENIGMA_UKW_B;
    k.n_wheels = 3;
    k.rotor[0] = ENIGMA_IV; k.rotor[1] = ENIGMA_III; k.rotor[2] = ENIGMA_II;  // wheels 4 3 2
    int ring[3]; letters_to_idx("RIT", ring);
    int pos[3];  letters_to_idx("VOR", pos);
    for (int i = 0; i < 3; i++) { k.ring[i] = ring[i]; k.pos[i] = pos[i]; }
    set_plugs(k.plug, "AH BO CG DP FL JQ KS MU TZ WY");

    int pt[64], ct[64], expect[64];
    int n  = letters_to_idx("EINSXEINSXVIERXNULLXNULL", pt);
    int ne = letters_to_idx("FVKFCDWRIICYFHVSKQOWQTTH", expect);
    enigma_encrypt(pt, n, &k, ct);
    CHECK(n == ne && arrays_equal(ct, expect, n), "KAT-1 (B432) ciphertext mismatch");

    // Self-reciprocal: encrypt the ciphertext with the same start key -> plaintext.
    int back[64];
    enigma_encrypt(ct, n, &k, back);
    CHECK(arrays_equal(back, pt, n), "KAT-1 (B432) decrypt did not recover plaintext");
}

// ---- KAT-2: Gillogly's 647-letter message -------------------------------------------
// Reflector B, wheels II I III, rings 1 23 4 (=A W D), msg key 2 7 9 (=B G I),
// plugboard EZ RW MV IU BL PX JO. Decrypting the ciphertext yields the German plaintext
// beginning AUFBEFEHLDESOBERSTEN... and ending ...DERKAMPFGEHTWEITERXDOENITZX.
static const char *GILLOGLY_CT =
    "QKRQWUQTZKFXZOMJFOYRHYZWVBXYSIWMMVWBLEBDMWUWBTVHMRFLKSDCCEXIYPAHRMPZIOVBBRV"
    "LNHZUPOSYEIPWJTUGYOSLAOXRHKVCHQOSVDTRBPDJEUKSBBXHTTGVHGFICACVGUVOQFAQWBKXZJ"
    "SQJFZPEVJROJTOESLBQHQTRAAHXVYAUHTNBGIBVCLBLXCYBDMQRTVPYKFFZXNDDPCCJBHQFDKXE"
    "EYWPBYQWDXDRDHNIGDXEUJJPVMHUKPCFHLLFERAZHZOHXDGBKOQXKTLDVDCWKAEDHCPHJIWZMMT"
    "UAMQENNFCHUIAWCCHNCFYPWUARBBNIEPHGDDKMDQLMSNMTWOHMAUHRHGCUMQPKQRKDVSWVMTYVN"
    "FFDDSKIISXONXQHHLIYQSDFHENCMCOMREZQDRPBMRVPQTVRSWZPGLPITRVIBPXXHPRFISZTPUEP"
    "LKOTTXNAZMHTJPCHAASFZLEFCEZUTPYBAOSKPZCJCYZOVAPZZVELBLLZEVDCHRMIOYEPFVUGNDL"
    "ENISXYCHKSJUWVXUSBITDEQTCNKRLSNXMXYZGCUPAWFULTZZSFAHMPXGLLNZRXYJNSKYNQAMZBU"
    "GFZJCURWGTQZCTLLOIEKAOISKHAAQFOPFUZIRTLWEVYWMDN";

static void test_kat_gillogly(void) {
    EnigmaKey k;
    k.reflector = ENIGMA_UKW_B;
    k.n_wheels = 3;
    k.rotor[0] = ENIGMA_II; k.rotor[1] = ENIGMA_I; k.rotor[2] = ENIGMA_III;  // II I III
    int ring[3]; letters_to_idx("AWD", ring);   // rings 1 23 4 (1-based) = A W D (0-based)
    int pos[3];  letters_to_idx("BGI", pos);     // msg  2 7 9 (1-based) = B G I (0-based)
    for (int i = 0; i < 3; i++) { k.ring[i] = ring[i]; k.pos[i] = pos[i]; }
    set_plugs(k.plug, "EZ RW MV IU BL PX JO");

    int ct[700], pt[700];
    int n = letters_to_idx(GILLOGLY_CT, ct);
    CHECK(n == 647, "KAT-2 (Gillogly) ciphertext length %d, expected 647", n);
    enigma_encrypt(ct, n, &k, pt);   // decrypt == encrypt

    int pref[64]; int lp = letters_to_idx("AUFBEFEHLDESOBERSTENBEFEHLSHABERS", pref);
    CHECK(arrays_equal(pt, pref, lp), "KAT-2 (Gillogly) plaintext prefix mismatch");
    int suf[80]; int ls = letters_to_idx("DERFUEHRERISTTOTXDERKAMPFGEHTWEITERXDOENITZX", suf);
    CHECK(arrays_equal(pt + (n - ls), suf, ls), "KAT-2 (Gillogly) plaintext suffix mismatch");
}

// ---- KAT-3: reference wiring strings ------------------------------------------------
static void test_wiring_strings(void) {
    struct { int id; const char *w; const char *notch; } R[] = {
        { ENIGMA_I,    "EKMFLGDQVZNTOWYHXUSPAIBRCJ", "Q"  },
        { ENIGMA_II,   "AJDKSIRUXBLHWTMCQGZNPYFVOE", "E"  },
        { ENIGMA_III,  "BDFHJLCPRTXVZNYEIWGAKMUSQO", "V"  },
        { ENIGMA_IV,   "ESOVPZJAYQUIRHXLNFTGKDCMWB", "J"  },
        { ENIGMA_V,    "VZBRGITYUPSDNHLXAWMJQOFECK", "Z"  },
        { ENIGMA_VI,   "JPGVOUMFYQBENHZRDKASXLICTW", "ZM" },
        { ENIGMA_VII,  "NZJHGRCXMYSWBOUFAIVLPEKQDT", "ZM" },
        { ENIGMA_VIII, "FKQHTLXOCBJSPDZRAMEWNIUYGV", "ZM" },
        { ENIGMA_BETA, "LEYJVCNIXWPBQMDRTAKZGFUHOS", ""   },
        { ENIGMA_GAMMA,"FSOKANUERHMBTIYCWLQPZXVGJD", ""   },
    };
    for (int i = 0; i < 10; i++) {
        CHECK(strcmp(enigma_rotor_wiring(R[i].id), R[i].w) == 0, "rotor %s wiring", enigma_rotor_name(R[i].id));
        CHECK(enigma_num_notches(R[i].id) == (int) strlen(R[i].notch), "rotor %s notch count", enigma_rotor_name(R[i].id));
    }
    CHECK(strcmp(enigma_reflector_wiring(ENIGMA_UKW_B), "YRUHQSLDPXNGOKMIEBFZCWVJAT") == 0, "UKW-B");
    CHECK(strcmp(enigma_reflector_wiring(ENIGMA_UKW_C), "FVPJIAOYEDRZXWGCTKUQSBNMHL") == 0, "UKW-C");
    CHECK(strcmp(enigma_reflector_wiring(ENIGMA_UKW_B_THIN), "ENKQAUYWJICOPBLMDXZVFTHRGS") == 0, "UKW-B thin");
    CHECK(strcmp(enigma_reflector_wiring(ENIGMA_UKW_C_THIN), "RDOBJNTKVEHMLFCWZAXGYIPSUQ") == 0, "UKW-C thin");

    // Reflectors are fixed-point-free involutions.
    for (int f = 0; f < ENIGMA_N_REFLECTORS; f++) {
        const char *w = enigma_reflector_wiring(f);
        for (int i = 0; i < 26; i++) {
            int j = w[i] - 'A';
            CHECK(j != i, "reflector %s has fixed point at %c", enigma_reflector_name(f), 'A' + i);
            CHECK(w[j] - 'A' == i, "reflector %s not involutory at %c", enigma_reflector_name(f), 'A' + i);
        }
    }
}

// ---- Double-step: rotors I II III, rings AAA, start (A,D,U). -------------------------
// step1 -> (A,D,V); step2 -> (A,E,W) [right notch V carries middle]; step3 -> (B,F,X)
// [middle now on notch E: double-steps itself AND carries the left].
static void test_double_step(void) {
    EnigmaKey k;
    k.reflector = ENIGMA_UKW_B; k.n_wheels = 3;
    k.rotor[0] = ENIGMA_I; k.rotor[1] = ENIGMA_II; k.rotor[2] = ENIGMA_III;
    for (int i = 0; i < 3; i++) k.ring[i] = 0;
    k.pos[0] = 0; k.pos[1] = 3; k.pos[2] = 20;   // A D U
    int exp[3][3] = { {0,3,21}, {0,4,22}, {1,5,23} };  // A D V / A E W / B F X
    for (int s = 0; s < 3; s++) {
        enigma_step(&k);
        CHECK(k.pos[0] == exp[s][0] && k.pos[1] == exp[s][1] && k.pos[2] == exp[s][2],
              "double-step %d: got (%d,%d,%d)", s + 1, k.pos[0], k.pos[1], k.pos[2]);
    }
}

static void random_key(EnigmaKey *k, int m4) {
    if (m4) {
        k->reflector = ENIGMA_UKW_B_THIN + rand_int(0, 2);   // thin B or C
        k->n_wheels = 4;
        k->rotor[0] = ENIGMA_BETA + rand_int(0, 2);          // Greek
    } else {
        k->reflector = rand_int(0, 2);                       // wide B or C
        k->n_wheels = 3;
    }
    int base = k->n_wheels - 3;                              // first stepping wheel index
    int pool[8] = {0,1,2,3,4,5,6,7};
    for (int i = 0; i < 3; i++) { int j = rand_int(i, 8); int t = pool[i]; pool[i] = pool[j]; pool[j] = t; }
    for (int i = 0; i < 3; i++) k->rotor[base + i] = pool[i];
    for (int i = 0; i < k->n_wheels; i++) { k->ring[i] = rand_int(0, 26); k->pos[i] = rand_int(0, 26); }
    enigma_plug_identity(k->plug);
    int avail[26]; for (int i = 0; i < 26; i++) avail[i] = i; int na = 26;
    int nplugs = rand_int(0, 11);
    for (int p = 0; p < nplugs && na >= 2; p++) {
        int i = rand_int(0, na); int a = avail[i]; avail[i] = avail[--na];
        int j = rand_int(0, na); int b = avail[j]; avail[j] = avail[--na];
        enigma_plug_set_pair(k->plug, a, b);
    }
}

// Self-reciprocity + no-fixed-point over random keys x lengths (M3 and M4).
static void test_roundtrip(void) {
    int recip_fail = 0, fixpt_fail = 0;
    for (int t = 0; t < 4000; t++) {
        EnigmaKey k; random_key(&k, t & 1);
        int n = rand_int(1, 300);
        int pt[300], ct[300], back[300];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        enigma_encrypt(pt, n, &k, ct);
        enigma_encrypt(ct, n, &k, back);   // same start key
        if (!arrays_equal(back, pt, n)) recip_fail++;
        for (int i = 0; i < n; i++) if (ct[i] == pt[i]) { fixpt_fail++; break; }
    }
    CHECK(recip_fail == 0, "self-reciprocity failed in %d/4000 trials", recip_fail);
    CHECK(fixpt_fail == 0, "fixed point produced in %d/4000 trials", fixpt_fail);
}

// M4 with Greek Beta at position A / ring A and the thin-B reflector reduces EXACTLY to
// M3 with UKW-B (the documented backward-compatibility mode).
static void test_m4_reduces_to_m3(void) {
    int fails = 0;
    for (int t = 0; t < 500; t++) {
        EnigmaKey m3;
        m3.reflector = ENIGMA_UKW_B; m3.n_wheels = 3;
        int pool[8] = {0,1,2,3,4,5,6,7};
        for (int i = 0; i < 3; i++) { int j = rand_int(i, 8); int tt = pool[i]; pool[i] = pool[j]; pool[j] = tt; }
        for (int i = 0; i < 3; i++) { m3.rotor[i] = pool[i]; m3.ring[i] = rand_int(0, 26); m3.pos[i] = rand_int(0, 26); }
        enigma_plug_identity(m3.plug);

        EnigmaKey m4 = m3;
        m4.reflector = ENIGMA_UKW_B_THIN; m4.n_wheels = 4;
        m4.rotor[0] = ENIGMA_BETA; m4.ring[0] = 0; m4.pos[0] = 0;   // Greek @ A, ring A
        for (int i = 0; i < 3; i++) { m4.rotor[i + 1] = m3.rotor[i]; m4.ring[i + 1] = m3.ring[i]; m4.pos[i + 1] = m3.pos[i]; }

        int n = rand_int(1, 200);
        int pt[200], c3[200], c4[200];
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        enigma_encrypt(pt, n, &m3, c3);
        enigma_encrypt(pt, n, &m4, c4);
        if (!arrays_equal(c3, c4, n)) fails++;
    }
    CHECK(fails == 0, "M4(Beta@A,B-thin) != M3(UKW-B) in %d/500 trials", fails);
}

int main(void) {
    seed_rand(20260907u);
    enigma_init();
    test_kat_b432();
    test_kat_gillogly();
    test_wiring_strings();
    test_double_step();
    test_roundtrip();
    test_m4_reduces_to_m3();
    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

//
//  Unit tests for the Chaocipher primitive (encrypt / decrypt).
//
//  Framework-free: build with `make test`, which links this against chaocipher.c + utils.c.
//  Exits non-zero if any check fails.
//
//  Strategy:
//   * KAT: Moshe Rubin's "Chaocipher Revealed: The Algorithm" (2010) worked example --
//     starting alphabets HXUCZVAMDSLKPEFJRIGTWOBNYQ / PTLNBQDEOYSFAVZKGJRIHWXUMC encipher
//     WELLDONEISBETTERTHANWELLSAID -> OAHQHCNYNXTSZJRRHJBYHQKSOUJY (both directions), plus the
//     paper's two structural checks (after each step the next left[0] == the emitted ct letter
//     and the next right[25] == the enciphered pt letter) and the fact that a letter CAN
//     encipher to itself (N->N in the KAT).
//   * An INDEPENDENT reference encrypt (a fresh implementation building each permuted alphabet
//     into a new buffer, sharing no code with chaocipher.c) checked byte-identical over random
//     starting alphabets x lengths -- the strongest guard for the dual-disk stepping offsets.
//   * encrypt->decrypt round-trips == identity over random alphabets x lengths (incl. 1 and >500).
//   * Invariants: both disks stay permutations of 0..25 after every step; and the rotation
//     equivalence (rotating BOTH starting alphabets by the same k gives identical ciphertext).
//   * Edge cases: length 1, repeated letters, all-same letter.
//

#include "colossus.h"
#include "chaocipher.h"

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

// A..Z string -> 0..25 indices (Chaocipher is always the full 26-letter alphabet).
static int str_to_idx(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        int c = toupper((unsigned char) s[i]);
        if (c >= 'A' && c <= 'Z') out[n++] = c - 'A';
    }
    return n;
}

static void idx_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = (char) ('A' + a[i]);
    out[len] = '\0';
}

static int is_permutation(const int a[26]) {
    int seen[26] = {0};
    for (int i = 0; i < 26; i++) {
        if (a[i] < 0 || a[i] >= 26 || seen[a[i]]) return 0;
        seen[a[i]] = 1;
    }
    return 1;
}

// ---------------------------------------------------------------------------
//  Independent reference implementation (shares NO code with chaocipher.c).
//  Each permutation builds a fresh alphabet buffer element-by-element, following the paper's
//  1-indexed description literally, rather than the primitive's in-place rotate + block shift.
// ---------------------------------------------------------------------------

// Rotate so that value `v` sits at index 0, writing into out[] (out != a).
static void ref_bring_to_zenith(const int a[26], int v, int out[26]) {
    int p = 0; while (a[p] != v) p++;
    for (int i = 0; i < 26; i++) out[i] = a[(p + i) % 26];
}

// Left (ciphertext) disk permutation, precondition: the used ct letter is at index 0.
// New = [ a0 | a2 a3 ... a13 | a1 | a14 ... a25 ]  (extract index 1, reinsert at index 13).
static void ref_permute_left(int a[26]) {
    int out[26];
    out[0] = a[0];
    for (int i = 1; i <= 12; i++) out[i] = a[i + 1];
    out[13] = a[1];
    for (int i = 14; i < 26; i++) out[i] = a[i];
    for (int i = 0; i < 26; i++) a[i] = out[i];
}

// Right (plaintext) disk permutation, precondition: the used pt letter is at index 0.
// First shift one more (successor to index 0), then extract index 2, reinsert at index 13.
static void ref_permute_right(int a[26]) {
    int r[26];
    for (int i = 0; i < 26; i++) r[i] = a[(i + 1) % 26];   // one extra shift left
    int out[26];
    out[0] = r[0]; out[1] = r[1];
    for (int i = 2; i <= 12; i++) out[i] = r[i + 1];
    out[13] = r[2];
    for (int i = 14; i < 26; i++) out[i] = r[i];
    for (int i = 0; i < 26; i++) a[i] = out[i];
}

// Reference encrypt; if trace_left0/trace_right25 non-NULL, records the NEXT state's left[0]
// and right[25] after each step (the paper's structural checks).
static void ref_encrypt(const int plain[], int n, const int left0[26], const int right0[26],
                        int out[], int trace_left0[], int trace_right25[]) {
    int left[26], right[26];
    for (int k = 0; k < 26; k++) { left[k] = left0[k]; right[k] = right0[k]; }
    for (int s = 0; s < n; s++) {
        int p = plain[s];
        int i = 0; while (right[i] != p) i++;
        out[s] = left[i];
        int zl[26], zr[26];
        ref_bring_to_zenith(left, left[i], zl);   // ct letter to zenith
        ref_bring_to_zenith(right, p, zr);        // pt letter to zenith
        for (int k = 0; k < 26; k++) { left[k] = zl[k]; right[k] = zr[k]; }
        ref_permute_left(left);
        ref_permute_right(right);
        if (trace_left0)   trace_left0[s]   = left[0];
        if (trace_right25) trace_right25[s] = right[25];
    }
}

// ---------------------------------------------------------------------------

static void rand_perm(int a[26]) {
    for (int i = 0; i < 26; i++) a[i] = i;
    shuffle(a, 26);
}

// ---- KAT ----

static void test_known_answer(void) {
    int left[26], right[26], pt[64], ct[64], back[64];
    int n = str_to_idx("WELLDONEISBETTERTHANWELLSAID", pt);
    str_to_idx("HXUCZVAMDSLKPEFJRIGTWOBNYQ", left);
    str_to_idx("PTLNBQDEOYSFAVZKGJRIHWXUMC", right);

    int expect[64];
    str_to_idx("OAHQHCNYNXTSZJRRHJBYHQKSOUJY", expect);

    chaocipher_encrypt(pt, n, left, right, ct);
    char got[64]; idx_to_str(ct, n, got);
    CHECK(arrays_equal(ct, expect, n), "KAT encrypt mismatch: got %s", got);

    chaocipher_decrypt(ct, n, left, right, back);
    CHECK(arrays_equal(back, pt, n), "KAT decrypt did not invert encrypt");

    // A letter CAN encipher to itself: pt[6]=N enciphers to ct[6]=N in this KAT.
    CHECK(pt[6] == ('N' - 'A') && ct[6] == ('N' - 'A'),
          "expected the self-enciphering N->N at position 6 (pt=%d ct=%d)", pt[6], ct[6]);

    // Paper structural checks: next left[0] == emitted ct, next right[25] == enciphered pt.
    int tl0[64], tr25[64], refct[64];
    ref_encrypt(pt, n, left, right, refct, tl0, tr25);
    CHECK(arrays_equal(refct, expect, n), "reference encrypt disagrees with the KAT");
    int ok_left0 = 1, ok_right25 = 1;
    for (int s = 0; s < n - 1; s++) {
        if (tl0[s]   != ct[s]) ok_left0 = 0;      // next state's left[0]  == this step's ct
        if (tr25[s]  != pt[s]) ok_right25 = 0;    // next state's right[25]== this step's pt
    }
    CHECK(ok_left0,  "structural check failed: next left[0] should equal the emitted ct");
    CHECK(ok_right25, "structural check failed: next right[25] should equal the enciphered pt");
}

// ---- independent-reference agreement over random keys x lengths ----

static void test_reference_agreement(void) {
    int left[26], right[26], pt[600], a[600], b[600];
    int bad = 0;
    for (int trial = 0; trial < 4000 && bad < 5; trial++) {
        rand_perm(left); rand_perm(right);
        int n = 1 + rand_int(0, 500);
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        chaocipher_encrypt(pt, n, left, right, a);
        ref_encrypt(pt, n, left, right, b, NULL, NULL);
        if (!arrays_equal(a, b, n)) {
            bad++;
            CHECK(0, "reference disagreement at trial %d (n=%d)", trial, n);
        }
    }
    CHECK(bad == 0, "primitive disagreed with the independent reference on %d trials", bad);
}

// ---- round-trip ----

static void test_roundtrip(void) {
    int left[26], right[26], pt[1200], ct[1200], back[1200];
    int bad = 0;
    for (int trial = 0; trial < 4000 && bad < 5; trial++) {
        rand_perm(left); rand_perm(right);
        int n = 1 + rand_int(0, 1000);
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        chaocipher_encrypt(pt, n, left, right, ct);
        chaocipher_decrypt(ct, n, left, right, back);
        if (!arrays_equal(back, pt, n)) { bad++; CHECK(0, "round-trip failed (trial %d n=%d)", trial, n); }
    }
    CHECK(bad == 0, "round-trip failed on %d trials", bad);
}

// ---- invariants: disks stay permutations; rotation equivalence ----

static void test_invariants(void) {
    // Disks remain permutations of 0..25 after every step. Re-run the primitive stepping via
    // the exposed disk steppers on independent copies and check permutation-ness each step.
    int left[26], right[26], pt[300];
    for (int trial = 0; trial < 200; trial++) {
        rand_perm(left); rand_perm(right);
        int n = 1 + rand_int(0, 250);
        for (int i = 0; i < n; i++) pt[i] = rand_int(0, 26);
        int L[26], R[26];
        for (int k = 0; k < 26; k++) { L[k] = left[k]; R[k] = right[k]; }
        int ok = 1;
        for (int s = 0; s < n; s++) {
            int p = pt[s], i = 0; while (R[i] != p) i++;
            chaocipher_step_left(L, i);
            chaocipher_step_right(R, i);
            if (!is_permutation(L) || !is_permutation(R)) { ok = 0; break; }
        }
        CHECK(ok, "a disk stopped being a permutation during stepping (trial %d)", trial);
    }

    // Rotation equivalence: rotating BOTH starting alphabets by the same k gives identical ct.
    int rl[26], rr[26], base_ct[600], rot_ct[600];
    for (int trial = 0; trial < 300; trial++) {
        rand_perm(left); rand_perm(right);
        int n = 1 + rand_int(0, 400);
        int ptn[600]; for (int i = 0; i < n; i++) ptn[i] = rand_int(0, 26);
        chaocipher_encrypt(ptn, n, left, right, base_ct);
        int k = rand_int(1, 26);
        for (int i = 0; i < 26; i++) { rl[i] = left[(i + k) % 26]; rr[i] = right[(i + k) % 26]; }
        chaocipher_encrypt(ptn, n, rl, rr, rot_ct);
        CHECK(arrays_equal(base_ct, rot_ct, n),
              "rotation equivalence broken (trial %d k=%d n=%d)", trial, k, n);
    }
}

// ---- edges ----

static void test_edges(void) {
    int left[26], right[26], ct[600], back[600];
    rand_perm(left); rand_perm(right);

    int one = 3;                                   // single letter
    chaocipher_encrypt(&one, 1, left, right, ct);
    chaocipher_decrypt(ct, 1, left, right, back);
    CHECK(back[0] == one, "single-letter round-trip failed");

    int same[400];                                 // all the same letter (heavy repetition)
    for (int i = 0; i < 400; i++) same[i] = 7;
    chaocipher_encrypt(same, 400, left, right, ct);
    chaocipher_decrypt(ct, 400, left, right, back);
    CHECK(arrays_equal(back, same, 400), "all-same-letter round-trip failed");

    int rep[400];                                  // a short cycle repeated
    const char *pat = "THETHETHE";
    int m = 0; for (int i = 0; i < 400; i++) { rep[i] = pat[m] - 'A'; m = (m + 1) % 9; }
    chaocipher_encrypt(rep, 400, left, right, ct);
    chaocipher_decrypt(ct, 400, left, right, back);
    CHECK(arrays_equal(back, rep, 400), "repeated-pattern round-trip failed");
}

int main(void) {
    seed_rand(20260908u);
    init_alphabet(NULL);                           // full 26-letter A..Z alphabet
    CHECK(g_alpha == 26, "alphabet size %d, expected 26", g_alpha);

    test_known_answer();
    test_reference_agreement();
    test_roundtrip();
    test_invariants();
    test_edges();

    printf("\n%d checks, %d failures\n", checks, failures);
    if (failures) { printf("TESTS FAILED\n"); return 1; }
    printf("ALL TESTS PASSED\n");
    return 0;
}

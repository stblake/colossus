//
//  Unit tests for the Period column order transposition primitive
//  (period_column_transform), AZdecrypt's "Period column order".
//
//  Framework-free: build with `make test`, which links this against period_column.c
//  + utils.c. Exits non-zero if any check fails.
//
//  Strategy:
//   1. A tiny HAND-COMPUTED known-answer vector pinning the exact column-visit /
//      fill / readout convention (len 6, dx 3, period 2: [0,1,2,3,4,5] -> [0,2,1,3,5,4]),
//      plus its utp==1 inverse.
//   2. The REAL AZdecrypt worked example as a two-stage composition KAT: the 168-char
//      spaced ciphertext decrypts, via [4x42,TP,P:3] then [56x3,UTP,P:2], to
//      "I LIKE KILLING PEOPLE ..." exactly -- and the inverse composition re-encrypts
//      the plaintext back to the ciphertext. This anchors both the single-stage
//      convention AND the multi-stage pipeline against an external solver, spaces
//      included (spaces are carried as genuine grid cells).
//   3. Exact round-trip (utp==0 then utp==1 == identity) over random COMPLETE grids
//      (dx | len) x random periods, including period 1, dx == len (single row), and
//      values that carry negative "space" sentinels.
//   4. Multiset preservation + exact length over random INCOMPLETE grids (dx does not
//      divide len): the transform is always a length-preserving bijection of the value
//      multiset (the solver restricts itself to complete grids, so exact per-stage
//      inversion is only required, and only asserted, there).
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

// Opaque char <-> int mapping (the transform only permutes cells, so any injective
// encoding exercises the exact same permutation; we use the raw byte value).
static int str_to_ints(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) out[n++] = (unsigned char) s[i];
    return n;
}
static void ints_to_str(const int a[], int len, char out[]) {
    for (int i = 0; i < len; i++) out[i] = (char) a[i];
    out[len] = '\0';
}

// Is `out` a permutation of `in` (same length, same multiset)?
static int is_permutation(const int in[], const int out[], int len) {
    // O(len^2) is fine for the test sizes; avoids assuming a value range.
    static char used[MAX_CIPHER_LENGTH];
    for (int i = 0; i < len; i++) used[i] = 0;
    for (int i = 0; i < len; i++) {
        int found = 0;
        for (int j = 0; j < len; j++) {
            if (!used[j] && out[i] == in[j]) { used[j] = 1; found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

int main(void) {
    init_alphabet(NULL);   // not strictly needed (opaque values), but keeps parity
    seed_rand(0xC01055u);

    static int a[MAX_CIPHER_LENGTH], b[MAX_CIPHER_LENGTH], c[MAX_CIPHER_LENGTH];
    static int ref[MAX_CIPHER_LENGTH];

    // --- 1. tiny hand-computed KAT -------------------------------------------------
    // len 6, dx 3 (dy 2), period 2, utp 0. Column visit order (period 2 over 3 cols) is
    // 0,2,1; grid rows [0,1,2]/[3,4,5] -> columns permuted to [0,2,1] -> readout
    // 0,2,1,3,5,4.
    {
        int in6[6] = {0, 1, 2, 3, 4, 5};
        int want[6] = {0, 2, 1, 3, 5, 4};
        period_column_transform(in6, a, 6, 3, 2, 0);
        CHECK(arrays_equal(a, want, 6), "hand KAT dx3 p2 utp0");
        period_column_transform(a, b, 6, 3, 2, 1);      // inverse
        CHECK(arrays_equal(b, in6, 6), "hand KAT dx3 p2 utp1 inverts");
    }

    // --- 2. AZdecrypt worked-example composition KAT (spaces are grid cells) --------
    // Ciphertext C decrypts, applying [4x42,TP,P:3] then [56x3,UTP,P:2], to the
    // plaintext PT; and the inverse composition re-encrypts PT back to C.
    {
        const char *C =
            "I E LTIIKI E  SKIOSLLM INCUG  HPEUFOP.NLEI  B TECSIAUM SOEMREI  F NUNHT T "
            "EHAOFN ERKITSLLB INCEG UAWIESLDM  GNAA  LISFO TA HELL M .OSOTT K DALING "
            "LEROSOUEMS HTANNIIM GA";
        const char *PT =
            "I LIKE KILLING PEOPLE BECAUSE IT IS SO MUCH FUN. IT IS MORE FUN THAN "
            "KILLING WILD GAME IN THE FOREST BECAUSE MAN IS THE MOST DANGEROUS ANIMAL "
            "OF ALL. TO KILL SOMETHING ";
        int nC = str_to_ints(C, a);
        int nPT = str_to_ints(PT, ref);
        CHECK(nC == 168 && nPT == 168, "AZ KAT lengths (%d, %d)", nC, nPT);
        // decrypt: F(4,3,TP) then F(56,2,UTP)
        period_column_transform(a, b, nC, 4, 3, 0);
        period_column_transform(b, c, nC, 56, 2, 1);
        char got[200]; ints_to_str(c, nC, got);
        CHECK(arrays_equal(c, ref, nC), "AZ KAT decrypt -> plaintext\n  got: %s", got);
        // re-encrypt: inverse in reverse -> F(56,2,TP) then F(4,3,UTP)
        period_column_transform(ref, b, nPT, 56, 2, 0);
        period_column_transform(b, c, nPT, 4, 3, 1);
        CHECK(arrays_equal(c, a, nPT), "AZ KAT re-encrypt plaintext -> ciphertext");
    }

    // --- 3. exact round-trip over random COMPLETE grids ----------------------------
    // utp==0 then utp==1 must be the identity when dx | len.
    for (int t = 0; t < 20000; t++) {
        int len = 4 + rand_int(0, 400);
        // pick a divisor dx of len in [1, len]
        int dx;
        do { dx = 1 + rand_int(0, len); } while (len % dx != 0);
        int period = (dx >= 2) ? (1 + rand_int(0, dx)) : 1;   // 1 .. dx
        // values include negative "space" sentinels to prove they ride through
        for (int i = 0; i < len; i++)
            a[i] = (rand_int(0, 10) == 0) ? -(33 + rand_int(0, 20)) : rand_int(0, 26);
        period_column_transform(a, b, len, dx, period, 0);
        period_column_transform(b, c, len, dx, period, 1);
        if (!arrays_equal(c, a, len)) {
            CHECK(0, "complete-grid round-trip len=%d dx=%d p=%d", len, dx, period);
            break;
        }
        // and the forward transform is a permutation of the input
        if (!is_permutation(a, b, len)) {
            CHECK(0, "complete-grid permutation len=%d dx=%d p=%d", len, dx, period);
            break;
        }
    }
    CHECK(failures == 0 || checks > 0, "complete-grid stress ran");

    // period 1 is the identity column order (complete grid).
    {
        int len = 60;
        for (int i = 0; i < len; i++) a[i] = rand_int(0, 26);
        period_column_transform(a, b, len, 12, 1, 0);
        CHECK(arrays_equal(a, b, len), "period 1 is identity");
    }

    // dx == len: single row, columns are the whole message.
    {
        int len = 30;
        for (int i = 0; i < len; i++) a[i] = rand_int(0, 26);
        period_column_transform(a, b, len, len, 7, 0);
        period_column_transform(b, c, len, len, 7, 1);
        CHECK(arrays_equal(a, c, len), "dx==len (single row) round-trips");
        CHECK(is_permutation(a, b, len), "dx==len is a permutation");
    }

    // --- 4. INCOMPLETE grids: length + multiset preserved (bijection) --------------
    int incomplete_roundtrips = 0, incomplete_total = 0;
    for (int t = 0; t < 20000; t++) {
        int len = 5 + rand_int(0, 400);
        int dx = 3 + rand_int(0, len - 3 > 0 ? len - 3 : 1);
        if (len % dx == 0) continue;                  // want an incomplete grid
        int period = 1 + rand_int(0, dx);
        for (int i = 0; i < len; i++) a[i] = rand_int(0, 26);
        period_column_transform(a, b, len, dx, period, 0);
        // Always: length preserved and a permutation of the input multiset.
        if (!is_permutation(a, b, len)) {
            CHECK(0, "incomplete-grid permutation len=%d dx=%d p=%d", len, dx, period);
            break;
        }
        // Characterize (not assert) how often utp==1 still inverts on incomplete grids.
        period_column_transform(b, c, len, dx, period, 1);
        incomplete_total++;
        if (arrays_equal(a, c, len)) incomplete_roundtrips++;
    }
    CHECK(failures == 0 || checks > 0, "incomplete-grid stress ran");
    printf("  [incomplete-grid utp inverse holds for %d/%d random cases]\n",
           incomplete_roundtrips, incomplete_total);

    printf("\n%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}

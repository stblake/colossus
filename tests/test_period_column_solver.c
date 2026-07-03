//
//  In-process capability tests for the Period column order solver
//  (period_column_search).
//
//  Framework-free: build with `make testopt`. The solver is DETERMINISTIC and
//  EXHAUSTIVE (it tries every complete-grid stage at depth 1 and every ordered pair
//  at depth 2, scoring each decrypt with the shared n-gram fitness and keeping the
//  global best), so -- unlike a stochastic climber -- recovery is guaranteed inside
//  the searched (width, period, depth<=2) space whenever the true plaintext is the
//  highest-scoring reachable arrangement. It is a pure transposition, so it rides the
//  reward-only quadgram table (no -logprob), like the rest of the transposition family.
//
//  Checks:
//    1. depth-1 exact recovery across widths / periods / directions;
//    2. depth-2 exact recovery (the AZdecrypt two-stage scenario);
//    3. the REAL AZdecrypt worked example end to end: the 168-char spaced ciphertext
//       -> "I LIKE KILLING PEOPLE ..." with the recovered stages [4x42,TP,P:3]
//       [56x3,UTP,P:2] -- an asserted KAT reproducing the headline result;
//    4. a length cliff (recovery vs length, depth 2) -- printed + a long-length floor.
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "scoring.h"              // load_ngrams
#include "period_column_solver.h" // period_column_search / period_column_transform

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
#define NGRAM_SIZE 4

static float *ngrams = NULL;

static int arrays_equal(const int a[], const int b[], int len) {
    for (int i = 0; i < len; i++) if (a[i] != b[i]) return 0;
    return 1;
}

// A long chunk of natural English (Pride and Prejudice, opening), letters only.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEINWANTOFAWIFE"
    "HOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGANEIGHBOURHOODTHISTRUTHIS"
    "SOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOR"
    "OTHEROFTHEIRDAUGHTERSMYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLET"
    "ATLASTMRBENNETREPLIEDTHATHEHADNOTBUTITISRETURNEDSHEFORMRSLONGHASJUSTBEENHEREANDSHETOLDMEALLABOUT";

// Encode a bare A..Z string to alphabet indices (0..25).
static int encode_letters(const char *s, int len, int out[]) {
    int n = 0;
    for (int i = 0; i < len && s[i]; i++) out[n++] = g_char_to_idx[toupper((unsigned char) s[i])];
    return n;
}

// Encode a spaced string the way the solver's decode_cipher/ord does: letters -> 0..25,
// anything else -> a reversible negative sentinel.
static int encode_spaced(const char *s, int out[]) {
    int n = 0;
    for (int i = 0; s[i]; i++) {
        unsigned char c = (unsigned char) s[i];
        int v = g_char_to_idx[toupper(c)];
        out[n++] = (v >= 0) ? v : (-(int) c - 1);
    }
    return n;
}

int main(void) {
    init_alphabet(NULL);                    // full 26-letter alphabet (no J->I)
    ngrams = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!ngrams) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    static int pt[MAX_CIPHER_LENGTH], cipher[MAX_CIPHER_LENGTH], got[MAX_CIPHER_LENGTH];
    PCStage stages[2];
    int nstg = 0;

    // --- 1. depth-1 exact recovery --------------------------------------------------
    // For each case: encrypt the plaintext with one forward stage, then the exhaustive
    // depth-1 search must recover it exactly (the inverse is a single forward stage).
    struct { int len, dx, p, utp; } d1[] = {
        {168, 4, 3, 0}, {168, 56, 2, 1}, {168, 24, 5, 0}, {168, 42, 7, 1},
        {180, 12, 5, 0}, {180, 60, 3, 1}, {200, 8, 3, 1}, {150, 30, 4, 0},
    };
    int d1_ok = 0, d1_n = (int) (sizeof d1 / sizeof d1[0]);
    for (int t = 0; t < d1_n; t++) {
        int len = encode_letters(PLAINTEXT, d1[t].len, pt);
        period_column_transform(pt, cipher, len, d1[t].dx, d1[t].p, d1[t].utp);
        period_column_search(cipher, len, 1, ngrams, NGRAM_SIZE,
                             NULL, NULL, 0, 1.0f, 0.0f, got, stages, &nstg);
        if (arrays_equal(got, pt, len)) d1_ok++;
        else CHECK(0, "depth-1 recovery len=%d dx=%d p=%d utp=%d", d1[t].len, d1[t].dx, d1[t].p, d1[t].utp);
    }
    CHECK(d1_ok == d1_n, "depth-1 exact recovery %d/%d", d1_ok, d1_n);

    // --- 2. depth-2 exact recovery (AZdecrypt two-stage scenario) -------------------
    struct { int len, dx1, p1, u1, dx2, p2, u2; } d2[] = {
        {168, 56, 2, 0, 4, 3, 1},     // inverse of AZ's ([4,3,TP],[56,2,UTP]) decrypt
        {168, 42, 4, 0, 8, 3, 1},
        {180, 60, 3, 0, 12, 5, 1},
        {240, 16, 5, 0, 15, 4, 1},
    };
    int d2_ok = 0, d2_n = (int) (sizeof d2 / sizeof d2[0]);
    for (int t = 0; t < d2_n; t++) {
        int len = encode_letters(PLAINTEXT, d2[t].len, pt);
        static int mid[MAX_CIPHER_LENGTH];
        period_column_transform(pt, mid, len, d2[t].dx1, d2[t].p1, d2[t].u1);
        period_column_transform(mid, cipher, len, d2[t].dx2, d2[t].p2, d2[t].u2);
        period_column_search(cipher, len, 2, ngrams, NGRAM_SIZE,
                             NULL, NULL, 0, 1.0f, 0.0f, got, stages, &nstg);
        if (arrays_equal(got, pt, len)) d2_ok++;
        else CHECK(0, "depth-2 recovery len=%d", d2[t].len);
    }
    CHECK(d2_ok == d2_n, "depth-2 exact recovery %d/%d", d2_ok, d2_n);

    // --- 3. the REAL AZdecrypt worked example (asserted KAT) ------------------------
    {
        const char *C =
            "I E LTIIKI E  SKIOSLLM INCUG  HPEUFOP.NLEI  B TECSIAUM SOEMREI  F NUNHT T "
            "EHAOFN ERKITSLLB INCEG UAWIESLDM  GNAA  LISFO TA HELL M .OSOTT K DALING "
            "LEROSOUEMS HTANNIIM GA";
        const char *PT =
            "I LIKE KILLING PEOPLE BECAUSE IT IS SO MUCH FUN. IT IS MORE FUN THAN "
            "KILLING WILD GAME IN THE FOREST BECAUSE MAN IS THE MOST DANGEROUS ANIMAL "
            "OF ALL. TO KILL SOMETHING ";
        int len = encode_spaced(C, cipher);
        static int want[MAX_CIPHER_LENGTH];
        int wlen = encode_spaced(PT, want);
        CHECK(len == 168 && wlen == 168, "AZ e2e lengths (%d,%d)", len, wlen);
        period_column_search(cipher, len, 2, ngrams, NGRAM_SIZE,
                             NULL, NULL, 0, 1.0f, 0.0f, got, stages, &nstg);
        CHECK(arrays_equal(got, want, len), "AZ e2e recovers 'I LIKE KILLING PEOPLE ...'");
        // The recovered stages must be AZdecrypt's exact pair.
        int stages_ok = (nstg == 2 &&
            stages[0].dx == 4  && stages[0].period == 3 && stages[0].utp == 0 &&
            stages[1].dx == 56 && stages[1].period == 2 && stages[1].utp == 1);
        CHECK(stages_ok, "AZ e2e stages == [4x42,TP,P:3][56x3,UTP,P:2] (got %d stages: "
              "[%dx?,%s,P:%d][%dx?,%s,P:%d])", nstg,
              stages[0].dx, stages[0].utp ? "UTP" : "TP", stages[0].period,
              nstg > 1 ? stages[1].dx : 0, (nstg > 1 && stages[1].utp) ? "UTP" : "TP",
              nstg > 1 ? stages[1].period : 0);
    }

    // --- 4. length cliff (depth-2), printed + a long-length floor -------------------
    printf("\n  depth-2 recovery vs length (single stage 60x?,TP,P:3 -> solved blind):\n");
    int lens[] = {60, 90, 120, 168, 240};
    int nlens = (int) (sizeof lens / sizeof lens[0]);
    int long_ok = 0, long_total = 0;
    for (int li = 0; li < nlens; li++) {
        int L = lens[li];
        int len = encode_letters(PLAINTEXT, L, pt);
        // pick a non-trivial divisor width: smallest divisor >= 4 (falls back to 2).
        int dx = 2;
        for (int d = 4; d <= len; d++) { if (len % d == 0) { dx = d; break; } }
        int period = (dx >= 4) ? 3 : 2;
        period_column_transform(pt, cipher, len, dx, period, 0);
        period_column_search(cipher, len, 2, ngrams, NGRAM_SIZE,
                             NULL, NULL, 0, 1.0f, 0.0f, got, stages, &nstg);
        int ok = arrays_equal(got, pt, len);
        printf("    len %4d (dx=%d): %s\n", len, dx, ok ? "recovered" : "MISS");
        if (L >= 120) { long_total++; if (ok) long_ok++; }
    }
    CHECK(long_ok == long_total, "long-length (>=120) recovery %d/%d", long_ok, long_total);

    printf("\n%d checks, %d failures\n", checks, failures);
    free(ngrams);
    return failures ? 1 : 0;
}

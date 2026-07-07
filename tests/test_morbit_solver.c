//
//  In-process capability + calibration tests for the Morbit solver (morbit_search).
//
//  Framework-free: build with `make testopt`. The solver is DETERMINISTIC and EXHAUSTIVE
//  (it decodes + scores every one of the 9! = 362,880 pair<->digit bijections and keeps the
//  global best), so there is no anneal/shotgun/pso schedule to tune -- the only knob is the
//  Morse-validity weight. The "stress test" is therefore the RECOVERY-VS-LENGTH cliff (the
//  fundamental n-gram-discrimination limit) plus a validity-weight / -logprob sweep that
//  CALIBRATES that single knob, the honest analogue of the other types' -method calibration.
//
//  Morbit's encode is DETERMINISTIC (a bijection, unlike Pollux's polyphonic map), so the
//  ciphertext for a given (plaintext, keyword) is fixed -- there is no RNG-seed dimension to
//  sweep. We vary over several keywords and lengths instead.
//
//  Checks:
//    1. exact recovery across keywords at a solid length (asserted 100%);
//    2. a recovery-vs-length cliff (printed) with an asserted floor;
//    3. a validity-weight sweep at a short length (printed) -- the default MORBIT_VALID_WEIGHT
//       must be no worse than weight 0;
//    4. reward-only vs -logprob at a short length (printed characterization);
//    5. determinism (same input -> identical recovery).
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "scoring.h"          // load_ngrams
#include "morbit.h"           // morbit_encrypt
#include "morbit_solver.h"    // morbit_search / MORBIT_VALID_WEIGHT

static int failures = 0;
static int checks = 0;

#define CHECK(cond, ...) do { \
    checks++; \
    if (!(cond)) { failures++; printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } \
} while (0)

#define NGRAM_FILE "english_quadgrams.txt"
#define NGRAM_SIZE 4

// A long chunk of natural English (Pride and Prejudice, opening), letters only.
static const char *PLAINTEXT =
    "ITISATRUTHUNIVERSALLYACKNOWLEDGEDTHATASINGLEMANINPOSSESSIONOFAGOODFORTUNEMUSTBEINWANTOFAWIFE"
    "HOWEVERLITTLEKNOWNTHEFEELINGSORVIEWSOFSUCHAMANMAYBEONHISFIRSTENTERINGANEIGHBOURHOODTHISTRUTHIS"
    "SOWELLFIXEDINTHEMINDSOFTHESURROUNDINGFAMILIESTHATHEISCONSIDEREDTHERIGHTFULPROPERTYOFSOMEONEOR"
    "OTHEROFTHEIRDAUGHTERSMYDEARMRBENNETSAIDHISLADYTOHIMONEDAYHAVEYOUHEARDTHATNETHERFIELDPARKISLET";

// 9-letter keywords -> pair<->digit bijections. (Any 9-letter string works; ranks of 9 items
// always form a permutation of 1..9.)
static const char *KEYWORDS[] = { "WISECRACK", "DUPLICATE", "BLACKENED", "CHIVALROU" };

// Derive the digit(1..9) -> pair(0..8) key from a 9-letter keyword by stable alphabetical
// rank (the ACA convention -- the same derivation as tools/morbit_gen.c).
static void key_from_keyword(const char *kw, int key[10]) {
    key[0] = 0;
    for (int i = 0; i < 9; i++) {
        int rank = 1;
        for (int j = 0; j < 9; j++)
            if (kw[j] < kw[i] || (kw[j] == kw[i] && j < i)) rank++;
        key[rank] = i;
    }
}

// Encipher the first L letters of PLAINTEXT under keyword kw, run the exhaustive search at
// validity weight vw over the n-gram table `ng`, and return the fraction of plaintext
// positions recovered (1.0 == exact). *full receives 1 iff every position (and the length) match.
static double try_solve(int L, const char *kw, double vw, float *ng, int *full) {
    static int pt[MAX_CIPHER_LENGTH];
    int n = 0;
    for (int i = 0; i < L && PLAINTEXT[i]; i++)
        pt[n++] = g_char_to_idx[toupper((unsigned char) PLAINTEXT[i])];

    int key[10]; key_from_keyword(kw, key);
    static int cipher[MAX_CIPHER_LENGTH];
    int clen = morbit_encrypt(pt, n, key, cipher);

    static int bpt[MAX_CIPHER_LENGTH];
    int bk[10], bn = 0, nt = 0, nv = 0;
    morbit_search(cipher, clen, ng, NGRAM_SIZE, 1.0f, vw, bk, bpt, &bn, &nt, &nv);

    int m = 0;
    if (bn == n) for (int i = 0; i < n; i++) if (bpt[i] == pt[i]) m++;
    if (full) *full = (bn == n && m == n);
    return n ? (double) m / n : 0.0;
}

int main(void) {
    init_alphabet(NULL);                    // full 26-letter alphabet (no J->I)
    float *ng = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!ng) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    int nkw = (int) (sizeof KEYWORDS / sizeof KEYWORDS[0]);

    // --- 1. exact recovery across keywords at L=120 ---------------------------------
    {
        int ok = 0, total = 0;
        for (int ki = 0; ki < nkw; ki++) {
            int full = 0;
            try_solve(120, KEYWORDS[ki], MORBIT_VALID_WEIGHT, ng, &full);
            total++; ok += full;
        }
        CHECK(ok == total, "exact recovery at L=120: %d/%d", ok, total);
    }

    // --- 2. recovery-vs-length cliff (printed down to the limit; asserted floor) -----
    // Morbit's keyspace is 9! (= ~6x Pollux's 3^10), so exhaustive n-gram scoring recovers
    // from short text but the cliff sits a little higher than Pollux's. We print it and
    // assert a solid floor from L=50 up.
    printf("\n  Morbit recovery vs plaintext length (mean over %d keywords,"
           " valid weight %.1f):\n", nkw, MORBIT_VALID_WEIGHT);
    {
        int lens[] = {16, 20, 24, 30, 40, 50, 70, 100};
        int nlens = (int) (sizeof lens / sizeof lens[0]);
        double acc_at_50 = 0.0, acc_at_100 = 0.0;
        for (int li = 0; li < nlens; li++) {
            int L = lens[li];
            double sum = 0.0; int full_ct = 0, cnt = 0;
            for (int ki = 0; ki < nkw; ki++) {
                int full = 0;
                sum += try_solve(L, KEYWORDS[ki], MORBIT_VALID_WEIGHT, ng, &full);
                full_ct += full; cnt++;
            }
            double mean = sum / cnt;
            printf("    len %4d: mean acc %.3f, exact %d/%d\n", L, mean, full_ct, cnt);
            if (L == 50) acc_at_50 = mean;
            if (L == 100) acc_at_100 = mean;
        }
        CHECK(acc_at_50 > 0.90, "recovery at L=50 mean acc %.3f (want > 0.90)", acc_at_50);
        CHECK(acc_at_100 > 0.99, "recovery at L=100 mean acc %.3f (want > 0.99)", acc_at_100);
    }

    // --- 3. validity-weight calibration sweep in the MARGINAL regime -----------------
    // At the cliff (short text) the validity reward should help; here we sweep it so the
    // default MORBIT_VALID_WEIGHT is calibrated to do at least as well as weight 0.
    printf("\n  validity-weight sweep at L=16 (marginal; mean acc over %d keywords):\n", nkw);
    {
        double weights[] = {0.0, 1.0, 3.0, 6.0};
        int nw = (int) (sizeof weights / sizeof weights[0]);
        double acc0 = -1.0, acc_default = -1.0;
        for (int wi = 0; wi < nw; wi++) {
            double sum = 0.0; int cnt = 0;
            for (int ki = 0; ki < nkw; ki++) { sum += try_solve(16, KEYWORDS[ki], weights[wi], ng, NULL); cnt++; }
            double mean = sum / cnt;
            printf("    valid weight %.1f: mean acc %.3f\n", weights[wi], mean);
            if (weights[wi] == 0.0) acc0 = mean;
            if (weights[wi] == MORBIT_VALID_WEIGHT) acc_default = mean;
        }
        CHECK(acc_default >= acc0 - 1e-9,
              "default validity weight %.1f no worse than weight 0 at L=16 (%.3f vs %.3f)",
              MORBIT_VALID_WEIGHT, acc_default, acc0);
    }

    // --- 4. reward-only vs -logprob at the cliff (characterization) -----------------
    printf("\n  reward-only vs -logprob at L=20 (mean acc over %d keywords):\n", nkw);
    {
        double ro = 0.0; int cnt = 0;
        for (int ki = 0; ki < nkw; ki++) { ro += try_solve(20, KEYWORDS[ki], MORBIT_VALID_WEIGHT, ng, NULL); cnt++; }
        printf("    reward-only: %.3f\n", ro / cnt);

        g_ngram_logprob = true;
        float *lp = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
        if (lp) {
            double s = 0.0; int c2 = 0;
            for (int ki = 0; ki < nkw; ki++) { s += try_solve(20, KEYWORDS[ki], MORBIT_VALID_WEIGHT, lp, NULL); c2++; }
            printf("    -logprob   : %.3f\n", s / c2);
            free(lp);
        }
        g_ngram_logprob = false;
    }

    // --- 5. determinism -------------------------------------------------------------
    {
        int f1 = 0, f2 = 0;
        double a1 = try_solve(120, KEYWORDS[0], MORBIT_VALID_WEIGHT, ng, &f1);
        double a2 = try_solve(120, KEYWORDS[0], MORBIT_VALID_WEIGHT, ng, &f2);
        CHECK(a1 == a2 && f1 == f2, "deterministic: (%.4f,%d) vs (%.4f,%d)", a1, f1, a2, f2);
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(ng);
    return failures ? 1 : 0;
}

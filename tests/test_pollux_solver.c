//
//  In-process capability + calibration tests for the Pollux solver (pollux_search).
//
//  Framework-free: build with `make testopt`. The solver is DETERMINISTIC and EXHAUSTIVE
//  (it decodes + scores every one of the 3^10 = 59049 digit->{dot,dash,x} keys and keeps the
//  global best), so there is no anneal/shotgun/pso schedule to tune -- the only knob is the
//  Morse-validity weight. The "stress test" is therefore the RECOVERY-VS-LENGTH cliff (the
//  fundamental n-gram-discrimination limit) plus a validity-weight / -logprob sweep that
//  CALIBRATES that single knob, the honest analogue of the other types' -method calibration.
//
//  Checks:
//    1. exact recovery across assignments x seeds at a solid length (asserted 100%);
//    2. a recovery-vs-length cliff (printed) with an asserted floor at the ACA ~80-100 band;
//    3. a validity-weight sweep at a short length (printed) -- the default POLLUX_VALID_WEIGHT
//       must recover the sweep set;
//    4. reward-only vs -logprob at a short length (printed characterization);
//    5. determinism (same input -> identical recovery).
//
//  Run from the source directory so the n-gram table is found in the cwd.
//

#include "colossus.h"
#include "scoring.h"          // load_ngrams
#include "pollux.h"           // pollux_encrypt
#include "pollux_solver.h"    // pollux_search / POLLUX_VALID_WEIGHT

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

// Assignment strings (each 4 dots / 3 dashes / 3 x, a valid Pollux key).
static const char *ASSIGN[] = {
    ".x-..x.--x", "x.-.x-.x-.", "-.x-.x-.x.", "..-x.-x.-x"
};

static void key_from_str(const char *a, int key[10]) {
    for (int d = 0; d < 10; d++)
        key[d] = (a[d] == '.') ? PX_DOT : (a[d] == '-') ? PX_DASH : PX_X;
}

// Encipher the first L letters of PLAINTEXT under assignment a (RNG seed s), run the exhaustive
// search at validity weight vw over the n-gram table `ng`, and return the fraction of plaintext
// positions recovered (1.0 == exact). *full receives 1 iff every position (and the length) match.
static double try_solve(int L, const char *a, unsigned s, double vw, float *ng, int *full) {
    static int pt[MAX_CIPHER_LENGTH];
    int n = 0;
    for (int i = 0; i < L && PLAINTEXT[i]; i++)
        pt[n++] = g_char_to_idx[toupper((unsigned char) PLAINTEXT[i])];

    int key[10]; key_from_str(a, key);
    seed_rand(s);
    static int cipher[MAX_CIPHER_LENGTH];
    int clen = pollux_encrypt(pt, n, key, cipher);

    static int bpt[MAX_CIPHER_LENGTH];
    int bk[10], bn = 0, nt = 0, nv = 0;
    pollux_search(cipher, clen, ng, NGRAM_SIZE, 1.0f, vw, bk, bpt, &bn, &nt, &nv);

    int m = 0;
    if (bn == n) for (int i = 0; i < n; i++) if (bpt[i] == pt[i]) m++;
    if (full) *full = (bn == n && m == n);
    return n ? (double) m / n : 0.0;
}

int main(void) {
    init_alphabet(NULL);                    // full 26-letter alphabet (no J->I)
    float *ng = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
    if (!ng) { printf("FAIL: could not load %s\n", NGRAM_FILE); return 1; }

    int nassign = (int) (sizeof ASSIGN / sizeof ASSIGN[0]);
    unsigned seeds[] = {1u, 7u, 42u};

    // --- 1. exact recovery across assignments x seeds at L=120 ----------------------
    {
        int ok = 0, total = 0;
        for (int ai = 0; ai < nassign; ai++)
            for (int si = 0; si < 2; si++) {
                int full = 0;
                try_solve(120, ASSIGN[ai], seeds[si], POLLUX_VALID_WEIGHT, ng, &full);
                total++; ok += full;
            }
        CHECK(ok == total, "exact recovery at L=120: %d/%d", ok, total);
    }

    // --- 2. recovery-vs-length cliff (printed down to the limit; asserted floor) ----
    // Pollux's keyspace is tiny (3^10), so exhaustive n-gram scoring recovers from very
    // short text -- the cliff sits well BELOW the ACA 80-100 range. We print it down to a
    // dozen letters (the real limit) and assert a solid floor from L=40 up.
    printf("\n  Pollux recovery vs plaintext length (mean over %d assignments x seeds,"
           " valid weight %.1f):\n", nassign, POLLUX_VALID_WEIGHT);
    {
        int lens[] = {12, 16, 20, 24, 30, 40, 60, 90};
        int nlens = (int) (sizeof lens / sizeof lens[0]);
        double acc_at_40 = 0.0, acc_at_90 = 0.0;
        for (int li = 0; li < nlens; li++) {
            int L = lens[li];
            double sum = 0.0; int full_ct = 0, cnt = 0;
            for (int ai = 0; ai < nassign; ai++)
                for (int si = 0; si < 3; si++) {
                    int full = 0;
                    sum += try_solve(L, ASSIGN[ai], seeds[si], POLLUX_VALID_WEIGHT, ng, &full);
                    full_ct += full; cnt++;
                }
            double mean = sum / cnt;
            printf("    len %4d: mean acc %.3f, exact %d/%d\n", L, mean, full_ct, cnt);
            if (L == 40) acc_at_40 = mean;
            if (L == 90) acc_at_90 = mean;
        }
        CHECK(acc_at_40 > 0.95, "recovery at L=40 mean acc %.3f (want > 0.95)", acc_at_40);
        CHECK(acc_at_90 > 0.99, "recovery at L=90 mean acc %.3f (want > 0.99)", acc_at_90);
    }

    // --- 3. validity-weight calibration sweep in the MARGINAL regime -----------------
    // At the cliff (short text) the validity reward should help; here we sweep it so the
    // default POLLUX_VALID_WEIGHT is calibrated to do at least as well as weight 0.
    printf("\n  validity-weight sweep at L=16 (marginal; mean acc over %d cases):\n",
           nassign * 3);
    {
        double weights[] = {0.0, 1.0, 3.0, 6.0};
        int nw = (int) (sizeof weights / sizeof weights[0]);
        double acc0 = -1.0, acc_default = -1.0;
        for (int wi = 0; wi < nw; wi++) {
            double sum = 0.0; int cnt = 0;
            for (int ai = 0; ai < nassign; ai++)
                for (int si = 0; si < 3; si++) {
                    sum += try_solve(16, ASSIGN[ai], seeds[si], weights[wi], ng, NULL); cnt++;
                }
            double mean = sum / cnt;
            printf("    valid weight %.1f: mean acc %.3f\n", weights[wi], mean);
            if (weights[wi] == 0.0) acc0 = mean;
            if (weights[wi] == POLLUX_VALID_WEIGHT) acc_default = mean;
        }
        CHECK(acc_default >= acc0 - 1e-9,
              "default validity weight %.1f no worse than weight 0 at L=16 (%.3f vs %.3f)",
              POLLUX_VALID_WEIGHT, acc_default, acc0);
    }

    // --- 4. reward-only vs -logprob at the cliff (characterization) -----------------
    printf("\n  reward-only vs -logprob at L=20 (mean acc over %d cases):\n", nassign * 3);
    {
        double ro = 0.0; int cnt = 0;
        for (int ai = 0; ai < nassign; ai++)
            for (int si = 0; si < 3; si++) { ro += try_solve(20, ASSIGN[ai], seeds[si], POLLUX_VALID_WEIGHT, ng, NULL); cnt++; }
        printf("    reward-only: %.3f\n", ro / cnt);

        g_ngram_logprob = true;
        float *lp = load_ngrams(NGRAM_FILE, NGRAM_SIZE, false);
        if (lp) {
            double s = 0.0; int c2 = 0;
            for (int ai = 0; ai < nassign; ai++)
                for (int si = 0; si < 3; si++) { s += try_solve(20, ASSIGN[ai], seeds[si], POLLUX_VALID_WEIGHT, lp, NULL); c2++; }
            printf("    -logprob   : %.3f\n", s / c2);
            free(lp);
        }
        g_ngram_logprob = false;
    }

    // --- 5. determinism -------------------------------------------------------------
    {
        int f1 = 0, f2 = 0;
        double a1 = try_solve(120, ASSIGN[0], 42u, POLLUX_VALID_WEIGHT, ng, &f1);
        double a2 = try_solve(120, ASSIGN[0], 42u, POLLUX_VALID_WEIGHT, ng, &f2);
        CHECK(a1 == a2 && f1 == f2, "deterministic: (%.4f,%d) vs (%.4f,%d)", a1, f1, a2, f2);
    }

    printf("\n%d checks, %d failures\n", checks, failures);
    free(ng);
    return failures ? 1 : 0;
}

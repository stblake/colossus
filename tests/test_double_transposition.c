// Unit tests for the double columnar transposition divide-and-conquer core:
// the bigram table (dct_load_bigrams) and, above all, the Index of Digraphic
// Potential (dct_idp) -- the fitness that scores the second key K2 WITHOUT knowing
// K1. We validate the paper's two key properties (Lasry/Kopal/Wacker 2014, Fig 2):
//   (1) the IDP is SELECTIVE  -- the true K2 out-scores random keys, and is the
//       per-edge maximum across L1 hypotheses at the true L1; and
//   (2) the IDP is MONOTONIC  -- it degrades as the key is progressively perturbed.
// Plus a round-trip (double encrypt then decrypt is the identity) and a greedy-
// reconstruction check (undoing the true K2, the IDP's own column chaining
// reconstructs the single K1 transposition, recovering the plaintext).
//
// Links the real solver core + primitives, so the test and the solver share code.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "colossus.h"
#include "double_transposition_solver.h"
#include "scoring.h"

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { \
    printf("  FAIL: %s\n", msg); failures++; } } while (0)

// ---- local double-transposition encryption (mirror of the generator) -----

static int keyword_order(const char *kw, int order[]) {
    int lets[MAX_COLS], K = 0;
    for (const char *k = kw; *k && K < MAX_COLS; k++) {
        int c = toupper((unsigned char) *k);
        if (c >= 'A' && c <= 'Z') lets[K++] = c - 'A';
    }
    int used[MAX_COLS];
    for (int i = 0; i < K; i++) used[i] = 0;
    for (int j = 0; j < K; j++) {
        int best = -1;
        for (int c = 0; c < K; c++)
            if (!used[c] && (best < 0 || lets[c] < lets[best])) best = c;
        used[best] = 1; order[j] = best;
    }
    return K;
}

static void columnar_encrypt(const int in[], int len, int K, const int order[], int out[]) {
    static int grid[MAX_CIPHER_LENGTH];
    int R = (len + K - 1) / K, rem = len % K, pos = 0;
    for (int r = 0; r < R; r++)
        for (int c = 0; c < K; c++) {
            int h = (rem == 0 || c < rem) ? R : R - 1;
            if (r < h) grid[r * K + c] = in[pos++];
        }
    int o = 0;
    for (int j = 0; j < K; j++) {
        int c = order[j], h = (rem == 0 || c < rem) ? R : R - 1;
        for (int r = 0; r < h; r++) out[o++] = grid[r * K + c];
    }
}

static int text_to_idx(const char *s, int idx[]) {
    int n = 0;
    for (const char *p = s; *p; p++) {
        int c = toupper((unsigned char) *p);
        if (c >= 'A' && c <= 'Z') idx[n++] = c - 'A';
    }
    return n;
}

// A chunk of natural English (double-transposition-challenge flavour, ~560 letters).
static const char *PLAIN =
    "THEDOUBLETRANSPOSITIONCIPHERHASBEENONEOFTHEMOSTPOPULARMANUALCIPHERSIT"
    "DIDNOTREQUIRETHEUSEOFADEVICEFORENCRYPTIONANDDECRYPTIONBECAUSEOFITSSIM"
    "PLICITYANDITSHIGHLEVELOFSECURITYITWASOFTENTHECIPHEROFCHOICEFORINTELLI"
    "GENCEANDSECRETOPERATIONSORGANIZATIONSTHEPROCESSOFENCRYPTIONANDDECRYPT"
    "IONISRELATIVELYSIMPLEFIRSTTWOTRANSPOSITIONKEYSMUSTBECHOSENANDAGREEDIN"
    "ADVANCEKEYSAREUSUALLYDERIVEDFROMKEYWORDSORKEYPHRASESWHICHINTURNARECON"
    "VERTEDTONUMERICALKEYSSUCHKEYWORDSOREXPRESSIONSAREOFTENTAKENFROMBOOKSO"
    "RNEWSPAPERSANDMUSTBEKEPTSECRETBYBOTHPARTIESATALLTIMESDURINGOPERATIONS";

int main(void) {
    printf("test_double_transposition:\n");
    init_alphabet(NULL);
    rng_state = 20140711ULL;

    static int P[MAX_CIPHER_LENGTH];
    int n = text_to_idx(PLAIN, P);
    printf("  plaintext length = %d\n", n);

    // --- bigram table sanity -------------------------------------------------
    double *bg = dct_load_bigrams("english_quadgrams.txt", 4);
    CHECK(bg != NULL, "dct_load_bigrams returned non-NULL");
    if (!bg) { printf("  (cannot continue without english_quadgrams.txt)\n"); return 1; }
    int TH = ('T' - 'A') * 26 + ('H' - 'A');
    int QZ = ('Q' - 'A') * 26 + ('Z' - 'A');
    int HE = ('H' - 'A') * 26 + ('E' - 'A');
    CHECK(bg[TH] > bg[QZ], "bigram TH scores higher than QZ");
    CHECK(bg[HE] > bg[QZ], "bigram HE scores higher than QZ");

    // --- build a double-transposition cipher --------------------------------
    int o1[MAX_COLS], o2[MAX_COLS];
    int L1 = keyword_order("SUBMARINE", o1);   // 9 columns
    int L2 = keyword_order("WATERFALLS", o2);  // 10 columns (distinct-letter positions)
    static int mid[MAX_CIPHER_LENGTH], C[MAX_CIPHER_LENGTH];
    columnar_encrypt(P,   n, L1, o1, mid);
    columnar_encrypt(mid, n, L2, o2, C);
    printf("  L1=%d L2=%d\n", L1, L2);

    // --- round-trip: decrypt(C) == P ----------------------------------------
    static int back_mid[MAX_CIPHER_LENGTH], back[MAX_CIPHER_LENGTH];
    decrypt_columnar(C, n, L2, o2, COL_READ_TB, back_mid);
    decrypt_columnar(back_mid, n, L1, o1, COL_READ_TB, back);
    int rt_ok = 1;
    for (int i = 0; i < n; i++) if (back[i] != P[i]) { rt_ok = 0; break; }
    CHECK(rt_ok, "double encrypt -> decrypt round-trips to the plaintext");

    static int work[MAX_CIPHER_LENGTH];

    // --- (1) SELECTIVITY: true K2 beats random, and is the max across L1 -----
    double idp_true = dct_idp(C, n, L1, L2, COL_READ_TB, o2, bg, work, NULL);
    double sum = 0.0, best_rand = -1e30;
    int NR = 400;
    for (int r = 0; r < NR; r++) {
        int p[MAX_COLS];
        for (int i = 0; i < L2; i++) p[i] = i;
        shuffle(p, L2);
        double v = dct_idp(C, n, L1, L2, COL_READ_TB, p, bg, work, NULL);
        sum += v; if (v > best_rand) best_rand = v;
    }
    double mean_rand = sum / NR;
    printf("  IDP(true)=%.4f  random mean=%.4f  random best(of %d)=%.4f\n",
        idp_true, mean_rand, NR, best_rand);
    CHECK(idp_true > mean_rand + 0.05, "IDP(true K2) clearly exceeds the random mean");
    CHECK(idp_true > best_rand,         "IDP(true K2) is the maximum over 400 random keys");

    // The per-edge normalization must make the TRUE L1 the best across L1 --
    // an unnormalized sum would trivially favour the smallest L1.
    double idp_at_trueL1 = idp_true;
    int best_L1 = L1; double best_L1_val = idp_true;
    printf("  IDP(true K2) vs L1: ");
    for (int Lx = L1 - 2; Lx <= L1 + 2; Lx++) {
        if (Lx < 2) continue;
        double v = dct_idp(C, n, Lx, L2, COL_READ_TB, o2, bg, work, NULL);
        printf("L1=%d:%.3f ", Lx, v);
        if (v > best_L1_val) { best_L1_val = v; best_L1 = Lx; }
    }
    printf("\n");
    CHECK(best_L1 == L1, "IDP(true K2) peaks at the true L1 (per-edge normalized)");
    (void) idp_at_trueL1;

    // --- (2) MONOTONICITY: IDP degrades as K2 is perturbed -------------------
    // Average the IDP of the true key with x random transposition-swaps applied,
    // for x = 0, 2, 6; the means must strictly decrease (paper Fig 2).
    double mono[3]; int xs[3] = {0, 2, 6};
    for (int t = 0; t < 3; t++) {
        int trials = 60; double acc = 0.0;
        for (int tr = 0; tr < trials; tr++) {
            int p[MAX_COLS];
            for (int i = 0; i < L2; i++) p[i] = o2[i];
            for (int sw = 0; sw < xs[t]; sw++) {
                int a = rand_int(0, L2), b = rand_int(0, L2);
                int tmp = p[a]; p[a] = p[b]; p[b] = tmp;
            }
            acc += dct_idp(C, n, L1, L2, COL_READ_TB, p, bg, work, NULL);
        }
        mono[t] = acc / trials;
    }
    printf("  IDP mean by perturbations: x=0:%.4f  x=2:%.4f  x=6:%.4f\n",
        mono[0], mono[1], mono[2]);
    CHECK(mono[0] > mono[1] && mono[1] > mono[2], "IDP degrades monotonically with perturbation");

    // --- greedy reconstruction: true K2 -> recover the K1 column order -------
    // The chain reconstructs the column adjacency, but a columnar transposition is
    // recoverable only up to a cyclic rotation (free path endpoints), so score the
    // BEST rotation of the chain -- which cuts it at the true column boundary.
    int chain[MAX_COLS];
    dct_idp(C, n, L1, L2, COL_READ_TB, o2, bg, work, chain);
    static int rec_pt[MAX_CIPHER_LENGTH];
    int order[MAX_COLS];
    double best_frac = 0.0;
    for (int rot = 0; rot < L1; rot++) {
        for (int g = 0; g < L1; g++) order[chain[(g + rot) % L1]] = g;
        decrypt_columnar(back_mid, n, L1, order, COL_READ_TB, rec_pt);
        int match = 0;
        for (int i = 0; i < n; i++) if (rec_pt[i] == P[i]) match++;
        double frac = (double) match / n;
        if (frac > best_frac) best_frac = frac;
    }
    printf("  greedy K1 reconstruction from true K2 (best rotation): %.1f%% plaintext match\n",
        100.0 * best_frac);
    CHECK(best_frac > 0.90, "greedy digraph chaining recovers the K1 order (up to rotation) from the true K2");

    // --- crib guidance of K2: the true K2's chain rotation-resolves to a crib-matching
    // K1, while random K2 keys do not -- the signal the crib-aware K2 fitness exploits to
    // steer the K2 climb toward the true key (a deterministic proof of crib guidance). --
    int crib_pos[64], crib_idx[64], n_cr = 0;               // first 40 plaintext letters
    for (int i = 0; i < 40 && i < n; i++) { crib_pos[n_cr] = i; crib_idx[n_cr] = P[i]; n_cr++; }
    // Best crib fraction over the chain rotations, for a given K2 (what dct_k2_fitness does).
    #define BEST_CRIB(K2ord) ({                                                            \
        int ch[MAX_COLS]; dct_idp(C, n, L1, L2, COL_READ_TB, (K2ord), bg, work, ch);        \
        int *II = work; double bc = 0.0; int ord[MAX_COLS];                                 \
        for (int rr = 0; rr < L1; rr++) {                                                   \
            for (int g = 0; g < L1; g++) ord[ch[(g + rr) % L1]] = g;                        \
            static int dd[MAX_CIPHER_LENGTH]; decrypt_columnar(II, n, L1, ord, COL_READ_TB, dd); \
            double cf = crib_score(dd, n, crib_idx, crib_pos, n_cr); if (cf > bc) bc = cf;  \
        } bc; })
    double crib_true = BEST_CRIB(o2);
    double crib_rand_sum = 0.0, crib_rand_max = 0.0;
    for (int r = 0; r < 50; r++) {
        int p[MAX_COLS]; for (int i = 0; i < L2; i++) p[i] = i; shuffle(p, L2);
        double v = BEST_CRIB(p); crib_rand_sum += v; if (v > crib_rand_max) crib_rand_max = v;
    }
    printf("  crib fraction (best rotation): true K2=%.3f  random mean=%.3f max=%.3f\n",
        crib_true, crib_rand_sum / 50, crib_rand_max);
    CHECK(crib_true > 0.90, "true K2 rotation-resolves to a crib-matching K1 (crib guides K2)");
    CHECK(crib_true > crib_rand_max + 0.3, "true K2's crib match dominates random K2 keys");

    free(bg);
    if (failures == 0) printf("  all double-transposition primitive tests passed.\n");
    else               printf("  %d CHECK(S) FAILED.\n", failures);
    return failures ? 1 : 0;
}

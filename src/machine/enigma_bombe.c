// =====================================================================
//  Enigma Turing-Welchman Bombe (known-plaintext / crib attack)
// =====================================================================
//
// Given a crib (known plaintext aligned to a ciphertext stretch), the Bombe recovers the
// rotor config (order + fast ring + start position) far more cheaply than a ciphertext-only
// IoC search, then hands the winner to the shared ring-refine + plugboard-completion path.
//
// Menu.  Each crib position j fixes PT=p_j and CT=c_j. Enigma applies the plugboard at both
// ends, so with S_j the (plugboard-excluded) scrambler permutation at that position:
//     stecker(c_j) = S_j( stecker(p_j) )
// S_j is an involution, so the relation is symmetric. The crib positions form a graph (the
// "menu") whose nodes are letters and whose edges are the scrambler relations.
//
// Search.  For each rotor order x fast-rotor ring x start position (the 60 x 26^4 "locations"
// of Gillogly / Ostwald-Weierud -- the fast ring is searched because a middle-rotor turnover
// inside the crib window depends on it):
//   * build the scramblers S_j (one forward pass over the crib span);
//   * for each hypothesis of the most-connected letter's stecker, PROPAGATE the relation
//     through the menu with the WELCHMAN DIAGONAL BOARD (the involution: setting
//     stecker(x)=y forces stecker(y)=x). A contradiction (a letter forced to two values, or
//     the involution violated) kills the hypothesis; a survivor is a "stop".
//   * score the stop by how many crib letters its partial stecker decrypts correctly.
// The best stops (top-M by crib match) become base keys; enigma_attack_from_bases() refines
// the remaining rings and completes/polishes the plugboard (seeded by the Bombe's steckers)
// under the n-gram (+ crib) fitness and reports the winner.
//
// This recovers all steckers "within seconds" from a good crib on a known/short wheel-order
// set (Ostwald-Weierud). A fully blind (60-order) search is ~26^4x60 locations, parallelised
// across -nthreads over the order x fast-ring work. M4 requires a pinned wheel order.

#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "enigma_solver.h"

#define BOMBE_MAX_EDGES 256
#define BOMBE_TOP_M     8      // stops carried to the ring/plugboard completion

typedef struct { int p, c, pos; } BombeEdge;      // menu edge: PT/CT letters at a position
typedef struct { EnigmaKey base; int matches; } BombeStop;

// Set stecker[x]=y enforcing the involution (diagonal board). Returns false on contradiction.
static inline bool steck_set(int steck[26], int x, int y) {
    if (steck[x] != -1) return steck[x] == y;
    if (steck[y] != -1 && steck[y] != x) return false;
    steck[x] = y;
    steck[y] = x;
    return true;
}

// Propagate the menu constraints from an initial assignment; false on contradiction.
static bool bombe_propagate(int steck[26], const BombeEdge *e, int ne, int (*scr)[26]) {
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int j = 0; j < ne; j++) {
            int a = e[j].p, b = e[j].c, va = steck[a], vb = steck[b];
            if (va != -1 && vb == -1)      { if (!steck_set(steck, b, scr[j][va])) return false; changed = 1; }
            else if (vb != -1 && va == -1) { if (!steck_set(steck, a, scr[j][vb])) return false; changed = 1; }
            else if (va != -1 && vb != -1) { if (scr[j][va] != vb) return false; }
        }
    }
    return true;
}

static void bombe_keep(BombeStop *top, int *ntop, const BombeStop *s, int M) {
    if (*ntop < M) { top[(*ntop)++] = *s; return; }
    int worst = 0;
    for (int i = 1; i < M; i++) if (top[i].matches < top[worst].matches) worst = i;
    if (s->matches > top[worst].matches) top[worst] = *s;
}

// Shared read-only search inputs + one worker's (order x fast-ring) unit sub-range. A unit u
// maps to order = u / n_rr and fast ring = rr_base + (u % n_rr); splitting over units lets a
// pinned wheel order (n_orders == 1) still parallelise across its 26 ring settings.
typedef struct {
    const ColossusConfig *cfg;
    const EnigmaKey *tmpl;
    const BombeEdge *edge;
    const int *want;        // want[pos] = edge index, else -1
    int ne, maxpos, test, cipher_len;
    int n_rr, rr_base;      // fast-ring count (26 or 1) and base value
    int u_lo, u_hi;         // unit range for this worker
    BombeStop top[BOMBE_TOP_M];
    int ntop;
} BombeWork;

static void *bombe_worker(void *arg) {
    BombeWork *w = (BombeWork *) arg;
    const ColossusConfig *cfg = w->cfg;
    w->ntop = 0;
    // Per-unit scrambler table (M3): S[window_triple*26 + c] = the identity-plug scrambler at
    // each (left,middle,fast) window. Built ONCE per (order x fast-ring) unit, so each start's
    // crib-span scramblers are gathered by a cheap row copy instead of ne*26 encipherments.
    // M4's 26^4 table (47 MB) is not built -- nw==4 (and pos-pinned) fall back to the direct build.
    int *S = NULL;

    for (int u = w->u_lo; u < w->u_hi; u++) {
        int o  = u / w->n_rr;
        int rr = w->rr_base + (u % w->n_rr);
        EnigmaKey kbase = w->tmpl[o];
        int nw = kbase.n_wheels;
        for (int i = 0; i < nw; i++) kbase.ring[i] = 0;
        if (cfg->enigma_ring_present) {
            int off = nw - 3;
            for (int i = 0; i < 3; i++) kbase.ring[off + i] = cfg->enigma_ring[i];
        }
        kbase.ring[nw - 1] = rr;
        enigma_plug_identity(kbase.plug);
        long total = 1;
        for (int i = 0; i < nw; i++) total *= 26;

        bool use_table = (nw == 3 && !cfg->enigma_pos_present);
        if (use_table) {
            if (!S) S = (int *) malloc((size_t) 26 * 26 * 26 * 26 * sizeof(int));
            EnigmaKey kt = kbase;         // identity plug
            for (int a = 0; a < 26; a++)
                for (int b = 0; b < 26; b++)
                    for (int c = 0; c < 26; c++) {
                        kt.pos[0] = a; kt.pos[1] = b; kt.pos[2] = c;
                        int base = ((a * 26 + b) * 26 + c) * 26;
                        for (int ch = 0; ch < 26; ch++) S[base + ch] = enigma_encipher_letter(&kt, ch);
                    }
        }

        for (long pp = 0; pp < total; pp++) {
            EnigmaKey k = kbase;
            long q = pp;
            for (int i = nw - 1; i >= 0; i--) { k.pos[i] = (int)(q % 26); q /= 26; }
            if (cfg->enigma_pos_present) {
                int off = nw - 3;
                if (nw == 4) k.pos[0] = 0;
                for (int i = 0; i < 3; i++) k.pos[off + i] = cfg->enigma_pos[i];
            }

            // One-pass scrambler gather across the crib span (row copy from S, or direct build).
            int scr[BOMBE_MAX_EDGES][26];
            EnigmaKey kk = k;              // plugboard identity
            for (int s = 0; s <= w->maxpos; s++) {
                enigma_step(&kk);          // stepped s+1 times == crib position s
                int e = w->want[s];
                if (e >= 0) {
                    if (use_table)
                        memcpy(scr[e], &S[((kk.pos[0] * 26 + kk.pos[1]) * 26 + kk.pos[2]) * 26],
                               26 * sizeof(int));
                    else
                        for (int c = 0; c < 26; c++) scr[e][c] = enigma_encipher_letter(&kk, c);
                }
            }

            int best_match = -1, best_steck[26];
            for (int h = 0; h < 26; h++) {
                int steck[26];
                for (int i = 0; i < 26; i++) steck[i] = -1;
                if (!steck_set(steck, w->test, h)) continue;
                if (!bombe_propagate(steck, w->edge, w->ne, scr)) continue;
                int plug[26];
                for (int i = 0; i < 26; i++) plug[i] = (steck[i] >= 0) ? steck[i] : i;
                int m = 0;
                for (int j = 0; j < w->ne; j++)
                    if (plug[scr[j][plug[w->edge[j].c]]] == w->edge[j].p) m++;
                if (m > best_match) { best_match = m; memcpy(best_steck, steck, sizeof(steck)); }
            }
            if (best_match > 0) {
                BombeStop s;
                s.base = k;
                for (int i = 0; i < 26; i++) s.base.plug[i] = (best_steck[i] >= 0) ? best_steck[i] : i;
                s.matches = best_match;
                bombe_keep(w->top, &w->ntop, &s, BOMBE_TOP_M);
            }
            if (cfg->enigma_pos_present) break;
        }
    }
    if (S) free(S);
    return NULL;
}

bool solve_enigma_bombe(ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    enigma_init();

    // Build the menu edges from the crib.
    BombeEdge edge[BOMBE_MAX_EDGES];
    int ne = 0, deg[26] = {0};
    for (int i = 0; i < n_cribs && ne < BOMBE_MAX_EDGES; i++) {
        int pos = crib_positions[i];
        if (pos < 0 || pos >= cipher_len) continue;
        edge[ne].p = crib_indices[i];
        edge[ne].c = cipher_indices[pos];
        edge[ne].pos = pos;
        deg[edge[ne].p]++; deg[edge[ne].c]++;
        ne++;
    }
    if (ne < 3) {
        printf("\n\nERROR: Enigma Bombe needs a crib of at least 3 letters (got %d).\n\n", ne);
        return false;
    }
    int test = 0;
    for (int i = 1; i < 26; i++) if (deg[i] > deg[test]) test = i;

    int maxplugs = (cfg->enigma_maxplugs > 0) ? cfg->enigma_maxplugs : 10;
    if (maxplugs > 13) maxplugs = 13;

    if (cfg->enigma_model == 4 && !cfg->enigma_rotors_present) {
        printf("\n\nERROR: M4 Bombe needs a pinned wheel order (-rotors); 26^4 x hundreds "
               "of orders is impractical.\n\n");
        return false;
    }

    static EnigmaKey tmpl[64];
    int n_orders = enigma_enumerate_wheel_orders(cfg, tmpl, 64);

    // Positions needing a scrambler (unique crib positions), for the one-pass build.
    static int want[MAX_CIPHER_LENGTH];
    int maxpos = 0;
    for (int i = 0; i < cipher_len; i++) want[i] = -1;
    for (int j = 0; j < ne; j++) { want[edge[j].pos] = j; if (edge[j].pos > maxpos) maxpos = edge[j].pos; }

    int n_rr = cfg->enigma_ring_present ? 1 : 26;
    int rr_base = cfg->enigma_ring_present ? cfg->enigma_ring[2] : 0;   // pinned fast ring
    int n_units = n_orders * n_rr;

    int nthreads = cfg->n_threads > 0 ? cfg->n_threads : 1;
    if (nthreads > n_units) nthreads = n_units;
    if (nthreads < 1) nthreads = 1;

    printf("\nenigma: Bombe (crib) attack, %d ciphertext letters, %d menu edges, "
           "test letter %c, %d wheel order(s) x %d fast-ring x 26^%d start, %d thread(s)\n",
           cipher_len, ne, 'A' + test, n_orders, n_rr, tmpl[0].n_wheels, nthreads);

    // Split the (order x fast-ring) units across workers.
    BombeWork work[64];
    pthread_t th[64];
    int per = (n_units + nthreads - 1) / nthreads;
    int nspawn = 0;
    for (int t = 0; t < nthreads; t++) {
        int lo = t * per, hi = lo + per;
        if (lo >= n_units) break;
        if (hi > n_units) hi = n_units;
        work[nspawn] = (BombeWork){ .cfg = cfg, .tmpl = tmpl, .edge = edge, .want = want,
            .ne = ne, .maxpos = maxpos, .test = test, .cipher_len = cipher_len,
            .n_rr = n_rr, .rr_base = rr_base, .u_lo = lo, .u_hi = hi, .ntop = 0 };
        nspawn++;
    }
    if (nspawn <= 1) {
        bombe_worker(&work[0]);
    } else {
        for (int t = 0; t < nspawn; t++) pthread_create(&th[t], NULL, bombe_worker, &work[t]);
        for (int t = 0; t < nspawn; t++) pthread_join(th[t], NULL);
    }

    // Merge the workers' stops.
    BombeStop top[BOMBE_TOP_M];
    int ntop = 0;
    for (int t = 0; t < nspawn; t++)
        for (int i = 0; i < work[t].ntop; i++) bombe_keep(top, &ntop, &work[t].top[i], BOMBE_TOP_M);

    if (ntop == 0) {
        printf("\n\nEnigma Bombe: no consistent stop found (crib may be wrong or too short).\n\n");
        if (result) result->solved = false;
        return false;
    }

    // Sort stops by crib match (descending).
    for (int i = 0; i < ntop; i++)
        for (int j = i + 1; j < ntop; j++)
            if (top[j].matches > top[i].matches) { BombeStop t = top[i]; top[i] = top[j]; top[j] = t; }

    if (cfg->verbose)
        for (int i = 0; i < ntop; i++) {
            int nw = top[i].base.n_wheels;
            printf("  stop %d: %d/%d crib edges, rotors", i, top[i].matches, ne);
            for (int r = 0; r < nw; r++) printf(" %s", enigma_rotor_name(top[i].base.rotor[r]));
            printf(", pos %c%c%c\n", 'A' + top[i].base.pos[nw - 3],
                   'A' + top[i].base.pos[nw - 2], 'A' + top[i].base.pos[nw - 1]);
        }

    // Complete the best stops: ring refine + plugboard (seeded by the Bombe steckers) + report.
    EnigmaKey bases[BOMBE_TOP_M];
    for (int i = 0; i < ntop; i++) bases[i] = top[i].base;
    enigma_attack_from_bases(cfg, shared, cipher_indices, cipher_len, /*cribtext*/ NULL,
        crib_indices, crib_positions, n_cribs, bases, ntop, maxplugs, result);
    return result ? result->solved : true;
}

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
// across -nthreads over the order x fast-ring work. The positioned-crib path keeps M4 to a
// pinned wheel order; the CRIB DRAG below (-bombe -cribdrag W1|W2|..., no -crib) runs blind
// M4 (336 naval orders x 26 Greek positions x 26 fast rings) with Turing's loop test.

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

// =====================================================================
//  Bombe CRIB DRAG (-bombe -cribdrag W1|W2|... with no positioned -crib)
// =====================================================================
//
// The Bletchley crib drag: a probable word is slid along the ciphertext; Enigma never
// enciphers a letter to itself, so every offset where a crib letter coincides with its
// ciphertext letter is rejected outright, and each surviving (word, offset) becomes a Bombe
// menu. All menus are tested in ONE pass over the key space: per location the crib-span
// scramblers are gathered once (row copies from the per-unit table) and every menu runs its
// 26 test-letter hypotheses through the diagonal-board propagation against them. A menu's
// best stops are then completed QUIETLY (enigma_complete_base: ring refine + seeded reswap
// plugboard climb) and all completions are ranked by the full-decrypt n-gram score into a
// leaderboard (a wrong crib/offset yields a flood of stops whose completions score like
// noise; a right one stands out and reads). The global best is finally handed to the ordinary
// reported completion path so the ">>>" summary / result are the usual ones.
//
// M4: the Greek wheel never steps, so its position is folded into the unit (order x Greek
// position x fast ring): the per-unit table is the M3-sized 26^3 x 26 over the three stepping
// wheels and each start sweeps 26^3 -- blind M4 costs 336 x 26 units per fast ring. Cost is
// dominated by per-(menu, location) propagation (~1 us), i.e. ~10 s per menu per Greek/
// reflector pair on 16 threads with the fast ring pinned (-ring A,A,A), 26x with it swept.

#define DRAG_MAX_MENUS 2048

typedef struct {
    int word, off, ne, test;
    BombeEdge edge[BOMBE_MAX_EDGES];
    int ncyc, cyc[BOMBE_MAX_EDGES];   // Turing loop: a closed walk of edge indices test->...->test
} DragMenu;

// Turing's loop test (the pre-diagonal-board Bombe logic). Composing the crib scramblers
// around a closed walk that starts and ends at the test letter yields a permutation P of
// which the true stecker(test) must be a FIXED POINT. A random permutation has ~1 fixed point
// on average, so evaluating P once (|walk| x 26 lookups) leaves ~1 hypothesis to push through
// the full diagonal-board propagation instead of 26 -- ~5-8x fewer propagations per
// location. A menu without a loop has no such filter (and floods anyway), so it is skipped.
// The walk is a BFS tree path test->u, a non-tree edge (u,v), then the tree path v->test;
// shared prefixes compose to identity (S.S = id), leaving the underlying cycle's constraint.
static void drag_find_loop(DragMenu *mn) {
    int parent_edge[26], parent[26], dist[26], q[26], qh = 0, qt = 0;
    for (int i = 0; i < 26; i++) { parent_edge[i] = -1; parent[i] = -1; dist[i] = -1; }
    mn->ncyc = 0;
    dist[mn->test] = 0; q[qt++] = mn->test;
    while (qh < qt) {
        int u = q[qh++];
        for (int j = 0; j < mn->ne; j++) {
            int a = mn->edge[j].p, b = mn->edge[j].c;
            if (a != u && b != u) continue;
            int v = (a == u) ? b : a;
            if (dist[v] < 0) { dist[v] = dist[u] + 1; parent[v] = u; parent_edge[v] = j; q[qt++] = v; }
        }
    }
    // Shortest closed walk: the non-tree edge (u,v) minimising dist[u]+dist[v]+1.
    int best = -1, bestlen = 1 << 30;
    for (int j = 0; j < mn->ne; j++) {
        int a = mn->edge[j].p, b = mn->edge[j].c;
        if (dist[a] < 0 || dist[b] < 0) continue;
        if (parent_edge[a] == j || parent_edge[b] == j) continue;     // tree edge
        int len = dist[a] + dist[b] + 1;
        if (len < bestlen) { bestlen = len; best = j; }
    }
    if (best < 0) return;                                             // no loop
    int a = mn->edge[best].p, b = mn->edge[best].c;
    // Walk: test -> a (tree path forward), edge best, b -> test (tree path backward).
    int path_a[26], na = 0;
    for (int x = a; x != mn->test; x = parent[x]) path_a[na++] = parent_edge[x];
    for (int i = na - 1; i >= 0; i--) mn->cyc[mn->ncyc++] = path_a[i];
    mn->cyc[mn->ncyc++] = best;
    for (int x = b; x != mn->test; x = parent[x]) mn->cyc[mn->ncyc++] = parent_edge[x];
}

typedef struct {
    const ColossusConfig *cfg;
    const EnigmaKey *tmpl;
    const DragMenu *menu; int nmenus;
    const int *cipher; int cipher_len, maxpos;
    int n_g, n_rr, rr_base;      // Greek positions per order (26 for M4, else 1); fast rings
    int u_lo, u_hi;
    BombeStop *top;              // [nmenus][BOMBE_TOP_M]
    int *ntop;                   // [nmenus]
    long *nstops;                // [nmenus] consistent stops seen (flood diagnostic)
} DragWork;

// bombe_propagate with the scramblers indexed by crib POSITION (shared across menus).
static bool drag_propagate(int steck[26], const BombeEdge *e, int ne, int (*scr)[26]) {
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int j = 0; j < ne; j++) {
            int a = e[j].p, b = e[j].c, va = steck[a], vb = steck[b];
            const int *S = scr[e[j].pos];
            if (va != -1 && vb == -1)      { if (!steck_set(steck, b, S[va])) return false; changed = 1; }
            else if (vb != -1 && va == -1) { if (!steck_set(steck, a, S[vb])) return false; changed = 1; }
            else if (va != -1 && vb != -1) { if (S[va] != vb) return false; }
        }
    }
    return true;
}

static void *drag_worker(void *arg) {
    DragWork *w = (DragWork *) arg;
    const ColossusConfig *cfg = w->cfg;
    int *S = (int *) malloc((size_t) 26 * 26 * 26 * 26 * sizeof(int));
    int (*scr)[26] = (int (*)[26]) malloc(sizeof(int[26]) * (size_t) (w->maxpos + 1));
    int per_order = w->n_g * w->n_rr;

    for (int u = w->u_lo; u < w->u_hi; u++) {
        int o   = u / per_order, rem = u % per_order;
        int g   = rem / w->n_rr;
        int rr  = w->rr_base + (rem % w->n_rr);
        EnigmaKey kbase = w->tmpl[o];
        int nw = kbase.n_wheels, off = nw - 3;
        for (int i = 0; i < nw; i++) kbase.ring[i] = 0;
        if (cfg->enigma_ring_present)
            for (int i = 0; i < 3; i++) kbase.ring[off + i] = cfg->enigma_ring[i];
        kbase.ring[nw - 1] = rr;
        enigma_plug_identity(kbase.plug);
        if (nw == 4) kbase.pos[0] = g;

        // Per-unit scrambler table over the three stepping wheels (Greek fixed at g).
        EnigmaKey kt = kbase;
        for (int a = 0; a < 26; a++)
            for (int b = 0; b < 26; b++)
                for (int c = 0; c < 26; c++) {
                    kt.pos[off] = a; kt.pos[off + 1] = b; kt.pos[off + 2] = c;
                    int base = ((a * 26 + b) * 26 + c) * 26;
                    for (int ch = 0; ch < 26; ch++) S[base + ch] = enigma_encipher_letter(&kt, ch);
                }

        for (long pp = 0; pp < 26L * 26 * 26; pp++) {
            EnigmaKey k = kbase;
            k.pos[off] = (int) (pp / 676); k.pos[off + 1] = (int) ((pp / 26) % 26); k.pos[off + 2] = (int) (pp % 26);
            if (cfg->enigma_pos_present)
                for (int i = 0; i < 3; i++) k.pos[off + i] = cfg->enigma_pos[i];

            // Gather every position's scrambler once (shared by all menus).
            EnigmaKey kk = k;
            for (int s = 0; s <= w->maxpos; s++) {
                enigma_step(&kk);
                memcpy(scr[s], &S[((kk.pos[off] * 26 + kk.pos[off + 1]) * 26 + kk.pos[off + 2]) * 26],
                       26 * sizeof(int));
            }

            for (int m = 0; m < w->nmenus; m++) {
                const DragMenu *mn = &w->menu[m];
                int best_match = -1, best_steck[26];
                // Turing loop test: keep only the fixed points of the loop composition.
                unsigned cand = 0;
                for (int h = 0; h < 26; h++) {
                    int y = h;
                    for (int t = 0; t < mn->ncyc; t++) y = scr[mn->edge[mn->cyc[t]].pos][y];
                    if (y == h) cand |= 1u << h;
                }
                if (!cand) continue;
                for (int h = 0; h < 26; h++) {
                    if (!(cand & (1u << h))) continue;
                    int steck[26];
                    for (int i = 0; i < 26; i++) steck[i] = -1;
                    if (!steck_set(steck, mn->test, h)) continue;
                    if (!drag_propagate(steck, mn->edge, mn->ne, scr)) continue;
                    int plug[26];
                    for (int i = 0; i < 26; i++) plug[i] = (steck[i] >= 0) ? steck[i] : i;
                    int mt = 0;
                    for (int j = 0; j < mn->ne; j++)
                        if (plug[scr[mn->edge[j].pos][plug[mn->edge[j].c]]] == mn->edge[j].p) mt++;
                    if (mt > best_match) { best_match = mt; memcpy(best_steck, steck, sizeof(steck)); }
                }
                if (best_match > 0) {
                    w->nstops[m]++;
                    BombeStop s;
                    s.base = k;
                    for (int i = 0; i < 26; i++) s.base.plug[i] = (best_steck[i] >= 0) ? best_steck[i] : i;
                    s.matches = best_match;
                    bombe_keep(&w->top[(size_t) m * BOMBE_TOP_M], &w->ntop[m], &s, BOMBE_TOP_M);
                }
            }
            if (cfg->enigma_pos_present) break;
        }
    }
    free(scr);
    free(S);
    return NULL;
}

typedef struct { double score; int menu; EnigmaKey key; } DragHit;

bool solve_enigma_bombe_drag(ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len, SolveResult *result) {

    enigma_init();
    const CribDrag *cd = &cfg->cribdrag;

    // Build the menus: every (word, offset) the no-self-encipherment test allows.
    static DragMenu menu[DRAG_MAX_MENUS];
    int nm = 0, maxpos = 0, n_offsets = 0, n_clash = 0, n_noloop = 0;
    for (int wi = 0; wi < cd->nwords; wi++) {
        int L = cd->wordlen[wi];
        if (L < 3 || L > cipher_len) continue;
        int off_hi = cipher_len - L;
        if (cfg->cribdrag_max_offset >= 0 && off_hi > cfg->cribdrag_max_offset) off_hi = cfg->cribdrag_max_offset;
        for (int off = 0; off <= off_hi; off++) {
            n_offsets++;
            bool clash = false;
            for (int i = 0; i < L; i++) if (cd->words[wi][i] == cipher_indices[off + i]) { clash = true; break; }
            if (clash) { n_clash++; continue; }
            if (nm >= DRAG_MAX_MENUS) break;
            DragMenu *mn = &menu[nm];
            int deg[26] = {0};
            mn->word = wi; mn->off = off; mn->ne = 0;
            for (int i = 0; i < L && mn->ne < BOMBE_MAX_EDGES; i++) {
                BombeEdge *e = &mn->edge[mn->ne++];
                e->p = cd->words[wi][i]; e->c = cipher_indices[off + i]; e->pos = off + i;
                deg[e->p]++; deg[e->c]++;
                if (e->pos > maxpos) maxpos = e->pos;
            }
            mn->test = 0;
            for (int i = 1; i < 26; i++) if (deg[i] > deg[mn->test]) mn->test = i;
            drag_find_loop(mn);
            if (mn->ncyc == 0) { n_noloop++; continue; }                 // no Turing loop: skip
            nm++;
        }
    }
    printf("\nenigma: Bombe CRIB DRAG, %d ciphertext letters, %d word(s), %d offsets tested, "
           "%d rejected by self-encipherment, %d skipped (no Turing loop), %d menus\n",
           cipher_len, cd->nwords, n_offsets, n_clash, n_noloop, nm);
    if (nm == 0) {
        printf("\n\nERROR: crib drag produced no admissible (word, offset) menu.\n\n");
        return false;
    }

    int maxplugs = (cfg->enigma_maxplugs > 0) ? cfg->enigma_maxplugs : 10;
    if (maxplugs > 13) maxplugs = 13;

    static EnigmaKey tmpl[512];
    int n_orders = enigma_enumerate_wheel_orders(cfg, tmpl, 512);
    int n_g  = (tmpl[0].n_wheels == 4) ? 26 : 1;
    int n_rr = cfg->enigma_ring_present ? 1 : 26;
    int rr_base = cfg->enigma_ring_present ? cfg->enigma_ring[2] : 0;
    int n_units = n_orders * n_g * n_rr;
    int nthreads = cfg->n_threads > 0 ? cfg->n_threads : 1;
    if (nthreads > n_units) nthreads = n_units;
    if (nthreads > 64) nthreads = 64;
    if (nthreads < 1) nthreads = 1;
    printf("enigma: %d wheel order(s) x %d Greek pos x %d fast-ring x 26^3 start = %d units, %d thread(s)\n",
           n_orders, n_g, n_rr, n_units, nthreads);

    DragWork work[64];
    pthread_t th[64];
    int per = (n_units + nthreads - 1) / nthreads, nspawn = 0;
    for (int t = 0; t < nthreads; t++) {
        int lo = t * per, hi = lo + per;
        if (lo >= n_units) break;
        if (hi > n_units) hi = n_units;
        work[nspawn] = (DragWork){ .cfg = cfg, .tmpl = tmpl, .menu = menu, .nmenus = nm,
            .cipher = cipher_indices, .cipher_len = cipher_len, .maxpos = maxpos,
            .n_g = n_g, .n_rr = n_rr, .rr_base = rr_base, .u_lo = lo, .u_hi = hi,
            .top = (BombeStop *) calloc((size_t) nm * BOMBE_TOP_M, sizeof(BombeStop)),
            .ntop = (int *) calloc((size_t) nm, sizeof(int)),
            .nstops = (long *) calloc((size_t) nm, sizeof(long)) };
        nspawn++;
    }
    if (nspawn <= 1) drag_worker(&work[0]);
    else {
        for (int t = 0; t < nspawn; t++) pthread_create(&th[t], NULL, drag_worker, &work[t]);
        for (int t = 0; t < nspawn; t++) pthread_join(th[t], NULL);
    }

    // Merge per-menu stops; complete each quietly; rank.
    DragHit *hits = (DragHit *) malloc((size_t) nm * BOMBE_TOP_M * sizeof(DragHit));
    long *nstops = (long *) calloc((size_t) nm, sizeof(long));
    int nh = 0;
    for (int m = 0; m < nm; m++) {
        BombeStop top[BOMBE_TOP_M]; int ntop = 0;
        for (int t = 0; t < nspawn; t++) {
            nstops[m] += work[t].nstops[m];
            for (int i = 0; i < work[t].ntop[m]; i++)
                bombe_keep(top, &ntop, &work[t].top[(size_t) m * BOMBE_TOP_M + i], BOMBE_TOP_M);
        }
        for (int i = 0; i < ntop; i++) {
            DragHit *h = &hits[nh];
            h->menu = m;
            h->score = enigma_complete_base(cfg, shared, cipher_indices, cipher_len, &top[i].base, maxplugs, &h->key);
            nh++;
        }
    }
    for (int t = 0; t < nspawn; t++) { free(work[t].top); free(work[t].ntop); free(work[t].nstops); }
    for (int i = 0; i < nh; i++)
        for (int j = i + 1; j < nh; j++)
            if (hits[j].score > hits[i].score) { DragHit t = hits[i]; hits[i] = hits[j]; hits[j] = t; }

    if (nh == 0) {
        printf("\nenigma: crib drag found NO consistent stop for any (word, offset) -- every crib is "
               "inconsistent with this key space (wrong words, wrong -greek/-reflector, or wrong model).\n\n");
        for (int t = 0; t < nspawn; t++) (void) t;
        free(hits); free(nstops);
        if (result) result->solved = false;
        return false;
    }

    // Leaderboard.
    int nshow = nh < 30 ? nh : 30;
    printf("\nenigma: crib-drag leaderboard (%d completed stops; score = full-decrypt n-gram)\n", nh);
    static int dec[MAX_CIPHER_LENGTH];
    for (int i = 0; i < nshow; i++) {
        const DragMenu *mn = &menu[hits[i].menu];
        char word[MAX_CRIBDRAG_LEN + 1];
        for (int c = 0; c < cd->wordlen[mn->word]; c++) word[c] = 'A' + cd->words[mn->word][c];
        word[cd->wordlen[mn->word]] = '\0';
        EnigmaKey *k = &hits[i].key;
        int nw = k->n_wheels;
        char ro[64] = {0}, plugs[128];
        for (int j = 0; j < nw; j++) { strcat(ro, enigma_rotor_name(k->rotor[j])); if (j + 1 < nw) strcat(ro, " "); }
        enigma_format_plugs(k->plug, plugs);
        enigma_encrypt(cipher_indices, cipher_len, k, dec);
        printf("  %2d. %.4f  %-24s @%-3d stops=%-7ld %s %s rings=", i + 1, hits[i].score, word, mn->off,
               nstops[hits[i].menu], enigma_reflector_name(k->reflector), ro);
        for (int j = 0; j < nw; j++) printf("%c", 'A' + k->ring[j]);
        printf(" pos=");
        for (int j = 0; j < nw; j++) printf("%c", 'A' + k->pos[j]);
        printf(" plugs=%s\n      ", plugs);
        for (int j = 0; j < cipher_len; j++) printf("%c", 'A' + dec[j]);
        printf("\n");
    }

    // Report the global best through the ordinary completion path (positioned crib = its menu).
    const DragMenu *bm = &menu[hits[0].menu];
    int crib_idx[BOMBE_MAX_EDGES], crib_pos[BOMBE_MAX_EDGES];
    for (int j = 0; j < bm->ne; j++) { crib_idx[j] = bm->edge[j].p; crib_pos[j] = bm->edge[j].pos; }
    EnigmaKey bases[BOMBE_TOP_M]; int nb = 0;
    for (int i = 0; i < nh && nb < BOMBE_TOP_M; i++)
        if (hits[i].menu == hits[0].menu) bases[nb++] = hits[i].key;
    enigma_attack_from_bases(cfg, shared, cipher_indices, cipher_len, NULL,
        crib_idx, crib_pos, bm->ne, bases, nb, maxplugs, result);
    free(hits); free(nstops);
    return result ? result->solved : true;
}

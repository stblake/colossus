// =====================================================================
//  Enigma solver (TYPE enigma)
// =====================================================================
//
// Ciphertext-only attack after Gillogly (1995) / Ostwald-Weierud (2017) / Mike Pound.
// The plugboard is a FIXED substitution that barely perturbs letter-frequency statistics,
// so the rotor config shows through the Index of Coincidence even with an empty plugboard.
// Three phases:
//
//   Phase 1  Wheel order + start positions.  For each ordered triple of rotors (default
//            pool I..V, 60 orders) and each of 26^k start positions (rings = AAA, empty
//            plugboard), decrypt and score IoC. Keep the best start per order and the top-K
//            orders.  This is the ~1M-decrypt "Bombe workload".
//   Phase 2  Ring settings.  For each surviving candidate, spin the fast (right) ring 0..25
//            moving that rotor's start in tandem (IoC), then the middle ring likewise. The
//            left ring is irrelevant. This converts Gillogly's rings-AAA hit into the true
//            ring setting (his "rings 1 23 4 / msg 2 7 9").
//   Phase 3  Plugboard.  A cheap greedy warm start ranks the candidates; the best few then get
//            a full deterministic reswap climb (enigma_plugboard_climb). The shared engine then
//            anneals the plugboard from that seed (one config per candidate) under the n-gram
//            (+ crib) fitness and keeps the global best.
//
// The plugboard phase runs through the standard CipherModel/run_solver engine, so -logprob,
// cribs (real positional cribs -- Enigma is length-preserving), -method, and -nthreads all
// apply. Phases 1-2 are deterministic IoC pre-passes.
//
// A known-key path (enough of the key pinned via -rotors/-ring/-startpos/-plugboard) skips
// the search and just decrypts -- the drag/verification capability. The Bombe (crib) attack
// lives in enigma_bombe.c and is selected by -bombe.

#include <string.h>
#include <pthread.h>
#include "enigma_solver.h"
#include "engine.h"
#include "scoring.h"

#define ENIGMA_MAX_CAND   16      // top-K wheel-order candidates carried to phase 3
#define ENIGMA_PHASE1_TOPP 32     // per-worker top positions kept in a phase-1 IoC pass
#define ENIGMA_PHASE1_ORDERS 12   // distinct orders from the coarse pass fed to the fast-ring pass
#define ENIGMA_POOL_SIZE  5       // default blind rotor pool = I..V (Services Enigma I)
static const int ENIGMA_POOL[ENIGMA_POOL_SIZE] = {
    ENIGMA_I, ENIGMA_II, ENIGMA_III, ENIGMA_IV, ENIGMA_V,
};

typedef struct {
    EnigmaKey key;      // fixed order/ring/pos; .plug is the greedy warm-start plugboard
    double    ioc;      // phase-2 IoC (diagnostic)
} EnigmaCand;

typedef struct {
    int        n_cand;
    EnigmaCand cand[ENIGMA_MAX_CAND];
    int        maxplugs;
} EnigmaScratch;

// ------------------------------------------------------------------ plugboard climb
//
// The plugboard is the ONLY smooth part of the Enigma keyspace (the rotor order/position/ring
// are a non-smooth needle that phase 1 EXHAUSTS by IoC; annealing them jointly with the plugs
// cannot climb toward the true rotors -- there is no gradient). Once the rotors are fixed the
// plugboard is a substitution and hill-climbs cleanly. This is the Pound / Ostwald-Weierud
// plugboard climb: RESWAP moves (each pass sweeps all C(26,2) "set plug (i,j)" moves + the 26
// removals and applies the single best-improving one) with a few RESTARTS to escape local
// optima. Strictly stronger than greedy-add-only, which can never revisit a wrong early plug
// and so leaves "one stecker short" near-solutions.

#define ENIGMA_PLUG_RESTARTS 3
#define ENIGMA_PLUG_CLIMB_TOP 6   // run the full reswap climb on this many best-ranked candidates

// Apply the "set plug (i,j)" move to p IN PLACE: disconnect i and j from their partners, then
// connect i<->j if there is room. i==j is a pure DISCONNECT (removal). Over all i<=j the move
// set is add + swap + remove.
static void enigma_plug_move(int p[26], int i, int j, int maxplugs) {
    int pi = p[i], pj = p[j];
    p[pi] = pi; p[i] = i; p[pj] = pj; p[j] = j;
    if (i != j) {
        int pairs = 0;
        for (int k = 0; k < 26; k++) if (p[k] > k) pairs++;
        if (pairs < maxplugs) { p[i] = j; p[j] = i; }
    }
}

static double enigma_plug_score(const EnigmaKey *key, const int plug[26], int cipher[], int len,
                                float *ngram_data, int ngram_size, int *dec) {
    EnigmaKey k = *key;
    memcpy(k.plug, plug, sizeof(int) * 26);
    enigma_encrypt(cipher, len, &k, dec);
    return ngram_score(dec, len, ngram_data, ngram_size);
}

// Fast greedy-ADD plugboard estimate (capped), for cheaply RANKING many rotor candidates
// before the full climb is spent on the best few. Returns the n-gram score.
static double enigma_quick_plug(const EnigmaKey *key, int cipher[], int cipher_len,
    float *ngram_data, int ngram_size, int maxplugs, int out_plug[26]) {
    static _Thread_local int dec[MAX_CIPHER_LENGTH];
    int plug[26]; enigma_plug_identity(plug);
    double best = enigma_plug_score(key, plug, cipher, cipher_len, ngram_data, ngram_size, dec);
    for (int n = 0; n < maxplugs; n++) {
        int ba = -1, bb = -1; double bgain = 1e-9;
        for (int a = 0; a < 26; a++) { if (plug[a] != a) continue;
            for (int b = a + 1; b < 26; b++) { if (plug[b] != b) continue;
                plug[a] = b; plug[b] = a;
                double s = enigma_plug_score(key, plug, cipher, cipher_len, ngram_data, ngram_size, dec);
                plug[a] = a; plug[b] = b;
                if (s - best > bgain) { bgain = s - best; ba = a; bb = b; }
            } }
        if (ba < 0) break;
        plug[ba] = bb; plug[bb] = ba; best += bgain;
    }
    memcpy(out_plug, plug, sizeof(int) * 26);
    return best;
}

// Full reswap plugboard climb (see above). `seed_plug` (or NULL) is restart 0's starting board
// -- pass the tier-A greedy result or a Bombe stecker seed so the climb refines it.
double enigma_plugboard_climb(const EnigmaKey *key, int cipher[], int cipher_len,
    float *ngram_data, int ngram_size, int maxplugs, const int *seed_plug, int out_plug[26]) {
    static _Thread_local int dec[MAX_CIPHER_LENGTH];
    int best_plug[26];
    double best_overall = -1e18;

    for (int r = 0; r < ENIGMA_PLUG_RESTARTS; r++) {
        int plug[26];
        if (r == 0 && seed_plug) memcpy(plug, seed_plug, sizeof(int) * 26);
        else enigma_plug_identity(plug);
        if (r > 0) {                                    // diversify: seed with r random plugs
            int avail[26]; for (int i = 0; i < 26; i++) avail[i] = i; int na = 26;
            int seedn = (r < maxplugs) ? r : maxplugs;
            for (int p = 0; p < seedn && na >= 2; p++) {
                int a = rand_int(0, na); int x = avail[a]; avail[a] = avail[--na];
                int b = rand_int(0, na); int y = avail[b]; avail[b] = avail[--na];
                plug[x] = y; plug[y] = x;
            }
        }
        double cur = enigma_plug_score(key, plug, cipher, cipher_len, ngram_data, ngram_size, dec);
        for (;;) {                                      // reswap hill-climb to a local optimum
            double bestgain = 1e-9; int bi = -1, bj = -1;
            for (int i = 0; i < 26; i++)
                for (int j = i; j < 26; j++) {
                    int p[26]; memcpy(p, plug, sizeof(p));
                    enigma_plug_move(p, i, j, maxplugs);
                    double s = enigma_plug_score(key, p, cipher, cipher_len, ngram_data, ngram_size, dec);
                    if (s - cur > bestgain) { bestgain = s - cur; bi = i; bj = j; }
                }
            if (bi < 0) break;
            enigma_plug_move(plug, bi, bj, maxplugs);
            cur += bestgain;
        }
        if (cur > best_overall) { best_overall = cur; memcpy(best_plug, plug, sizeof(best_plug)); }
    }
    memcpy(out_plug, best_plug, sizeof(int) * 26);
    return best_overall;
}

// ------------------------------------------------------------------ IoC search helpers

// Decrypt with an empty plugboard and return the IoC of the result.
static double enigma_ioc_of(const EnigmaKey *key, int cipher[], int cipher_len) {
    static _Thread_local int dec[MAX_CIPHER_LENGTH];
    EnigmaKey k = *key;
    enigma_plug_identity(k.plug);
    enigma_encrypt(cipher, cipher_len, &k, dec);
    return index_of_coincidence(dec, cipher_len);
}

// Phase-1 coarse IoC search worker. Over its (wheel order x fast-ring) unit range it sweeps
// all 26^n_wheels start positions (middle/left ring 0, fast ring = the unit's rr) and keeps
// the top-P configs by IoC. Searching the FAST ring here (not just position) makes the
// decrypt correct all the way to the much-later MIDDLE-rotor turnover -- not just the first
// right-rotor turnover -- so the true config's IoC stands out even on short messages (the
// 60 x 26^4 "Bombe workload" of Gillogly / Ostwald-Weierud). The caller then middle-ring-
// refines the merged top-P and ranks by that.
typedef struct {
    const ColossusConfig *cfg;
    const EnigmaKey *tmpl;
    int *cipher, cipher_len;
    int  n_rr, rr_base, u_lo, u_hi;    // (order x fast-ring) unit range
    EnigmaKey posP[ENIGMA_PHASE1_TOPP];
    double    iocP[ENIGMA_PHASE1_TOPP];
    int nP;
} Phase1Work;

static void *enigma_phase1_worker(void *arg) {
    Phase1Work *w = (Phase1Work *) arg;
    const ColossusConfig *cfg = w->cfg;
    w->nP = 0;
    for (int u = w->u_lo; u < w->u_hi; u++) {
        int o = u / w->n_rr, rr = w->rr_base + (u % w->n_rr);
        EnigmaKey kb = w->tmpl[o];
        int nw = kb.n_wheels;
        for (int i = 0; i < nw; i++) kb.ring[i] = 0;
        if (cfg->enigma_ring_present) {
            int off = nw - 3;
            for (int i = 0; i < 3; i++) kb.ring[off + i] = cfg->enigma_ring[i];
        }
        kb.ring[nw - 1] = rr;
        long total = 1;
        for (int i = 0; i < nw; i++) total *= 26;
        for (long pp = 0; pp < total; pp++) {
            EnigmaKey k = kb;
            long q = pp;
            for (int i = nw - 1; i >= 0; i--) { k.pos[i] = (int)(q % 26); q /= 26; }
            if (cfg->enigma_pos_present) {
                int off = nw - 3;
                if (nw == 4) k.pos[0] = 0;
                for (int i = 0; i < 3; i++) k.pos[off + i] = cfg->enigma_pos[i];
            }
            double ioc = enigma_ioc_of(&k, w->cipher, w->cipher_len);
            if (w->nP < ENIGMA_PHASE1_TOPP) { w->posP[w->nP] = k; w->iocP[w->nP] = ioc; w->nP++; }
            else {
                int worst = 0;
                for (int j = 1; j < ENIGMA_PHASE1_TOPP; j++) if (w->iocP[j] < w->iocP[worst]) worst = j;
                if (ioc > w->iocP[worst]) { w->posP[worst] = k; w->iocP[worst] = ioc; }
            }
            if (cfg->enigma_pos_present) break;
        }
    }
    return NULL;
}

// Phase 2: refine the ring of stepping wheel `w` (moving its start in tandem so the WIRING
// offset (pos - ring) stays fixed and only the turnover timing sweeps) by IoC. Works for any
// starting ring (base ring 0 for the ciphertext-only path; the Bombe's recovered fast ring
// otherwise) -- for a ring-0 base it is identical to the naive pos0+d sweep.
static void enigma_refine_ring(EnigmaKey *key, int w, int cipher[], int cipher_len) {
    int offset = ((key->pos[w] - key->ring[w]) % 26 + 26) % 26;
    double best = -1.0; int best_d = 0;
    for (int d = 0; d < 26; d++) {
        key->ring[w] = d;
        key->pos[w] = (offset + d) % 26;
        double ioc = enigma_ioc_of(key, cipher, cipher_len);
        if (ioc > best) { best = ioc; best_d = d; }
    }
    key->ring[w] = best_d;
    key->pos[w] = (offset + best_d) % 26;
}

// Threaded phase-1 IoC search over tmpl[0..n_orders). n_rr=1 => a rings-coarse pass (cheap,
// for ranking wheel orders); n_rr=26 => the full fast-ring search (robust position recovery).
// Merges every worker's top positions into out[]/outioc[] (top-`cap`), returns the count.
static int enigma_phase1_search(const ColossusConfig *cfg, const EnigmaKey *tmpl, int n_orders,
    int n_rr, int rr_base, int *cipher, int cipher_len,
    EnigmaKey out[], double outioc[], int cap) {
    int n_units = n_orders * n_rr;
    int nthreads = cfg->n_threads > 0 ? cfg->n_threads : 1;
    if (nthreads > n_units) nthreads = n_units;
    if (nthreads < 1) nthreads = 1;

    static Phase1Work work[64];   // static: each holds a top-P buffer; keep off the stack
    pthread_t th[64];
    int per = (n_units + nthreads - 1) / nthreads, nspawn = 0;
    for (int t = 0; t < nthreads; t++) {
        int lo = t * per, hi = lo + per;
        if (lo >= n_units) break;
        if (hi > n_units) hi = n_units;
        work[nspawn].cfg = cfg; work[nspawn].tmpl = tmpl; work[nspawn].cipher = cipher;
        work[nspawn].cipher_len = cipher_len; work[nspawn].n_rr = n_rr; work[nspawn].rr_base = rr_base;
        work[nspawn].u_lo = lo; work[nspawn].u_hi = hi; work[nspawn].nP = 0;
        nspawn++;
    }
    if (nspawn <= 1) {
        enigma_phase1_worker(&work[0]);
    } else {
        for (int t = 0; t < nspawn; t++) pthread_create(&th[t], NULL, enigma_phase1_worker, &work[t]);
        for (int t = 0; t < nspawn; t++) pthread_join(th[t], NULL);
    }

    int n = 0;
    for (int t = 0; t < nspawn; t++)
        for (int i = 0; i < work[t].nP; i++) {
            if (n < cap) { out[n] = work[t].posP[i]; outioc[n] = work[t].iocP[i]; n++; }
            else {
                int worst = 0;
                for (int j = 1; j < cap; j++) if (outioc[j] < outioc[worst]) worst = j;
                if (work[t].iocP[i] > outioc[worst]) { out[worst] = work[t].posP[i]; outioc[worst] = work[t].iocP[i]; }
            }
        }
    return n;
}

// Distinct wheel orders (by rotor[]) among the phase-1a top positions -> tmpl2[], for the
// fast-ring pass. Returns the count (<= cap).
static int enigma_distinct_orders(const EnigmaKey posA[], int nA, EnigmaKey tmpl2[], int cap) {
    int n = 0;
    for (int i = 0; i < nA && n < cap; i++) {
        int dup = 0;
        for (int j = 0; j < n; j++) {
            int same = 1;
            for (int r = 0; r < posA[i].n_wheels; r++)
                if (tmpl2[j].rotor[r] != posA[i].rotor[r]) { same = 0; break; }
            if (same && tmpl2[j].reflector == posA[i].reflector) { dup = 1; break; }
        }
        if (dup) continue;
        EnigmaKey t = posA[i];
        for (int w = 0; w < t.n_wheels; w++) { t.ring[w] = 0; t.pos[w] = 0; }
        enigma_plug_identity(t.plug);
        tmpl2[n++] = t;
    }
    return n;
}

// Insert (order, key) into the top-K candidate list kept sorted by IoC descending.
static void enigma_topk_insert(EnigmaScratch *s, const EnigmaKey *key, double ioc, int K) {
    int slot = s->n_cand;
    if (slot >= K) {
        // replace the weakest if this is better
        int worst = 0;
        for (int i = 1; i < s->n_cand; i++) if (s->cand[i].ioc < s->cand[worst].ioc) worst = i;
        if (ioc <= s->cand[worst].ioc) return;
        slot = worst;
    } else {
        s->n_cand++;
    }
    s->cand[slot].key = *key;
    s->cand[slot].ioc = ioc;
}

// ------------------------------------------------------------------ CipherModel (phase 3)

static int enigma_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    const EnigmaScratch *s = (const EnigmaScratch *) ctx->model_scratch;
    int n = s->n_cand; if (n > cap) n = cap;
    for (int i = 0; i < n; i++) {
        out[i].period = 0; out[i].j = 0; out[i].k = 0;
        out[i].aux[0] = i; out[i].aux[1] = 0;
    }
    return n;
}

static int enigma_key_len(const SolverCtx *ctx, const SolverConfig *cc) {
    (void) ctx; (void) cc;
    return 26;   // the climbed key is the 26-entry plugboard involution
}

static void enigma_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    const EnigmaScratch *s = (const EnigmaScratch *) ctx->model_scratch;
    memcpy(st->key, s->cand[cc->aux[0]].key.plug, sizeof(int) * 26);  // greedy warm start
    st->key_len = 26;
}

static void enigma_perturb(const SolverCtx *ctx, const SolverConfig *cc,
                           SolverState *st, bool *force_primary) {
    (void) cc; (void) force_primary;
    const EnigmaScratch *s = (const EnigmaScratch *) ctx->model_scratch;
    int *plug = st->key;
    double r = frand();
    if (r < 0.15) {                                     // remove a random plug
        int a = rand_int(0, 26);
        if (plug[a] != a) { int pa = plug[a]; plug[a] = a; plug[pa] = pa; }
        return;
    }
    int a = rand_int(0, 26);
    int b = rand_int(0, 25); if (b >= a) b++;           // b != a
    int pa = plug[a]; plug[a] = a; plug[pa] = pa;       // unplug a (+ partner)
    int pb = plug[b]; plug[b] = b; plug[pb] = pb;       // unplug b (+ partner)
    int pairs = 0;
    for (int i = 0; i < 26; i++) if (plug[i] > i) pairs++;
    if (pairs < s->maxplugs) { plug[a] = b; plug[b] = a; }  // plug a<->b if room
}

static void enigma_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    (void) cc;
    memcpy(dst->key, src->key, sizeof(int) * 26);
    dst->key_len = src->key_len;
}

static void enigma_decrypt_hook(const SolverCtx *ctx, const SolverConfig *cc,
                                SolverState *st, int *out, double *score_adjust) {
    const EnigmaScratch *s = (const EnigmaScratch *) ctx->model_scratch;
    EnigmaKey k = s->cand[cc->aux[0]].key;
    memcpy(k.plug, st->key, sizeof(int) * 26);
    enigma_encrypt(ctx->cipher, ctx->cipher_len, &k, out);
    *score_adjust = 0.0;
}

static void enigma_report_hook(const SolverCtx *ctx, const SolverConfig *cc,
                               const SolverState *st, double score, int *decrypted) {
    const EnigmaScratch *s = (const EnigmaScratch *) ctx->model_scratch;
    EnigmaKey k = s->cand[cc->aux[0]].key;
    memcpy(k.plug, st->key, sizeof(int) * 26);
    enigma_emit_report(ctx->cfg, ctx->shared, ctx->cipher, ctx->cipher_len, &k,
                       decrypted, ctx->cribtext, score, ctx->result);
}

static const CipherModel ENIGMA_MODEL = {
    .name = "enigma", .shape = SHAPE_ANNEAL, .needs_hist = false,
    .enumerate_configs = enigma_enumerate, .key_len = enigma_key_len,
    .seed = enigma_seed, .perturb = enigma_perturb, .copy_state = enigma_copy,
    .decrypt = enigma_decrypt_hook, .report = enigma_report_hook,
};

// ------------------------------------------------------------------ report

// Render the plugboard involution as "AB CD EF" pairs; returns chars written.
static void enigma_format_plugs(const int plug[26], char *buf) {
    char *p = buf; int first = 1;
    for (int i = 0; i < 26; i++) if (plug[i] > i) {
        if (!first) *p++ = ' ';
        *p++ = 'A' + i; *p++ = 'A' + plug[i];
        first = 0;
    }
    if (first) { strcpy(buf, "(none)"); return; }
    *p = '\0';
}

void enigma_emit_report(ColossusConfig *cfg, SharedData *shared,
    int cipher[], int cipher_len, const EnigmaKey *key, int decrypted[], char *cribtext,
    double score, SolveResult *result) {

    int n = key->n_wheels;
    char rotors[64] = {0}, rings[16] = {0}, pos[16] = {0}, plugs[128];
    for (int i = 0; i < n; i++) {
        strcat(rotors, enigma_rotor_name(key->rotor[i]));
        if (i + 1 < n) strcat(rotors, " ");
        rings[i] = 'A' + key->ring[i];
        pos[i]   = 'A' + key->pos[i];
    }
    enigma_format_plugs(key->plug, plugs);

    int n_words = 0;
    char pt[MAX_CIPHER_LENGTH];
    for (int i = 0; i < cipher_len; i++) pt[i] = index_to_char(decrypted[i]);
    pt[cipher_len] = '\0';
    if (cfg->dictionary_present && shared->dict != NULL)
        n_words = find_dictionary_words(pt, shared->dict, shared->n_dict_words,
                                        shared->max_dict_word_len);

    printf("\nResult Score: %.4f | Words: %d | M%d | reflector %s | rotors %s | rings %s | pos %s | plugs %s\n",
        score, n_words, n, enigma_reflector_name(key->reflector), rotors, rings, pos, plugs);
    print_cipher(cipher, cipher_len, NULL);
    printf("\n");
    print_text(decrypted, cipher_len);
    printf("\n");
    print_solution_check(decrypted, cipher_len);
    if (cribtext) printf("%s\n", cribtext);

    if (result) {
        result->solved = true;
        result->cipher_type = cfg->cipher_type;
        result->score = score;
        result->n_words = n_words;
        vec_copy(decrypted, result->decrypted, cipher_len);
        result->decrypted_len = cipher_len;
    }

    printf(">>> %.4f, %d, refl=%s, rotors=%s, rings=%s, pos=%s, plugs=%s, ",
        score, cfg->cipher_type, enigma_reflector_name(key->reflector),
        rotors, rings, pos, plugs);
    printf("%s, ", cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
    print_cipher(cipher, cipher_len, NULL);
    printf(", ");
    print_text(decrypted, cipher_len);
    printf("\n");
}

// ------------------------------------------------------------------ wheel-order enumeration

// Enumerate the wheel orders phase 1 will try, writing full EnigmaKey templates (reflector,
// n_wheels, rotor[], greek for M4) into tmpl[] and returning the count. A pinned -rotors
// yields exactly one; otherwise all ordered triples of the default pool (60).
int enigma_enumerate_wheel_orders(const ColossusConfig *cfg, EnigmaKey *tmpl, int cap) {
    int n = 0;
    EnigmaKey base;
    memset(&base, 0, sizeof(base));
    base.reflector = cfg->enigma_reflector;
    base.n_wheels = cfg->enigma_model;   // 3 or 4
    enigma_plug_identity(base.plug);
    int off = base.n_wheels - 3;         // stepping-wheel base index (1 for M4)

    if (cfg->enigma_rotors_present) {
        base.rotor[0] = (base.n_wheels == 4) ? cfg->enigma_greek : cfg->enigma_rotors[0];
        for (int i = 0; i < 3; i++) base.rotor[off + i] = cfg->enigma_rotors[i];
        tmpl[n++] = base;
        return n;
    }
    // Blind: ordered triples of the pool. (M4 blind is impractical -- guarded by the caller.)
    if (base.n_wheels == 4) base.rotor[0] = cfg->enigma_greek;
    for (int a = 0; a < ENIGMA_POOL_SIZE; a++)
      for (int b = 0; b < ENIGMA_POOL_SIZE; b++) {
        if (b == a) continue;
        for (int c = 0; c < ENIGMA_POOL_SIZE; c++) {
            if (c == a || c == b) continue;
            if (n >= cap) return n;
            base.rotor[off + 0] = ENIGMA_POOL[a];
            base.rotor[off + 1] = ENIGMA_POOL[b];
            base.rotor[off + 2] = ENIGMA_POOL[c];
            tmpl[n++] = base;
        }
      }
    return n;
}

// ------------------------------------------------------------------ phases 2-3 driver

// Ring refinement + plugboard climb over a list of candidate base keys (each carrying a
// wheel order + start positions, rings AAA). Shared by the ciphertext-only path (bases from
// the phase-1 IoC search) and the Bombe (bases from the surviving stops). Rings/plugboard
// pins on cfg are honoured. Reports and fills *result.
void enigma_attack_from_bases(ColossusConfig *cfg, SharedData *shared,
    int cipher[], int cipher_len, char *cribtext,
    int crib_indices[], int crib_positions[], int n_cribs,
    const EnigmaKey bases[], int n_bases, int maxplugs, SolveResult *result) {

    EnigmaScratch scratch;
    scratch.maxplugs = maxplugs;
    scratch.n_cand = (n_bases > ENIGMA_MAX_CAND) ? ENIGMA_MAX_CAND : n_bases;
    for (int i = 0; i < scratch.n_cand; i++) { scratch.cand[i].key = bases[i]; scratch.cand[i].ioc = 0.0; }

    // Phase 2: ring settings (skip if pinned).
    for (int i = 0; i < scratch.n_cand; i++) {
        EnigmaKey *k = &scratch.cand[i].key;
        int r = k->n_wheels - 1, m = k->n_wheels - 2;
        if (cfg->enigma_ring_present) {
            int off = k->n_wheels - 3;
            for (int j = 0; j < 3; j++) k->ring[off + j] = cfg->enigma_ring[j];
        } else {
            enigma_refine_ring(k, r, cipher, cipher_len);   // fast ring
            enigma_refine_ring(k, m, cipher, cipher_len);   // middle ring
        }
        scratch.cand[i].ioc = enigma_ioc_of(k, cipher, cipher_len);
    }

    // Phase 3 plugboard, two tiers. TIER A ranks every candidate CHEAPLY: a fast greedy-add
    // plugboard (enigma_quick_plug) and its n-gram score. A pinned plugboard is scored as-is
    // and left fixed; a base that already carries steckers (the Bombe's recovered plugboard) is
    // scored as-is but still refined by the tier-B climb.
    static _Thread_local int qdec[MAX_CIPHER_LENGTH];
    double qscore[ENIGMA_MAX_CAND];
    bool   plug_fixed[ENIGMA_MAX_CAND];
    for (int i = 0; i < scratch.n_cand; i++) {
        EnigmaKey *k = &scratch.cand[i].key;
        int plug[26];
        bool have_seed = false;
        for (int j = 0; j < 26; j++) if (k->plug[j] != j) { have_seed = true; break; }
        plug_fixed[i] = cfg->enigma_plug_present;
        if (cfg->enigma_plug_present) {                     // pinned plugboard: score, don't climb
            memcpy(plug, cfg->enigma_plug, sizeof(int) * 26);
            qscore[i] = enigma_plug_score(k, plug, cipher, cipher_len,
                            shared->ngram_data, cfg->ngram_size, qdec);
        } else if (have_seed) {                             // Bombe stecker seed: score, refine below
            memcpy(plug, k->plug, sizeof(int) * 26);
            qscore[i] = enigma_plug_score(k, plug, cipher, cipher_len,
                            shared->ngram_data, cfg->ngram_size, qdec);
        } else {                                            // Gillogly greedy warm start
            qscore[i] = enigma_quick_plug(k, cipher, cipher_len,
                            shared->ngram_data, cfg->ngram_size, maxplugs, plug);
        }
        memcpy(k->plug, plug, sizeof(int) * 26);
    }

    // TIER B: spend the full reswap climb only on the best-ranked few. Ranking with the cheap
    // greedy first keeps the expensive climb off the no-hope orders that survived phase 1's IoC
    // cut (matters when -ntopk is large); by default all candidates fall within the cap.
    int rank[ENIGMA_MAX_CAND];
    for (int i = 0; i < scratch.n_cand; i++) rank[i] = i;
    for (int a = 0; a < scratch.n_cand; a++)               // selection sort by qscore desc (n small)
        for (int b = a + 1; b < scratch.n_cand; b++)
            if (qscore[rank[b]] > qscore[rank[a]]) { int t = rank[a]; rank[a] = rank[b]; rank[b] = t; }
    int n_climb = (scratch.n_cand < ENIGMA_PLUG_CLIMB_TOP) ? scratch.n_cand : ENIGMA_PLUG_CLIMB_TOP;
    for (int r = 0; r < n_climb; r++) {
        int i = rank[r];
        if (plug_fixed[i]) continue;                        // pinned board stays as pinned
        EnigmaKey *k = &scratch.cand[i].key;
        int plug[26];
        qscore[i] = enigma_plugboard_climb(k, cipher, cipher_len, shared->ngram_data,
                        cfg->ngram_size, maxplugs, k->plug, plug);   // seed from the tier-A board
        memcpy(k->plug, plug, sizeof(int) * 26);
    }

    if (cfg->verbose) {
        for (int i = 0; i < scratch.n_cand; i++) {
            EnigmaKey *k = &scratch.cand[i].key;
            char pl[128]; enigma_format_plugs(k->plug, pl);
            char ro[64] = {0}; for (int j = 0; j < k->n_wheels; j++) { strcat(ro, enigma_rotor_name(k->rotor[j])); strcat(ro, " "); }
            printf("  cand %d: %s ioc=%.4f score=%.4f plugs=%s\n", i, ro, scratch.cand[i].ioc, qscore[i], pl);
        }
    }

    // Phase 3: engine plugboard climb (n-gram + crib fitness), keeps the global best.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext,
        cipher, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.model_scratch = &scratch;
    ctx.result = result;
    run_solver(&ENIGMA_MODEL, &ctx);
}

// ------------------------------------------------------------------ entry point

void solve_enigma(char *ciphertext_str, char *cribtext_str,
    ColossusConfig *cfg, SharedData *shared,
    int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;
    enigma_init();

    if (g_alpha != ALPHABET_SIZE) {
        printf("\n\nERROR: Enigma needs the full 26-letter alphabet (got %d).\n\n", g_alpha);
        return;
    }
    if (cipher_len < 1) { printf("\n\nERROR: empty Enigma ciphertext.\n\n"); return; }
    for (int i = 0; i < cipher_len; i++)
        if (cipher_indices[i] < 0 || cipher_indices[i] >= 26) {
            printf("\n\nERROR: Enigma ciphertext must be solid letters (bad symbol at %d).\n\n", i);
            return;
        }

    // --- Bombe (crib) attack ---------------------------------------------------------
    if (cfg->enigma_bombe) {
        if (n_cribs <= 0) {
            printf("\n\nERROR: -bombe requires a crib (-crib). No crib supplied.\n\n");
            return;
        }
        solve_enigma_bombe(cfg, shared, cipher_indices, cipher_len,
                           crib_indices, crib_positions, n_cribs, result);
        return;
    }

    int maxplugs = (cfg->enigma_maxplugs > 0) ? cfg->enigma_maxplugs : 10;
    if (maxplugs > 13) maxplugs = 13;

    // --- Known-key decrypt: rotors + ring + pos + plug all pinned ---------------------
    if (cfg->enigma_rotors_present && cfg->enigma_ring_present &&
        cfg->enigma_pos_present && cfg->enigma_plug_present) {
        EnigmaKey k;
        memset(&k, 0, sizeof(k));
        k.reflector = cfg->enigma_reflector;
        k.n_wheels = cfg->enigma_model;
        int off = k.n_wheels - 3;
        if (k.n_wheels == 4) { k.rotor[0] = cfg->enigma_greek; k.ring[0] = 0; k.pos[0] = 0; }
        for (int i = 0; i < 3; i++) {
            k.rotor[off + i] = cfg->enigma_rotors[i];
            k.ring[off + i] = cfg->enigma_ring[i];
            k.pos[off + i]  = cfg->enigma_pos[i];
        }
        memcpy(k.plug, cfg->enigma_plug, sizeof(int) * 26);
        static int dec[MAX_CIPHER_LENGTH];
        enigma_encrypt(cipher_indices, cipher_len, &k, dec);
        printf("\nenigma: known-key decrypt (all settings pinned)\n");
        enigma_emit_report(cfg, shared, cipher_indices, cipher_len, &k, dec, cribtext_str,
                           0.0, result);
        return;
    }

    if (cfg->enigma_model == 4 && !cfg->enigma_rotors_present) {
        printf("\n\nERROR: blind M4 (4-rotor) ciphertext-only search is impractical "
               "(26^4 x hundreds of orders). Pin the wheel order with -rotors, or use "
               "-bombe with a crib.\n\n");
        return;
    }

    // --- Phase 1: wheel order + fast ring + start positions by IoC (threaded) ---------
    // A threaded (order x fast-ring x position) IoC sweep keeps the top-P configs; each is
    // then middle-ring-refined and the global top-K by that IoC feed the plugboard climb.
    static EnigmaKey tmpl[64];
    int n_orders = enigma_enumerate_wheel_orders(cfg, tmpl, 64);
    int K = (cfg->enigma_ntopk > 0) ? cfg->enigma_ntopk : 6;
    if (K > ENIGMA_MAX_CAND) K = ENIGMA_MAX_CAND;

    EnigmaScratch scratch;
    scratch.n_cand = 0;
    scratch.maxplugs = maxplugs;

    int rr_base = cfg->enigma_ring_present ? cfg->enigma_ring[2] : 0;
    int nthreads = cfg->n_threads > 0 ? cfg->n_threads : 1;

    printf("\nenigma: ciphertext-only attack, %d ciphertext letters, %d wheel order(s), "
           "top-K=%d, maxplugs=%d, %d thread(s)\n",
           cipher_len, n_orders, K, maxplugs, nthreads);

    static EnigmaKey posA[512], posB[512];
    double iocA[512], iocB[512];

    // Pass 1a: cheap rings-coarse position sweep to rank wheel orders (one fast-ring value).
    int nA = enigma_phase1_search(cfg, tmpl, n_orders, 1, rr_base,
                                  cipher_indices, cipher_len, posA, iocA, 512);

    // Pass 1b: the full fast-ring search, but only on the distinct top orders from 1a (unless
    // the ring is pinned, in which case 1a already used it -- no fast-ring sweep needed).
    int nB;
    if (cfg->enigma_ring_present) {
        nB = nA;
        for (int i = 0; i < nA; i++) { posB[i] = posA[i]; iocB[i] = iocA[i]; }
    } else {
        static EnigmaKey tmpl2[ENIGMA_PHASE1_ORDERS];
        int n2 = enigma_distinct_orders(posA, nA, tmpl2, ENIGMA_PHASE1_ORDERS);
        nB = enigma_phase1_search(cfg, tmpl2, n2, 26, 0,
                                  cipher_indices, cipher_len, posB, iocB, 512);
    }

    // Middle-ring-refine each candidate and keep the global top-K by the refined IoC.
    for (int i = 0; i < nB; i++) {
        EnigmaKey k = posB[i];
        if (!cfg->enigma_ring_present)
            enigma_refine_ring(&k, k.n_wheels - 2, cipher_indices, cipher_len);   // middle ring
        double rioc = enigma_ioc_of(&k, cipher_indices, cipher_len);
        enigma_topk_insert(&scratch, &k, rioc, K);
    }

    // --- Phases 2-3: ring refine + plugboard climb on the top-K candidates ------------
    EnigmaKey bases[ENIGMA_MAX_CAND];
    for (int i = 0; i < scratch.n_cand; i++) bases[i] = scratch.cand[i].key;
    enigma_attack_from_bases(cfg, shared, cipher_indices, cipher_len, cribtext_str,
        crib_indices, crib_positions, n_cribs, bases, scratch.n_cand, maxplugs, result);
}

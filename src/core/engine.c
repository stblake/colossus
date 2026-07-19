#include "engine.h"
#include "scoring.h"
#include <pthread.h>




// =====================================================================
//  Cipher-type-agnostic search engine
// =====================================================================
//
// run_solver() drives every cipher type through one skeleton; the per-type
// CipherModel (colossus.h) supplies the cipher specifics as hooks. See the header
// for the interface contract. The big SolverState buffers are file-static (the
// program is single-threaded -- rng_state is a global) so unifying them here does
// not grow the stack.

// Build the optimal-cycleword per-column ciphertext histogram for a given period
// into ctx->hist_by_col (laid out as hist_by_col[col*ALPHABET_SIZE + c]). Depends
// only on the fixed ciphertext and the period, so the engine builds it once per
// config rather than on every derive_optimal_cycleword call.
static void engine_build_hist(SolverCtx *ctx, int period) {
    if (ctx->hist_by_col == NULL || period <= 0) return;
    for (int i = 0; i < period * ALPHABET_SIZE; i++) ctx->hist_by_col[i] = 0;
    int col = 0;
    for (int i = 0; i < ctx->cipher_len; i++) {
        int c = ctx->cipher[i];
        if (c >= 0) ctx->hist_by_col[col * ALPHABET_SIZE + c]++;
        if (++col == period) col = 0;
    }
}

static double engine_score(const SolverCtx *ctx, int *decrypted, double adjust) {
    ColossusConfig *cfg = ctx->cfg;
    return state_score(decrypted, ctx->cipher_len,
        ctx->crib_indices, ctx->crib_positions, ctx->n_cribs,
        ctx->ngram_data, cfg->ngram_size,
        cfg->weight_ngram, cfg->weight_crib, cfg->weight_ioc, cfg->weight_entropy) + adjust;
}

// =====================================================================
//  Restart-loop parallelism (-nthreads)
// =====================================================================
//
// Each of the three drivers below (generic / incremental / PSO) splits its restart
// loop into per-thread "restart ranges", each run on a PRIVATE workspace and a
// THREAD-LOCAL RNG (rng_state, seeded per worker). The workers merge into one shared
// EngineGlobalBest under a mutex, which also serialises -verbose best-improvement
// logging so the lines never interleave and the screen shows monotonic GLOBAL bests.
//
// With -nthreads <= 1 no threads are spawned: the driver runs a single restart range
// on the calling (main) thread, over the file-static workspace, with the main-thread
// RNG and mtx == NULL -- i.e. the original sequential path, so every fixed-seed solve
// is bit-identical to before. The -nrestarts budget is SPLIT across the workers (same
// total search, ~T x faster), each worker taking a contiguous slice.

typedef struct {
    SolverState state;
    int decrypted[MAX_CIPHER_LENGTH];
    double score;
    bool have;
} EngineGlobalBest;

// Merge a candidate into the shared global best (serialised when mtx != NULL). Emits
// the -verbose report ONLY when the candidate improves the global best AND the caller
// permits it (allow_report is false for PSO's swarm-seed pass, which historically did
// not log), so the on-screen best sequence is unchanged from the sequential path.
static void engine_publish_best(const CipherModel *m, SolverCtx *ctx,
        const SolverConfig *cfg_c, const SolverState *cand, double score,
        int *decrypted, EngineStats *st, EngineGlobalBest *gb,
        pthread_mutex_t *mtx, bool verbose, bool allow_report) {
    if (mtx) pthread_mutex_lock(mtx);
    if (!gb->have || score > gb->score) {
        gb->have = true;
        gb->score = score;
        m->copy_state(cfg_c, cand, &gb->state);
        vec_copy(decrypted, gb->decrypted, ctx->cipher_len);
        if (verbose && allow_report && m->report_verbose)
            m->report_verbose(ctx, cfg_c, &gb->state, score, decrypted, st);
    }
    if (mtx) pthread_mutex_unlock(mtx);
}

// Worker-count for the restart split: clamp to [1, n_restarts]; deterministic shapes
// (which early-exit on the first best) stay single-threaded.
static int engine_plan_threads(int n_threads, int n_restarts, bool deterministic) {
    int T = n_threads;
    if (T < 1) T = 1;
    if (deterministic) T = 1;
    if (T > n_restarts) T = n_restarts;
    if (T < 1) T = 1;
    return T;
}

// Incremental variant of run_one_config (same restart / annealing / backtrack /
// best-tracking skeleton) for models that supply the score_neighbor/commit/sync
// hooks. The current state's decryption (cur_dec) is kept live; each neighbour is
// scored as a delta and only committed into cur_dec + the model caches on accept,
// so the per-iteration cost is O(positions the move touched) rather than O(N).
typedef struct { SolverState cur, loc, best;
                 int cur_dec[MAX_CIPHER_LENGTH], best_dec[MAX_CIPHER_LENGTH]; } IncWorkspace;

// One incremental restart range [rs_begin, rs_end) on a private workspace, merging
// its bests into gb. Mirrors the original run_one_config_incremental restart body.
static void inc_restart_range(const CipherModel *m, SolverCtx *ctx,
        const SolverConfig *cfg_c, SearchShape shape, bool deterministic,
        double temp_start, double cooling, int rs_begin, int rs_end,
        IncWorkspace *ws, EngineGlobalBest *gb, pthread_mutex_t *mtx) {
    ColossusConfig *cfg = ctx->cfg;
    double best_score = 0.0, cur_score, loc_score, adjust;
    bool have_best = false, force_primary = true, done = false;

    EngineStats st;
    memset(&st, 0, sizeof st);
    st.start_time = wall_time_sec();

    for (int rs = rs_begin; rs < rs_end && !done; rs++) {
        st.n_restarts = rs;

        if (have_best && frand() < cfg->backtracking_probability) {
            m->copy_state(cfg_c, &ws->best, &ws->cur);
            vec_copy(ws->best_dec, ws->cur_dec, ctx->cipher_len);
            cur_score = best_score;
            st.n_backtracks++;
        } else {
            m->seed(ctx, cfg_c, &ws->cur);
            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &ws->cur, ws->cur_dec, &adjust);
            cur_score = engine_score(ctx, ws->cur_dec, adjust);
        }
        // (Re)build the model caches so they describe the current decryption.
        m->sync_caches(ctx, cfg_c, ws->cur_dec);

        force_primary = true;

        double temp = temp_start;
        for (int it = 0; it < cfg->n_hill_climbs && !done; it++) {
            st.n_iterations++;

            m->copy_state(cfg_c, &ws->cur, &ws->loc);
            m->perturb(ctx, cfg_c, &ws->loc, &force_primary);

            loc_score = m->score_neighbor(ctx, cfg_c, &ws->cur, &ws->loc, ws->cur_dec, cur_score);

            bool accept;
            if (loc_score > cur_score) {
                accept = true;
            } else if (shape == SHAPE_ANNEAL) {
                accept = frand() < exp((loc_score - cur_score) / temp);
            } else {
                accept = frand() < cfg->slip_probability;
            }
            if (accept) {
                if (loc_score <= cur_score) st.n_slips++;
                m->commit_neighbor(ctx, cfg_c, ws->cur_dec);   // advance cur_dec + caches
                m->copy_state(cfg_c, &ws->loc, &ws->cur);
                cur_score = loc_score;
            }
            temp *= cooling;

            if (!have_best || cur_score > best_score) {
                best_score = cur_score; have_best = true;
                m->copy_state(cfg_c, &ws->cur, &ws->best);
                vec_copy(ws->cur_dec, ws->best_dec, ctx->cipher_len);
                engine_publish_best(m, ctx, cfg_c, &ws->cur, cur_score,
                                    ws->cur_dec, &st, gb, mtx, cfg->verbose, true);
                if (deterministic) done = true;
            }
        }
    }
}

typedef struct {
    const CipherModel *m; const SolverConfig *cfg_c;
    SearchShape shape; bool deterministic; double temp_start, cooling;
    int rs_begin, rs_end; uint32_t seed;
    IncWorkspace *ws; EngineGlobalBest *gb; pthread_mutex_t *mtx;
    SolverCtx ctx_copy;   // shallow copy of ctx, repointed at this worker's scratch
} IncThreadArg;

static void *inc_thread_main(void *p) {
    IncThreadArg *a = (IncThreadArg *)p;
    rng_state = a->seed;
    inc_restart_range(a->m, &a->ctx_copy, a->cfg_c, a->shape, a->deterministic,
                      a->temp_start, a->cooling, a->rs_begin, a->rs_end,
                      a->ws, a->gb, a->mtx);
    return NULL;
}

// The incremental fast-path supports -nthreads only for models that ALSO supply the
// scratch_clone/scratch_free hooks: its live neighbour caches live in
// ctx->model_scratch (a single shared pointer), so each worker needs its own private
// scratch (a shallow ctx copy repointed at a per-thread clone), exactly like the
// generic/PSO drivers give each worker its own workspace. A model that leaves the
// clone hooks NULL (or -nthreads <= 1) takes the verbatim sequential path -- so every
// fixed-seed solve stays bit-identical there.
static double run_one_config_incremental(const CipherModel *m, SolverCtx *ctx,
                                         const SolverConfig *cfg_c,
                                         SolverState *out_best, int *out_decrypted) {
    ColossusConfig *cfg = ctx->cfg;

    SearchShape shape = m->shape;
    if (cfg->method == METHOD_SHOTGUN) shape = SHAPE_SHOTGUN;
    else if (cfg->method == METHOD_ANNEAL) shape = SHAPE_ANNEAL;
    bool deterministic = (shape == SHAPE_DETERMINISTIC);

    double temp_start = cfg->init_temp;
    double cooling;
    if (cfg->cooling_rate > 0.0) {
        cooling = cfg->cooling_rate;
    } else {
        cooling = 1.0;
        if (cfg->n_hill_climbs > 1)
            cooling = pow(cfg->min_temp / temp_start, 1.0 / (double)(cfg->n_hill_climbs - 1));
    }

    static EngineGlobalBest gb;
    gb.have = false; gb.score = 0.0;

    int T = engine_plan_threads(cfg->n_threads, cfg->n_restarts, deterministic);
    if (T <= 1 || !m->scratch_clone || !m->scratch_free) {
        // Sequential path: single restart range on the calling thread, over the
        // file-static workspace + the owner scratch, main-thread RNG, no mutex ->
        // bit-identical to before.
        static IncWorkspace ws;
        inc_restart_range(m, ctx, cfg_c, shape, deterministic, temp_start, cooling,
                          0, cfg->n_restarts, &ws, &gb, NULL);
    } else {
        uint32_t base_seed = fast_rand();
        pthread_mutex_t mtx; pthread_mutex_init(&mtx, NULL);
        pthread_t *tids = malloc((size_t)T * sizeof *tids);
        IncThreadArg *args = malloc((size_t)T * sizeof *args);
        IncWorkspace *wss = malloc((size_t)T * sizeof *wss);
        bool ok = tids && args && wss;
        if (ok)
            for (int t = 0; t < T; t++) args[t].ctx_copy.model_scratch = NULL;  // = "no clone yet"
        if (ok)
            for (int t = 0; t < T; t++) {
                void *sc = m->scratch_clone(ctx);
                if (!sc) { ok = false; break; }
                args[t].ctx_copy = *ctx;                 // shallow copy of the shared ctx
                args[t].ctx_copy.model_scratch = sc;     // ... repointed at this worker's scratch
            }
        if (ok) {
            int per = cfg->n_restarts / T, rem = cfg->n_restarts % T, start = 0;
            for (int t = 0; t < T; t++) {
                int cnt = per + (t < rem ? 1 : 0);   // contiguous split of [0, n_restarts)
                args[t].m = m; args[t].cfg_c = cfg_c; args[t].shape = shape;
                args[t].deterministic = deterministic;
                args[t].temp_start = temp_start; args[t].cooling = cooling;
                args[t].rs_begin = start; args[t].rs_end = start + cnt;
                args[t].seed = rng_seed_thread(base_seed, t);
                args[t].ws = &wss[t]; args[t].gb = &gb; args[t].mtx = &mtx;
                start += cnt;
                pthread_create(&tids[t], NULL, inc_thread_main, &args[t]);
            }
            for (int t = 0; t < T; t++) pthread_join(tids[t], NULL);
        } else {
            // Allocation/clone failed: fall back to one sequential range (owner scratch).
            static IncWorkspace fb_ws;
            inc_restart_range(m, ctx, cfg_c, shape, deterministic, temp_start, cooling,
                              0, cfg->n_restarts, &fb_ws, &gb, NULL);
        }
        if (args)
            for (int t = 0; t < T; t++)
                if (args[t].ctx_copy.model_scratch) m->scratch_free(args[t].ctx_copy.model_scratch);
        free(tids); free(args); free(wss);
        pthread_mutex_destroy(&mtx);
    }

    if (!gb.have) return 0.0;
    m->copy_state(cfg_c, &gb.state, out_best);
    vec_copy(gb.decrypted, out_decrypted, ctx->cipher_len);
    return gb.score;
}

// =====================================================================
//  Particle swarm optimisation (SHAPE_PSO / -method pso)
// =====================================================================
//
// A third optimisation method, sibling to the shotgun / anneal climbers above,
// and -- like them -- COMPLETELY CIPHER-AGNOSTIC: it drives every cipher type
// through nothing but the model's existing hooks (seed / perturb / copy_state /
// decrypt) plus a generic Hamming distance over the raw SolverState lanes. It
// never interprets the representation, so a permutation stays a permutation, a
// keyword stays a keyword, a homophone map stays a map.
//
// Discrete swarm: a particle's "position" IS a SolverState. The two PSO
// primitives are built from the model's own neighbour operator:
//   * "pull toward an attractor (pbest/gbest)" = apply perturb() moves and keep
//     only those that do not increase the Hamming distance to the attractor
//     (pull_toward). Every kept move is the model's own move, so validity is
//     automatic.
//   * "inertia / momentum" = a few plain random perturb() moves.
// Each particle then does a short greedy local refinement (memetic PSO) before
// its decrypt+score updates the personal best (pbest) and the global best
// (gbest). The swarm reuses the engine budget knobs: n_particles particles run
// for n_hill_climbs iterations, relaunched n_restarts times.

#define MAX_PSO_PARTICLES 128   // swarm-size cap (n_particles is clamped to this)
#define PSO_PULL_ATTEMPTS  12   // tries to find a distance-reducing move per pull step

// Generic Hamming distance between two states: count of differing entries over the
// active lanes the engine already tracks (the fixed 26-entry keyword lanes, the
// cycleword to its length, the key lane to key_len). For polyalphabetic types the
// cycleword length is the config's `period`; non-polyalpha types never touch the
// cycleword lane (it stays zero across all particles) so comparing it is a no-op.
// Dead lanes are seeded identically and never moved by perturb(), so they
// contribute nothing -- this only needs to be a monotone proxy, not an exact metric.
static int state_distance(const SolverConfig *cfg_c,
                          const SolverState *a, const SolverState *b) {
    int d = 0;
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (a->pt_keyword[i] != b->pt_keyword[i]) d++;
        if (a->ct_keyword[i] != b->ct_keyword[i]) d++;
    }
    int cw = cfg_c->period;
    if (cw < 0) cw = 0;
    if (cw > MAX_CYCLEWORD_LEN) cw = MAX_CYCLEWORD_LEN;
    for (int i = 0; i < cw; i++)
        if (a->cycleword[i] != b->cycleword[i]) d++;
    int kl = a->key_len < b->key_len ? a->key_len : b->key_len;
    for (int i = 0; i < kl; i++)
        if (a->key[i] != b->key[i]) d++;
    return d;
}

// Move `st` up to n_moves of the model's own perturb() moves toward `target`,
// keeping only strictly-closer moves (a few attempts each); stops early once it
// can find no closer neighbour. `trial` is caller-provided scratch. Pulls never
// decrypt -- they only compare Hamming distances -- so they are cheap.
static void pull_toward(const CipherModel *m, const SolverCtx *ctx,
                        const SolverConfig *cfg_c, SolverState *st,
                        const SolverState *target, int n_moves, SolverState *trial) {
    bool force_primary = false;
    int dist = state_distance(cfg_c, st, target);
    for (int mv = 0; mv < n_moves && dist > 0; mv++) {
        bool accepted = false;
        for (int attempt = 0; attempt < PSO_PULL_ATTEMPTS; attempt++) {
            m->copy_state(cfg_c, st, trial);
            m->perturb(ctx, cfg_c, trial, &force_primary);
            int nd = state_distance(cfg_c, trial, target);
            if (nd < dist) {
                m->copy_state(cfg_c, trial, st);
                dist = nd;
                accepted = true;
                break;
            }
        }
        if (!accepted) break;   // stuck: no closer neighbour found, stop pulling
    }
}

// Run one outer config under particle-swarm optimisation. Same contract as
// run_one_config: writes the best state + its decryption to the out params and
// returns the best score. Always uses the full decrypt+score path (the optional
// incremental fast-path hooks are not used here).
typedef struct {
    SolverState *particle;    // [np]
    SolverState *pbest;       // [np]
    double      *pbest_score; // [np]
    SolverState gbest, trial, loc;
    int decrypted[MAX_CIPHER_LENGTH], gbest_dec[MAX_CIPHER_LENGTH];
} PsoWorkspace;

// One PSO swarm-relaunch range [rs_begin, rs_end) on a private workspace, merging its
// global bests into gb. Mirrors the original run_one_config_pso relaunch body. The
// swarm-seed pass publishes to gb but does NOT log (allow_report=false), matching the
// original which only reported on the social-update improvement.
static void pso_restart_range(const CipherModel *m, SolverCtx *ctx,
        const SolverConfig *cfg_c, int np, int refine, int rs_begin, int rs_end,
        PsoWorkspace *ws, EngineGlobalBest *gb, pthread_mutex_t *mtx) {
    ColossusConfig *cfg = ctx->cfg;
    double gbest_score = 0.0, adjust;
    bool have_gbest = false, force_primary = false;

    EngineStats st;
    memset(&st, 0, sizeof st);
    st.start_time = wall_time_sec();

    for (int rs = rs_begin; rs < rs_end; rs++) {
        st.n_restarts = rs;

        // (Re)seed the swarm and initialise personal / global bests.
        for (int p = 0; p < np; p++) {
            m->seed(ctx, cfg_c, &ws->particle[p]);
            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &ws->particle[p], ws->decrypted, &adjust);
            double sc = engine_score(ctx, ws->decrypted, adjust);
            m->copy_state(cfg_c, &ws->particle[p], &ws->pbest[p]);
            ws->pbest_score[p] = sc;
            if (!have_gbest || sc > gbest_score) {
                gbest_score = sc; have_gbest = true;
                m->copy_state(cfg_c, &ws->particle[p], &ws->gbest);
                vec_copy(ws->decrypted, ws->gbest_dec, ctx->cipher_len);
                engine_publish_best(m, ctx, cfg_c, &ws->gbest, gbest_score,
                                    ws->gbest_dec, &st, gb, mtx, cfg->verbose, false);
            }
        }

        for (int it = 0; it < cfg->n_hill_climbs; it++) {
            st.n_iterations++;
            for (int p = 0; p < np; p++) {
                // 1. inertia: a little random momentum (exploration).
                int n_inertia = (int)(cfg->inertia * (1.0 + frand()) + 0.5);
                for (int k = 0; k < n_inertia; k++)
                    m->perturb(ctx, cfg_c, &ws->particle[p], &force_primary);

                // 2. cognitive pull toward this particle's personal best.
                int n_cog = (int)(cfg->cognitive * frand() *
                                  state_distance(cfg_c, &ws->particle[p], &ws->pbest[p]) + 0.5);
                if (n_cog > 0)
                    pull_toward(m, ctx, cfg_c, &ws->particle[p], &ws->pbest[p], n_cog, &ws->trial);

                // 3. social pull toward the global best.
                int n_soc = (int)(cfg->social * frand() *
                                  state_distance(cfg_c, &ws->particle[p], &ws->gbest) + 0.5);
                if (n_soc > 0)
                    pull_toward(m, ctx, cfg_c, &ws->particle[p], &ws->gbest, n_soc, &ws->trial);

                // 4. memetic local refinement: a short greedy hill-climb.
                adjust = 0.0;
                m->decrypt(ctx, cfg_c, &ws->particle[p], ws->decrypted, &adjust);
                double sc = engine_score(ctx, ws->decrypted, adjust);
                for (int k = 0; k < refine; k++) {
                    m->copy_state(cfg_c, &ws->particle[p], &ws->loc);
                    m->perturb(ctx, cfg_c, &ws->loc, &force_primary);
                    double adj2 = 0.0;
                    m->decrypt(ctx, cfg_c, &ws->loc, ws->decrypted, &adj2);
                    double ls = engine_score(ctx, ws->decrypted, adj2);
                    if (ls > sc) { m->copy_state(cfg_c, &ws->loc, &ws->particle[p]); sc = ls; }
                }

                // 5. update personal and global bests.
                if (sc > ws->pbest_score[p]) {
                    m->copy_state(cfg_c, &ws->particle[p], &ws->pbest[p]);
                    ws->pbest_score[p] = sc;
                    if (sc > gbest_score) {
                        gbest_score = sc;
                        m->copy_state(cfg_c, &ws->particle[p], &ws->gbest);
                        adjust = 0.0;
                        m->decrypt(ctx, cfg_c, &ws->gbest, ws->gbest_dec, &adjust);
                        engine_publish_best(m, ctx, cfg_c, &ws->gbest, gbest_score,
                                            ws->gbest_dec, &st, gb, mtx, cfg->verbose, true);
                    }
                }
            }
        }
    }
}

typedef struct {
    const CipherModel *m; SolverCtx *ctx; const SolverConfig *cfg_c;
    int np, refine, rs_begin, rs_end; uint32_t seed;
    PsoWorkspace *ws; EngineGlobalBest *gb; pthread_mutex_t *mtx;
} PsoThreadArg;

static void *pso_thread_main(void *p) {
    PsoThreadArg *a = (PsoThreadArg *)p;
    rng_state = a->seed;
    pso_restart_range(a->m, a->ctx, a->cfg_c, a->np, a->refine,
                      a->rs_begin, a->rs_end, a->ws, a->gb, a->mtx);
    return NULL;
}

static double run_one_config_pso(const CipherModel *m, SolverCtx *ctx,
                                 const SolverConfig *cfg_c,
                                 SolverState *out_best, int *out_decrypted) {
    ColossusConfig *cfg = ctx->cfg;

    int np = cfg->n_particles;
    if (np < 1) np = 1;
    if (np > MAX_PSO_PARTICLES) np = MAX_PSO_PARTICLES;
    int refine = cfg->refine_steps < 0 ? 0 : cfg->refine_steps;

    static EngineGlobalBest gb;
    gb.have = false; gb.score = 0.0;

    // PSO is never a model default (reached only via -method pso), so shape is PSO,
    // never deterministic -> the restart split is purely over swarm relaunches.
    int T = engine_plan_threads(cfg->n_threads, cfg->n_restarts, false);
    if (T <= 1) {
        static SolverState s_particle[MAX_PSO_PARTICLES], s_pbest[MAX_PSO_PARTICLES];
        static double s_pbest_score[MAX_PSO_PARTICLES];
        static PsoWorkspace ws;
        ws.particle = s_particle; ws.pbest = s_pbest; ws.pbest_score = s_pbest_score;
        pso_restart_range(m, ctx, cfg_c, np, refine, 0, cfg->n_restarts, &ws, &gb, NULL);
    } else {
        uint32_t base_seed = fast_rand();
        pthread_mutex_t mtx; pthread_mutex_init(&mtx, NULL);
        pthread_t *tids = malloc((size_t)T * sizeof *tids);
        PsoThreadArg *args = malloc((size_t)T * sizeof *args);
        PsoWorkspace *wss = malloc((size_t)T * sizeof *wss);
        for (int t = 0; t < T; t++) {   // per-thread swarm arrays
            wss[t].particle = malloc((size_t)np * sizeof(SolverState));
            wss[t].pbest = malloc((size_t)np * sizeof(SolverState));
            wss[t].pbest_score = malloc((size_t)np * sizeof(double));
        }
        int per = cfg->n_restarts / T, rem = cfg->n_restarts % T, start = 0;
        for (int t = 0; t < T; t++) {
            int cnt = per + (t < rem ? 1 : 0);
            args[t] = (PsoThreadArg){ .m = m, .ctx = ctx, .cfg_c = cfg_c, .np = np,
                .refine = refine, .rs_begin = start, .rs_end = start + cnt,
                .seed = rng_seed_thread(base_seed, t), .ws = &wss[t], .gb = &gb, .mtx = &mtx };
            start += cnt;
            pthread_create(&tids[t], NULL, pso_thread_main, &args[t]);
        }
        for (int t = 0; t < T; t++) pthread_join(tids[t], NULL);
        for (int t = 0; t < T; t++) { free(wss[t].particle); free(wss[t].pbest); free(wss[t].pbest_score); }
        free(tids); free(args); free(wss);
        pthread_mutex_destroy(&mtx);
    }

    if (!gb.have) return 0.0;
    m->copy_state(cfg_c, &gb.state, out_best);
    double adjust = 0.0;
    m->decrypt(ctx, cfg_c, &gb.state, out_decrypted, &adjust);
    return gb.score;
}

// Hill-climb one outer config: shotgun restarts + per-iteration neighbour move,
// with SHOTGUN slip / ANNEAL Metropolis acceptance and best-state tracking. Writes
// the config's best state to *out_best and its decryption to out_decrypted, and
// returns the best score.
typedef struct { SolverState cur, loc, best; int decrypted[MAX_CIPHER_LENGTH]; } GenWorkspace;

// One generic (non-incremental) restart range [rs_begin, rs_end) on a private
// workspace, merging its bests into gb. Mirrors the original run_one_config restart
// body; the "restart-0 backtrack draws no RNG" invariant holds because have_best
// starts false (best is only recorded inside the inner loop).
static void gen_restart_range(const CipherModel *m, SolverCtx *ctx,
        const SolverConfig *cfg_c, SearchShape shape, bool deterministic,
        double temp_start, double cooling, int rs_begin, int rs_end,
        GenWorkspace *ws, EngineGlobalBest *gb, pthread_mutex_t *mtx) {
    ColossusConfig *cfg = ctx->cfg;
    double best_score = 0.0, cur_score, loc_score, adjust;
    bool have_best = false, force_primary = true, done = false;

    EngineStats st;
    memset(&st, 0, sizeof st);
    st.start_time = wall_time_sec();

    for (int rs = rs_begin; rs < rs_end && !done; rs++) {
        st.n_restarts = rs;

        if (have_best && frand() < cfg->backtracking_probability) {
            m->copy_state(cfg_c, &ws->best, &ws->cur);
            cur_score = best_score;
            st.n_backtracks++;
        } else {
            m->seed(ctx, cfg_c, &ws->cur);
            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &ws->cur, ws->decrypted, &adjust);
            cur_score = engine_score(ctx, ws->decrypted, adjust);
        }

        // force_primary is the per-restart "must perturb the primary lane" flag
        // (polyalpha's perturbate_keyword_p); the model reads and updates it.
        force_primary = true;

        double temp = temp_start;
        for (int it = 0; it < cfg->n_hill_climbs && !done; it++) {
            st.n_iterations++;

            m->copy_state(cfg_c, &ws->cur, &ws->loc);

            m->perturb(ctx, cfg_c, &ws->loc, &force_primary);

            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &ws->loc, ws->decrypted, &adjust);
            loc_score = engine_score(ctx, ws->decrypted, adjust);

            bool accept;
            if (loc_score > cur_score) {
                accept = true;
            } else if (shape == SHAPE_ANNEAL) {
                accept = frand() < exp((loc_score - cur_score) / temp);
            } else {
                accept = frand() < cfg->slip_probability;
            }
            if (accept) {
                if (loc_score <= cur_score) st.n_slips++;
                m->copy_state(cfg_c, &ws->loc, &ws->cur);
                cur_score = loc_score;
            }
            temp *= cooling;

            if (!have_best || cur_score > best_score) {
                best_score = cur_score; have_best = true;
                m->copy_state(cfg_c, &ws->cur, &ws->best);
                engine_publish_best(m, ctx, cfg_c, &ws->cur, cur_score,
                                    ws->decrypted, &st, gb, mtx, cfg->verbose, true);
                if (deterministic) done = true;
            }
        }
    }
}

typedef struct {
    const CipherModel *m; SolverCtx *ctx; const SolverConfig *cfg_c;
    SearchShape shape; bool deterministic; double temp_start, cooling;
    int rs_begin, rs_end; uint32_t seed;
    GenWorkspace *ws; EngineGlobalBest *gb; pthread_mutex_t *mtx;
} GenThreadArg;

static void *gen_thread_main(void *p) {
    GenThreadArg *a = (GenThreadArg *)p;
    rng_state = a->seed;
    gen_restart_range(a->m, a->ctx, a->cfg_c, a->shape, a->deterministic,
                      a->temp_start, a->cooling, a->rs_begin, a->rs_end,
                      a->ws, a->gb, a->mtx);
    return NULL;
}

static double run_one_config(const CipherModel *m, SolverCtx *ctx,
                             const SolverConfig *cfg_c,
                             SolverState *out_best, int *out_decrypted) {

    // Particle swarm (-method pso) is a fully separate, cipher-agnostic driver.
    if (ctx->cfg->method == METHOD_PSO)
        return run_one_config_pso(m, ctx, cfg_c, out_best, out_decrypted);

    // Models exposing the incremental hooks take the fast path; all others fall
    // through to the unchanged generic climber below.
    if (m->score_neighbor && m->commit_neighbor && m->sync_caches)
        return run_one_config_incremental(m, ctx, cfg_c, out_best, out_decrypted);

    ColossusConfig *cfg = ctx->cfg;

    // -method overrides the model's built-in shape on EVERY cipher type (the engine
    // is acceptance-strategy agnostic); METHOD_DEFAULT keeps the model's own shape.
    SearchShape shape = m->shape;
    if (cfg->method == METHOD_SHOTGUN) shape = SHAPE_SHOTGUN;
    else if (cfg->method == METHOD_ANNEAL) shape = SHAPE_ANNEAL;
    bool deterministic = (shape == SHAPE_DETERMINISTIC);

    // Geometric Metropolis annealing schedule (used only by SHAPE_ANNEAL). The
    // start temperature and cooling come from cfg (-inittemp / -coolingrate); when
    // no cooling rate is given it is derived to cool init_temp -> min_temp over the
    // hill-climb. Defaults reproduce the previously hardcoded 0.10 -> 0.001 schedule.
    double temp_start = cfg->init_temp;
    double cooling;
    if (cfg->cooling_rate > 0.0) {
        cooling = cfg->cooling_rate;
    } else {
        cooling = 1.0;
        if (cfg->n_hill_climbs > 1)
            cooling = pow(cfg->min_temp / temp_start, 1.0 / (double)(cfg->n_hill_climbs - 1));
    }

    static EngineGlobalBest gb;
    gb.have = false; gb.score = 0.0;

    int T = engine_plan_threads(cfg->n_threads, cfg->n_restarts, deterministic);
    if (T <= 1) {
        // Sequential path: single restart range on the calling thread, over the
        // file-static workspace, main-thread RNG, no mutex -> bit-identical to before.
        static GenWorkspace ws;
        gen_restart_range(m, ctx, cfg_c, shape, deterministic, temp_start, cooling,
                          0, cfg->n_restarts, &ws, &gb, NULL);
    } else {
        uint32_t base_seed = fast_rand();
        pthread_mutex_t mtx; pthread_mutex_init(&mtx, NULL);
        pthread_t *tids = malloc((size_t)T * sizeof *tids);
        GenThreadArg *args = malloc((size_t)T * sizeof *args);
        GenWorkspace *wss = malloc((size_t)T * sizeof *wss);
        int per = cfg->n_restarts / T, rem = cfg->n_restarts % T, start = 0;
        for (int t = 0; t < T; t++) {
            int cnt = per + (t < rem ? 1 : 0);   // contiguous split of [0, n_restarts)
            args[t] = (GenThreadArg){ .m = m, .ctx = ctx, .cfg_c = cfg_c, .shape = shape,
                .deterministic = deterministic, .temp_start = temp_start, .cooling = cooling,
                .rs_begin = start, .rs_end = start + cnt, .seed = rng_seed_thread(base_seed, t),
                .ws = &wss[t], .gb = &gb, .mtx = &mtx };
            start += cnt;
            pthread_create(&tids[t], NULL, gen_thread_main, &args[t]);
        }
        for (int t = 0; t < T; t++) pthread_join(tids[t], NULL);
        free(tids); free(args); free(wss);
        pthread_mutex_destroy(&mtx);
    }

    if (!gb.have) return 0.0;
    m->copy_state(cfg_c, &gb.state, out_best);
    double adjust = 0.0;
    m->decrypt(ctx, cfg_c, &gb.state, out_decrypted, &adjust);
    return gb.score;
}

double run_solver(const CipherModel *m, SolverCtx *ctx) {

    static SolverConfig configs[MAX_SOLVER_CONFIGS];
    static SolverState best_state, cand_state;
    static int best_decrypted[MAX_CIPHER_LENGTH], cand_decrypted[MAX_CIPHER_LENGTH];
    static int sweep_best_decrypted[MAX_CIPHER_LENGTH];
    // Optimal-cycleword histogram scratch, owned by the engine (single-threaded).
    static int hist_by_col[MAX_CYCLEWORD_LEN * ALPHABET_SIZE];
    ctx->hist_by_col = hist_by_col;

    int nconf = m->enumerate_configs(ctx, configs, MAX_SOLVER_CONFIGS);
    if (nconf <= 0) return 0.0;

    double best_score = 0.0;
    bool have_best = false;
    SolverConfig best_cfg;

    for (int c = 0; c < nconf; c++) {
        SolverConfig *cc = &configs[c];

        if (m->needs_hist) engine_build_hist(ctx, cc->period);

        double sc;
        if (m->key_len && m->key_len(ctx, cc) == 0) {
            // SWEEP cell: the config carries no searched key, so a single
            // decrypt+score fully evaluates it. But colossus is a stochastic
            // hill climber WITH RESTARTS, and -nrestarts must never be inert
            // for any model, so the requested restart budget drives this branch
            // too: each restart re-seeds and re-evaluates, and the best sample
            // is kept. (route/railfence seed deterministically, so their
            // complete sweep is re-confirmed per restart; a model with a
            // stochastic seed would draw an independent sample each time.)
            int nrs = ctx->cfg->n_restarts;
            if (nrs < 1) nrs = 1;
            sc = 0.0;
            for (int rs = 0; rs < nrs; rs++) {
                double adjust = 0.0;
                m->seed(ctx, cc, &cand_state);
                m->decrypt(ctx, cc, &cand_state, cand_decrypted, &adjust);
                double s = engine_score(ctx, cand_decrypted, adjust);
                if (rs == 0 || s > sc) {
                    sc = s;
                    vec_copy(cand_decrypted, sweep_best_decrypted, ctx->cipher_len);
                }
            }
            vec_copy(sweep_best_decrypted, cand_decrypted, ctx->cipher_len);
        } else {
            sc = run_one_config(m, ctx, cc, &cand_state, cand_decrypted);
        }

        if (!have_best || sc > best_score) {
            best_score = sc; have_best = true;
            best_cfg = *cc;
            m->copy_state(cc, &cand_state, &best_state);
            vec_copy(cand_decrypted, best_decrypted, ctx->cipher_len);
        }
    }

    if (!have_best) return 0.0;

    m->report(ctx, &best_cfg, &best_state, best_score, best_decrypted);
    return best_score;
}

// Assemble the invariant SolverCtx the engine and every model hook read from.
// hist_by_col is left NULL here; run_solver() points it at its own scratch buffer.
SolverCtx make_solver_ctx(ColossusConfig *cfg, SharedData *shared, char *cribtext,
    int cipher[], int cipher_len, int crib_indices[], int crib_positions[], int n_cribs) {

    SolverCtx ctx;
    ctx.cfg = cfg;
    ctx.shared = shared;
    ctx.cipher = cipher;
    ctx.cipher_len = cipher_len;
    ctx.crib_indices = crib_indices;
    ctx.crib_positions = crib_positions;
    ctx.n_cribs = n_cribs;
    ctx.cribtext = cribtext;
    ctx.ngram_data = shared->ngram_data;
    ctx.hist_by_col = NULL;
    ctx.model_scratch = NULL;
    ctx.result = NULL;
    return ctx;
}


// Tuned per-cipher-type search schedules. See SearchDefaults (colossus.h). Only types
// whose ideal schedule differs materially from the init_config globals appear here;
// every other type is absent and so keeps the global defaults bit-for-bit.
//
//   Playfair: the score is a MEAN log-probability (so the natural temperature lives on
//   a much smaller scale than the polyalphabetic/transposition reward score), and the
//   digraph landscape is riddled with local optima. The profile below (single-letter-
//   swap-dominated anneal at inittemp 0.08, several backtracking restarts) reliably
//   recovers ~600+ character ciphers; below that Playfair is genuinely near the limit
//   of a quadgram attack (see tests/test_playfair_solver.c).
//   Bifid: the same square-anneal landscape as Playfair, but several candidate periods
//   are each annealed (run_solver keeps the global best by n-gram score), so the per-
//   period budget is smaller than Playfair's single-config budget to keep the whole
//   solve in the same ballpark. Same small-scale temperature (mean log-probability).
//   Trifid: the same fractionation-anneal as Bifid but over a 27-cell 3x3x3 cube (a
//   larger permutation space than Bifid's 25-cell square), so it gets a larger per-
//   period budget; otherwise identical small-scale temperature and per-period scheme.
//   Hill: the state is a k x k decryption matrix (every entry 0..25) hill-climbed /
//   annealed with the same mean-log-probability fitness, so it shares the small-scale
//   temperature. The matrices are small (k=2 is only 26^4 keys) and greedy climbs
//   converge fast on a rugged landscape, so -- unlike the fractionation types -- the
//   lever is RESTARTS, not iterations: the profile is many short restarts (250x8000),
//   run once per swept block size (k = 2..5).
static const SearchDefaults g_search_defaults[] = {
    { .cipher_type = PLAYFAIR, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 400000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = BIFID, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 4, .a_n_hill_climbs = 200000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 20, .s_n_hill_climbs = 200000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = TRIFID, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 300000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 24, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Digrafid: two KEYED 27-symbol alphabets (keyword + ascending tail) annealed per swept
    // period -- a structured search over short keywords, not the free 54-cell permutation (see
    // digrafid_solver.c). The keyed-alphabet prior collapses the keyspace by orders of
    // magnitude, so RESTARTS are the lever (each samples / refines a keyword length, and short
    // ciphers need many basins tried), and the keyword moves are coarser than a cell swap so
    // they want a HIGHER temperature -- hence MANY restarts at a warm 0.30 (48x150000) rather
    // than a few cold long climbs. This recovers reliably from ~300 letters (vs the old free-
    // permutation square break's ~700-800 cliff -- the keyed-alphabet prior is the whole win),
    // tuned against test_digrafid_solver. The budget is per period, so a blind sweep multiplies it.
    { .cipher_type = DIGRAFID, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 48, .a_n_hill_climbs = 150000,
      .a_init_temp = 0.30, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 200, .s_n_hill_climbs = 150000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // CM Bifid: Bifid fractionation over TWO keyed 5x5 squares (50 cells, double Bifid's 25),
    // searched as a JOINT two-square anneal (no square-independent decoupling reward exists --
    // both squares are entangled in the n-gram fitness; see cm_bifid_solver.c), so it inherits
    // Two-Square's larger/rougher two-square landscape but, like Bifid, has the period SWEPT --
    // the budget is PER period, so a blind sweep multiplies it. The long per-restart climb is the
    // critical lever (more restarts with shorter climbs does WORSE near the cliff), hence
    // 8x400000 not many-short-restarts. Recovers reliably from ~480 letters at an ODD period
    // (~400 cliff); EVEN periods are a documented ciphertext-only DEGENERACY -- the rows-then-cols
    // re-pairing splits into pure-row / pure-col output pairs, leaving a transpose-like square
    // ambiguity no budget escapes (~noise floor), so the schedule is tuned/asserted on odd P.
    // Tuned against test_cm_bifid_solver.
    { .cipher_type = CM_BIFID, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 400000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 20, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = HILL, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 250, .a_n_hill_climbs = 8000,
      .a_init_temp = 0.10, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.25,
      .s_n_restarts = 250, .s_n_hill_climbs = 8000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Phillips (and its column / row-column variants): the same 5x5-square anneal as
    // Playfair (one config, the base grid is the only unknown), so it shares Playfair's
    // small-scale temperature and backtracking. But Phillips is MONOGRAPHIC (every letter
    // is independently substituted, period 40), so it carries more signal per character
    // than digraphic Playfair and recovers reliably from ~200 characters at a leaner
    // budget -- 4x250000 lands a 760-char solve at ~100% in ~16s (see
    // tests/test_phillips_solver.c). Same profile for all three variants.
    { .cipher_type = PHILLIPS, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 4, .a_n_hill_climbs = 250000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 20, .s_n_hill_climbs = 250000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = PHILLIPS_C, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 4, .a_n_hill_climbs = 250000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 20, .s_n_hill_climbs = 250000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = PHILLIPS_RC, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 4, .a_n_hill_climbs = 250000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 20, .s_n_hill_climbs = 250000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Two-Square / Four-Square: the same digraphic square-anneal as Playfair, but the state
    // is a PAIR of 5x5 squares (50 cells, double Playfair's 25), so the landscape is larger
    // and rougher and the budget is bigger (more restarts and longer climbs). One config
    // (no period to estimate). Same small-scale temperature (mean log-probability) and
    // backtracking. Both two-square arrangements share the profile; Four-Square's two
    // independent keyed squares make it the hardest of the three, so it gets the most.
    { .cipher_type = TWO_SQUARE, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 600000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 500000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = TWO_SQUARE_V, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 600000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 500000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = FOUR_SQUARE, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 700000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 40, .s_n_hill_climbs = 600000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Tri-Square: the same digraphic square-anneal, but the state is THREE independent 5x5
    // squares (75 cells, the largest square state of the family) with no square-independent
    // decoupling reward -- jointly annealed, one config (no period to estimate). Despite the
    // bigger state it recovers MORE easily than Four-Square: the polyphonic c0/c2 cipher
    // letters spread the full alphabet across every position, giving the n-gram gradient sharp
    // signal, so it is reliable from ~500 plaintext letters (~750 cipher) at a modest budget
    // (a 300-400 cliff). This carries headroom over the proven-sufficient 8x400000. Same
    // small-scale temperature (mean log-probability) and backtracking as the square family.
    { .cipher_type = TRI_SQUARE, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 500000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 40, .s_n_hill_climbs = 400000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // ADFGX / ADFGVX: a COUPLED search -- a keyed Polybius square AND a keyed columnar
    // column order, jointly annealed (per swept column count K). The square anneal is
    // Bifid's; the column-order moves ride a structural IoC reward (independent of the
    // square) folded into score_adjust, which decouples the two halves. The landscape is
    // the hardest of the polygraphic family, so the budget is the largest. ADFGVX's
    // 36-cell square (vs ADFGX's 25) gets more. Same small-scale temperature (mean
    // log-probability) and backtracking as the other square types.
    { .cipher_type = ADFGX, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 600000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 40, .s_n_hill_climbs = 500000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = ADFGVX, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 800000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 50, .s_n_hill_climbs = 700000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Nihilist Substitution (and its no-carry / mod-100 variants): a COUPLED search -- a keyed
    // 5x5 square AND a periodic additive key, jointly annealed (one config per candidate
    // period). The square anneal is Bifid's; the additive-key moves ride a square-independent
    // "validity" reward (folded into score_adjust) that decouples the two halves, exactly like
    // ADFGVX's IoC term. Same small-scale temperature (mean log-probability) and backtracking
    // as the other square types; an 8x300000 budget per period (between Bifid and ADFGX -- the
    // additive adds a coupled lane but the square is only 25 cells). All three conventions
    // share the profile.
    { .cipher_type = NIHILIST_SUB, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 300000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = NIHILIST_SUB_NC, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 300000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = NIHILIST_SUB_M100, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 300000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Gromark / Periodic Gromark: a primer PRE-PASS (gromark_rank_primers) ranks the finite
    // primer space and emits one config per top-K primer; each config then anneals the keyed
    // 26-letter alphabet (a simple-substitution anneal -- easier than a Polybius square, so a
    // lean per-config budget suffices, and the pre-pass warm-starts the right primer's sigma).
    // Periodic also anneals the P group offsets jointly, so it gets a little more. Same
    // small-scale temperature (mean log-probability) as the other substitution types.
    { .cipher_type = GROMARK, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 3, .a_n_hill_climbs = 120000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 12, .s_n_hill_climbs = 120000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = GROMARK_PERIODIC, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 4, .a_n_hill_climbs = 160000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 16, .s_n_hill_climbs = 160000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Nicodemus (+ Variant / Beaufort): a per-(P,H) COLUMN-ORDER anneal (the P shifts are
    // derived deterministically per order, so the climbed state is just a short permutation of
    // length P). One config per (period P, block height H), so a lean per-config budget of MANY
    // short restarts is repeated across the sweep -- the small permutation climbs converge fast
    // and the landscape is rugged, so restarts (independent draws) are the robustness lever, not
    // climbs (tuned in tests/test_nicodemus_solver.c). Same small-scale temperature (mean
    // log-probability) as the other -logprob substitution types.
    { .cipher_type = NICODEMUS, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 20000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 16, .s_n_hill_climbs = 20000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = NICODEMUS_VARIANT, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 20000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 16, .s_n_hill_climbs = 20000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = NICODEMUS_BEAUFORT, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 20000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 16, .s_n_hill_climbs = 20000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },

    // Bazeries: the climbed state is the key NUMBER's decimal digits (one config per digit
    // count D in 1..6), a tiny rugged < 10^6 keyspace, so RESTARTS are the robustness lever
    // (each restart reseeds a fresh random number; the square-quality monogram reward then
    // pulls the climb toward the right square). Many restarts x modest climbs, per D config.
    // Tuned in tests/test_bazeries_solver.c.
    { .cipher_type = BAZERIES, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 40, .a_n_hill_climbs = 20000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 40, .s_n_hill_climbs = 20000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },

    // Portax: the climbed state is the P per-column Porta shifts (0..12), a tiny per-period key
    // that the monogram-fit warm start gets mostly right on seed; the anneal/n-gram pass only
    // needs to correct a few columns. One config per swept period P, so MANY short restarts (the
    // robustness lever) x modest climbs. Same small-scale (mean log-probability) temperature as
    // the other Porta-family / -logprob types. Tuned in tests/test_portax_solver.c.
    { .cipher_type = PORTAX, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 20000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 12, .s_n_hill_climbs = 20000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },

    // Progressive Key (Vigenere / Variant / Beaufort base). The climbed state is the P per-column
    // base shifts (0..25); the per-column monogram-fit warm start gets most of them right on seed,
    // so the anneal only needs to correct a few columns. Many (P, prog) configs are enumerated
    // (period brute-forced x progression 0..25, since IoC fails through the drift), so each config
    // gets a LEAN budget: a few restarts x modest climbs. The reward-only quadgram table suffices
    // (Vigenere family), so the same small-scale temperature as the other -optimalcycle types.
    // Tuned in tests/test_progkey_solver.c.
    { .cipher_type = PROGKEY, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 3, .a_n_hill_climbs = 2500,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 3, .s_n_hill_climbs = 2500,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = PROGKEY_VAR, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 3, .a_n_hill_climbs = 2500,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 3, .s_n_hill_climbs = 2500,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = PROGKEY_BEAU, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 3, .a_n_hill_climbs = 2500,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 3, .s_n_hill_climbs = 2500,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },

    // Slidefair (periodic digraphic Vigenere / Variant / Beaufort). Like Portax, the climbed state
    // is the P per-column key letters (0..25) that the monogram-fit warm start gets mostly right on
    // seed; the anneal/n-gram pass only corrects a few columns -- recovery is so strong (100% from
    // ~50 letters in tests/test_slidefair_solver.c) that a LEAN budget suffices: a few short restarts
    // (the robustness lever) x modest climbs, one config per swept period P. Reward-only quadgram
    // table (Vigenere family, no -logprob), same small-scale temperature. Tuned in that test.
    { .cipher_type = SLIDEFAIR, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 10000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 8, .s_n_hill_climbs = 10000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = SLIDEFAIR_VAR, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 10000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 8, .s_n_hill_climbs = 10000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = SLIDEFAIR_BEAU, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 8, .a_n_hill_climbs = 10000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 8, .s_n_hill_climbs = 10000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },

    // Seriated Playfair: a single 5x5 keyed square (no per-column decoupling), so the
    // climbed state is the whole grid and each config is a FULL Playfair-scale grid anneal
    // -- but enumerated once per swept seriation period P, so the budget is per-P and the
    // blind sweep multiplies the cost (document it). Playfair's small-scale (mean log-
    // probability) temperature; effectively needs -logprob. Uses Playfair's proven
    // 6x400000 -- the 400000-climb cooling schedule is what reliably reaches the optimum
    // (300000 leaves it seed-sensitive). Tuned in tests/test_seriated_playfair_solver.c.
    { .cipher_type = SERIATED_PLAYFAIR, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 400000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 20, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },

    // Interrupted Key (periodic Vigenere / Variant / Beaufort keyword reset at break points). The
    // climbed state is just the P per-column key letters (0..25); a Gromark-style pre-pass derives
    // the monogram-best keyword per (period, interruptor) and keeps the top-K, so the seed is close
    // and only a few short restarts x modest climbs are needed to correct a column or two -- the
    // same lean profile as Progressive Key / Slidefair. One config per kept pre-pass candidate (CT
    // recovers on seed like a Vigenere; PT and the JOINT/breaks random modes need the anneal).
    // Reward-only quadgram table (Vigenere family, no -logprob), same small-scale temperature.
    // Tuned in tests/test_intkey_solver.c.
    { .cipher_type = INTERRUPTED_KEY, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 8000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 6, .s_n_hill_climbs = 8000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = INTERRUPTED_KEY_VAR, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 8000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 6, .s_n_hill_climbs = 8000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = INTERRUPTED_KEY_BEAU, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 8000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 6, .s_n_hill_climbs = 8000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Condi: a free-permutation sigma anneal per swept starter (26 configs). NOTE the plaintext-
    // feedback cascade makes the true sigma an ISOLATED NEEDLE (no basin: even the best single-swap
    // neighbour scores ~0.9 below the true key, the mean neighbour is near the random floor), so no
    // local-search budget cracks it blind -- this is a documented structural limitation, not a
    // tuning target (see tests/test_condi_solver.c). The budget is therefore kept MODEST (a bounded
    // honest attempt that terminates), not large.
    { .cipher_type = CONDI, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 6, .a_n_hill_climbs = 60000,
      .a_init_temp = 0.08, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 6, .s_n_hill_climbs = 60000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Fractionated Morse: a single KEYED 26-letter alphabet (keyword + ascending tail) annealed
    // as a keyed-alphabet SEQUENCE (fracmorse_move_seq), not a free 26! permutation -- the same
    // structured search as Digrafid, so the coarse keyword moves want a WARM temperature and short
    // ACA ciphers want many basins tried: MANY warm restarts (24x200000, inittemp 0.30). There is
    // NO period (one config), so unlike Digrafid the budget is not multiplied by a sweep. The tiled
    // decode + Morse-validity reward make it effectively need -logprob. Tuned against
    // test_fracmorse_solver.
    { .cipher_type = FRAC_MORSE, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 120000,
      .a_init_temp = 0.30, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 120, .s_n_hill_climbs = 120000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Straddling Checkerboard: the per-config SA MINI-SOLVE pre-pass (in the solver) does the
    // heavy lifting and warm-starts each kept indicator-pair config near the solution, so the
    // engine anneal is a short warm polish (a few restarts). inittemp 0.30 like the other
    // keyed digit/alphabet types; global-best tracking preserves the warm map through it.
    { .cipher_type = STRADDLING_CHECKERBOARD, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 4, .a_n_hill_climbs = 20000,
      .a_init_temp = 0.30, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 30, .s_n_hill_climbs = 40000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Ragbaby: a single keyed-alphabet anneal (24-letter KA; keyword prefix + ordered tail, the
    // ragbaby_move_seq twin of fracmorse). The per-letter shift is KNOWN, so the KA is heavily
    // constrained -- the coarse keyword moves want a warm temperature and RESTARTS are the lever
    // (no period, a single config). Same profile as Fractionated Morse; tuned vs the solver test.
    { .cipher_type = RAGBABY, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 120000,
      .a_init_temp = 0.30, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 120, .s_n_hill_climbs = 120000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Monome-Dinome: the Straddling Checkerboard's letter-only cousin. The solver PRE-FILTERS
    // the 45 indicator pairs by structural token validity and keeps only the fully-valid set
    // (1-6 configs, always incl. the true pair), so the engine anneals very few configs -- and
    // can therefore afford a bigger per-config budget. Near the ~120-letter floor a single
    // long climb occasionally lands in a near-miss basin, so RESTARTS are the lever (20x30000):
    // more basins tried removes the seed-dependent dips to ~80-85%, landing a stable ~98%. The
    // free code->letter substitution over a mis-segmented stream games a quadgram score, so
    // Monome-Dinome effectively NEEDS -logprob with QUINTGRAMS; recovery is reliable from ~120
    // letters (the ACA 60-100 low end is below the blind floor). Tuned against
    // test_monome_dinome_solver.
    { .cipher_type = MONOME_DINOME, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 20, .a_n_hill_climbs = 30000,
      .a_init_temp = 0.30, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 60, .s_n_hill_climbs = 40000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Tridigital: anneal a 26-letter -> 9-group partition, scored by an inner beam-Viterbi decode
    // (each eval decodes the whole stream, so climbs are pricier than the other digit types -> a
    // smaller climb budget, but restart-heavy). Each of the kept separator configs is mini-solve
    // warm-started, so the engine only polishes. Best with -logprob + quintgrams + a -dictionary
    // (the dense polyphonic decode games raw n-gram; coverage selection needs the dictionary).
    { .cipher_type = TRIDIGITAL, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 8000,
      .a_init_temp = 0.30, .a_min_temp = 0.001, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 40, .s_n_hill_climbs = 12000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Aristocrat / Patristocrat: a simple 26-letter monoalphabetic substitution, climbed as a free
    // 26-permutation by n-gram score with the homophonic-style INCREMENTAL fast path (each swap is
    // scored as a delta over the touched windows), so iterations are very cheap and the budget can
    // be large. A frequency rank-match warm start puts most cells right; the anneal corrects the
    // rest. Aristocrats are short (~80-150 letters) so RESTARTS + a warm temperature are the levers.
    // Best with -logprob (a free 26-sub is weaker signal than Ragbaby's known-shift alphabet). Both
    // variants share the schedule (identical search; only the report differs). Tuned vs the solver test.
    { .cipher_type = ARISTOCRAT, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 200000,
      .a_init_temp = 0.15, .a_min_temp = 0.0005, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 60, .s_n_hill_climbs = 200000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    { .cipher_type = PATRISTOCRAT, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 200000,
      .a_init_temp = 0.15, .a_min_temp = 0.0005, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 60, .s_n_hill_climbs = 200000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Checkerboard: once the per-axis label grouping is fixed, the search is a free 25-code -> 25-
    // letter bijection over the merged codes -- an ARISTOCRAT over 25 symbols -- so it shares the
    // aristocrat profile (warm freq-rank seed + the homophonic incremental fast path; every move is
    // a two-code swap, aristocrat granularity, NOT digrafid's coarse keyword move -> temp 0.15, not
    // 0.30). Restarts are the lever: the simple case is short (60-90 ACA letters) and the complex
    // case is strongly bimodal (finds the pairing basin or lands in garbage). Best with -logprob.
    { .cipher_type = CHECKERBOARD, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 12, .a_n_hill_climbs = 200000,
      .a_init_temp = 0.15, .a_min_temp = 0.0005, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 60, .s_n_hill_climbs = 200000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Grandpre: decoding is code -> letter, so the search is a HOMOPHONIC map over <= N^2
    // numeric codes -> 26 letters (single-symbol reassignment moves on the incremental fast
    // path). Shares the homophonic/aristocrat annealed-square profile; restarts are the lever
    // (undersampled ~64-code map over the ACA range is bimodal). Best with -logprob.
    { .cipher_type = GRANDPRE, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 16, .a_n_hill_climbs = 200000,
      .a_init_temp = 0.15, .a_min_temp = 0.0005, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 60, .s_n_hill_climbs = 200000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
    // Syllabary: a 100-token composite-map substitution with a length-changing (tiled) decode.
    // Two-code swap moves on the generic path; more restarts (the 100-token bijection over ~150
    // codes is severely undersampled and strongly bimodal). Same aristocrat temp scale. -logprob.
    { .cipher_type = SYLLABARY, .default_shape = SHAPE_ANNEAL,
      .a_n_restarts = 24, .a_n_hill_climbs = 300000,
      .a_init_temp = 0.15, .a_min_temp = 0.0005, .a_cooling_rate = 0.0,
      .a_backtracking_probability = 0.30,
      .s_n_restarts = 80, .s_n_hill_climbs = 300000,
      .s_slip_probability = 0.0005, .s_backtracking_probability = 0.20 },
};

bool apply_cipher_defaults(ColossusConfig *cfg, bool announce) {
    const SearchDefaults *d = NULL;
    for (size_t i = 0; i < sizeof(g_search_defaults) / sizeof(g_search_defaults[0]); i++)
        if (g_search_defaults[i].cipher_type == cfg->cipher_type) { d = &g_search_defaults[i]; break; }
    if (d == NULL) return false;

    // Effective shape: an explicit -method wins, else the type's own default shape.
    SearchShape shape = d->default_shape;
    if (cfg->method == METHOD_SHOTGUN) shape = SHAPE_SHOTGUN;
    else if (cfg->method == METHOD_ANNEAL) shape = SHAPE_ANNEAL;
    else if (cfg->method == METHOD_PSO) shape = SHAPE_PSO;

    if (shape == SHAPE_PSO) {
        // PSO reuses n_restarts (swarm relaunches) / n_hill_climbs (iterations) plus
        // its own swarm knobs. p_n_particles == 0 => this type has no tuned PSO
        // profile, so keep the init_config globals (mirrors a type with no entry).
        if (d->p_n_particles > 0) {
            cfg->n_restarts = d->p_n_restarts;
            cfg->n_hill_climbs = d->p_n_hill_climbs;
            cfg->n_particles = d->p_n_particles;
            cfg->inertia = d->p_inertia;
            cfg->cognitive = d->p_cognitive;
            cfg->social = d->p_social;
            cfg->refine_steps = d->p_refine_steps;
            if (announce)
                printf("-type defaults: pso schedule %dx%d swarm %d (inertia %.2f, cog %.2f, soc %.2f, refine %d)\n",
                    d->p_n_restarts, d->p_n_hill_climbs, d->p_n_particles,
                    d->p_inertia, d->p_cognitive, d->p_social, d->p_refine_steps);
        } else if (announce) {
            printf("-type defaults: pso (global schedule %dx%d swarm %d)\n",
                cfg->n_restarts, cfg->n_hill_climbs, cfg->n_particles);
        }
    } else if (shape == SHAPE_SHOTGUN) {
        cfg->n_restarts = d->s_n_restarts;
        cfg->n_hill_climbs = d->s_n_hill_climbs;
        cfg->slip_probability = d->s_slip_probability;
        cfg->backtracking_probability = d->s_backtracking_probability;
        if (announce)
            printf("-type defaults: shotgun schedule %dx%d (slipprob %.4f, backtrack %.2f)\n",
                d->s_n_restarts, d->s_n_hill_climbs, d->s_slip_probability,
                d->s_backtracking_probability);
    } else {
        cfg->n_restarts = d->a_n_restarts;
        cfg->n_hill_climbs = d->a_n_hill_climbs;
        cfg->init_temp = d->a_init_temp;
        cfg->min_temp = d->a_min_temp;
        cfg->cooling_rate = d->a_cooling_rate;
        cfg->backtracking_probability = d->a_backtracking_probability;
        if (announce)
            printf("-type defaults: anneal schedule %dx%d (inittemp %.3f, backtrack %.2f)\n",
                d->a_n_restarts, d->a_n_hill_climbs, d->a_init_temp,
                d->a_backtracking_probability);
    }
    return true;
}




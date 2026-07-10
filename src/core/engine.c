#include "engine.h"
#include "scoring.h"




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

// Incremental variant of run_one_config (same restart / annealing / backtrack /
// best-tracking skeleton) for models that supply the score_neighbor/commit/sync
// hooks. The current state's decryption (cur_dec) is kept live; each neighbour is
// scored as a delta and only committed into cur_dec + the model caches on accept,
// so the per-iteration cost is O(positions the move touched) rather than O(N).
static double run_one_config_incremental(const CipherModel *m, SolverCtx *ctx,
                                         const SolverConfig *cfg_c,
                                         SolverState *out_best, int *out_decrypted) {

    static SolverState cur, loc, best;
    static int cur_dec[MAX_CIPHER_LENGTH], best_dec[MAX_CIPHER_LENGTH];
    ColossusConfig *cfg = ctx->cfg;

    double best_score = 0.0, cur_score, loc_score, adjust;
    bool have_best = false;
    bool force_primary = true;

    SearchShape shape = m->shape;
    if (cfg->method == METHOD_SHOTGUN) shape = SHAPE_SHOTGUN;
    else if (cfg->method == METHOD_ANNEAL) shape = SHAPE_ANNEAL;
    bool done = false;

    EngineStats st;
    memset(&st, 0, sizeof st);
    st.start_time = clock();

    double temp_start = cfg->init_temp;
    double cooling;
    if (cfg->cooling_rate > 0.0) {
        cooling = cfg->cooling_rate;
    } else {
        cooling = 1.0;
        if (cfg->n_hill_climbs > 1)
            cooling = pow(cfg->min_temp / temp_start, 1.0 / (double)(cfg->n_hill_climbs - 1));
    }

    for (int rs = 0; rs < cfg->n_restarts && !done; rs++) {
        st.n_restarts = rs;

        if (have_best && frand() < cfg->backtracking_probability) {
            m->copy_state(cfg_c, &best, &cur);
            vec_copy(best_dec, cur_dec, ctx->cipher_len);
            cur_score = best_score;
            st.n_backtracks++;
        } else {
            m->seed(ctx, cfg_c, &cur);
            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &cur, cur_dec, &adjust);
            cur_score = engine_score(ctx, cur_dec, adjust);
        }
        // (Re)build the model caches so they describe the current decryption.
        m->sync_caches(ctx, cfg_c, cur_dec);

        force_primary = true;

        double temp = temp_start;
        for (int it = 0; it < cfg->n_hill_climbs && !done; it++) {
            st.n_iterations++;

            m->copy_state(cfg_c, &cur, &loc);
            m->perturb(ctx, cfg_c, &loc, &force_primary);

            loc_score = m->score_neighbor(ctx, cfg_c, &cur, &loc, cur_dec, cur_score);

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
                m->commit_neighbor(ctx, cfg_c, cur_dec);   // advance cur_dec + caches
                m->copy_state(cfg_c, &loc, &cur);
                cur_score = loc_score;
            }
            temp *= cooling;

            if (!have_best || cur_score > best_score) {
                best_score = cur_score; have_best = true;
                m->copy_state(cfg_c, &cur, &best);
                vec_copy(cur_dec, best_dec, ctx->cipher_len);
                if (cfg->verbose && m->report_verbose)
                    m->report_verbose(ctx, cfg_c, &best, best_score, cur_dec, &st);
                if (shape == SHAPE_DETERMINISTIC) done = true;
            }
        }
    }

    m->copy_state(cfg_c, &best, out_best);
    vec_copy(best_dec, out_decrypted, ctx->cipher_len);
    return best_score;
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
static double run_one_config_pso(const CipherModel *m, SolverCtx *ctx,
                                 const SolverConfig *cfg_c,
                                 SolverState *out_best, int *out_decrypted) {

    static SolverState particle[MAX_PSO_PARTICLES];
    static SolverState pbest[MAX_PSO_PARTICLES];
    static double      pbest_score[MAX_PSO_PARTICLES];
    static SolverState gbest, trial, loc;
    static int decrypted[MAX_CIPHER_LENGTH], gbest_dec[MAX_CIPHER_LENGTH];
    ColossusConfig *cfg = ctx->cfg;

    int np = cfg->n_particles;
    if (np < 1) np = 1;
    if (np > MAX_PSO_PARTICLES) np = MAX_PSO_PARTICLES;
    int refine = cfg->refine_steps < 0 ? 0 : cfg->refine_steps;

    double gbest_score = 0.0;
    bool have_gbest = false;
    bool force_primary = false;
    double adjust;

    EngineStats st;
    memset(&st, 0, sizeof st);
    st.start_time = clock();

    for (int rs = 0; rs < cfg->n_restarts; rs++) {
        st.n_restarts = rs;

        // (Re)seed the swarm and initialise personal / global bests.
        for (int p = 0; p < np; p++) {
            m->seed(ctx, cfg_c, &particle[p]);
            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &particle[p], decrypted, &adjust);
            double sc = engine_score(ctx, decrypted, adjust);
            m->copy_state(cfg_c, &particle[p], &pbest[p]);
            pbest_score[p] = sc;
            if (!have_gbest || sc > gbest_score) {
                gbest_score = sc; have_gbest = true;
                m->copy_state(cfg_c, &particle[p], &gbest);
                vec_copy(decrypted, gbest_dec, ctx->cipher_len);
            }
        }

        for (int it = 0; it < cfg->n_hill_climbs; it++) {
            st.n_iterations++;
            for (int p = 0; p < np; p++) {
                // 1. inertia: a little random momentum (exploration).
                int n_inertia = (int)(cfg->inertia * (1.0 + frand()) + 0.5);
                for (int k = 0; k < n_inertia; k++)
                    m->perturb(ctx, cfg_c, &particle[p], &force_primary);

                // 2. cognitive pull toward this particle's personal best.
                int n_cog = (int)(cfg->cognitive * frand() *
                                  state_distance(cfg_c, &particle[p], &pbest[p]) + 0.5);
                if (n_cog > 0)
                    pull_toward(m, ctx, cfg_c, &particle[p], &pbest[p], n_cog, &trial);

                // 3. social pull toward the global best.
                int n_soc = (int)(cfg->social * frand() *
                                  state_distance(cfg_c, &particle[p], &gbest) + 0.5);
                if (n_soc > 0)
                    pull_toward(m, ctx, cfg_c, &particle[p], &gbest, n_soc, &trial);

                // 4. memetic local refinement: a short greedy hill-climb.
                adjust = 0.0;
                m->decrypt(ctx, cfg_c, &particle[p], decrypted, &adjust);
                double sc = engine_score(ctx, decrypted, adjust);
                for (int k = 0; k < refine; k++) {
                    m->copy_state(cfg_c, &particle[p], &loc);
                    m->perturb(ctx, cfg_c, &loc, &force_primary);
                    double adj2 = 0.0;
                    m->decrypt(ctx, cfg_c, &loc, decrypted, &adj2);
                    double ls = engine_score(ctx, decrypted, adj2);
                    if (ls > sc) { m->copy_state(cfg_c, &loc, &particle[p]); sc = ls; }
                }

                // 5. update personal and global bests.
                if (sc > pbest_score[p]) {
                    m->copy_state(cfg_c, &particle[p], &pbest[p]);
                    pbest_score[p] = sc;
                    if (sc > gbest_score) {
                        gbest_score = sc;
                        m->copy_state(cfg_c, &particle[p], &gbest);
                        adjust = 0.0;
                        m->decrypt(ctx, cfg_c, &gbest, gbest_dec, &adjust);
                        if (cfg->verbose && m->report_verbose)
                            m->report_verbose(ctx, cfg_c, &gbest, gbest_score, gbest_dec, &st);
                    }
                }
            }
        }
    }

    m->copy_state(cfg_c, &gbest, out_best);
    adjust = 0.0;
    m->decrypt(ctx, cfg_c, &gbest, out_decrypted, &adjust);
    return gbest_score;
}

// Hill-climb one outer config: shotgun restarts + per-iteration neighbour move,
// with SHOTGUN slip / ANNEAL Metropolis acceptance and best-state tracking. Writes
// the config's best state to *out_best and its decryption to out_decrypted, and
// returns the best score.
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

    static SolverState cur, loc, best;
    static int decrypted[MAX_CIPHER_LENGTH];
    ColossusConfig *cfg = ctx->cfg;

    double best_score = 0.0, cur_score, loc_score, adjust;
    bool have_best = false;
    bool force_primary = true;               // first iteration forces a primary (keyword) move

    // -method overrides the model's built-in shape on EVERY cipher type (the engine
    // is acceptance-strategy agnostic); METHOD_DEFAULT keeps the model's own shape.
    SearchShape shape = m->shape;
    if (cfg->method == METHOD_SHOTGUN) shape = SHAPE_SHOTGUN;
    else if (cfg->method == METHOD_ANNEAL) shape = SHAPE_ANNEAL;
    bool deterministic = (shape == SHAPE_DETERMINISTIC);
    bool done = false;

    EngineStats st;
    memset(&st, 0, sizeof st);
    st.start_time = clock();

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

    for (int rs = 0; rs < cfg->n_restarts && !done; rs++) {
        st.n_restarts = rs;

        if (have_best && frand() < cfg->backtracking_probability) {
            m->copy_state(cfg_c, &best, &cur);
            cur_score = best_score;
            st.n_backtracks++;
        } else {
            m->seed(ctx, cfg_c, &cur);
            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &cur, decrypted, &adjust);
            cur_score = engine_score(ctx, decrypted, adjust);
        }

        // Best is recorded only inside the hill-climb loop below (matching every
        // original climber), so have_best transitions to true on the first
        // iteration -- the restart-0 backtrack check then draws no RNG, keeping
        // the draw sequence identical to the per-cipher climbers this replaces.

        // force_primary is the per-restart "must perturb the primary lane" flag
        // (polyalpha's perturbate_keyword_p); the model reads and updates it.
        force_primary = true;

        double temp = temp_start;
        for (int it = 0; it < cfg->n_hill_climbs && !done; it++) {
            st.n_iterations++;

            m->copy_state(cfg_c, &cur, &loc);

            m->perturb(ctx, cfg_c, &loc, &force_primary);

            adjust = 0.0;
            m->decrypt(ctx, cfg_c, &loc, decrypted, &adjust);
            loc_score = engine_score(ctx, decrypted, adjust);

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
                m->copy_state(cfg_c, &loc, &cur);
                cur_score = loc_score;
            }
            temp *= cooling;

            if (!have_best || cur_score > best_score) {
                best_score = cur_score; have_best = true;
                m->copy_state(cfg_c, &cur, &best);
                if (cfg->verbose && m->report_verbose)
                    m->report_verbose(ctx, cfg_c, &best, best_score, decrypted, &st);
                if (deterministic) done = true;
            }
        }
    }

    m->copy_state(cfg_c, &best, out_best);
    adjust = 0.0;
    m->decrypt(ctx, cfg_c, &best, out_decrypted, &adjust);
    return best_score;
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




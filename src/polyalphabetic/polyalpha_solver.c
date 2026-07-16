#include "polyalpha_solver.h"
#include "engine.h"
#include "scoring.h"
#include "trans_common.h"

// ====================================================================
//  Polyalphabetic model (Vigenere / Quagmire I-IV / Beaufort / Porta /
//  Autokey*; cipher-agnostic engine)
// ====================================================================
//
// All 14 polyalphabetic types share one model and one uniform state (a pt and ct
// keyed alphabet plus the periodic cycleword). The per-type differences -- which
// alphabet is straight vs keyed, which is perturbed, the optimal-cycleword refine,
// the crib constraint, and -samekey coupling -- stay as the explicit
// switch(cipher_type) ladders moved verbatim from the old shotgun_hill_climber
// into the seed and perturb hooks (per CLAUDE.md: don't genericise the dispatch).

static inline bool poly_is_autokey(int t) {
    return (t >= AUTOKEY_0 && t <= AUTOKEY_4) || t == AUTOKEY_BEAU || t == AUTOKEY_PORTA;
}

// The polyalphabetic model carries a single verbose-only contradiction counter
// through ctx->model_scratch (a long*); solve_cipher owns the storage on its stack.

// Outer configs: cycleword length(s) x (pt_keyword_len j, ct_keyword_len k), with
// the per-type validity pruning and crib gate. Reproduces solve_cipher's old
// period-estimation + keyword-bound + (j,k) loop, printing the same diagnostics.
static int polyalpha_enumerate(const SolverCtx *ctx, SolverConfig *out, int cap) {
    ColossusConfig *cfg = ctx->cfg;
    int *cipher_indices = ctx->cipher;
    int cipher_len = ctx->cipher_len;
    int n_cribs = ctx->n_cribs;
    int *crib_indices = ctx->crib_indices, *crib_positions = ctx->crib_positions;

    int n_cycleword_lengths, cycleword_lengths[MAX_CIPHER_LENGTH];

    // --- cycleword / primer length setup ---
    if (cfg->cycleword_len_present) {
        n_cycleword_lengths = 1;
        cycleword_lengths[0] = cfg->cycleword_len;
    } else if ((cfg->cipher_type >= AUTOKEY_0 && cfg->cipher_type <= AUTOKEY_4) ||
        cfg->cipher_type == AUTOKEY_BEAU || cfg->cipher_type == AUTOKEY_PORTA ||
        cfg->transperoffset_present) {
        // Autokey (aperiodic) or polyalphabetic+transposition: IoC fails, brute-force.
        n_cycleword_lengths = 0;
        for (int len = 1; len <= cfg->max_cycleword_len; len++) {
            cycleword_lengths[n_cycleword_lengths++] = len;
        }
    } else {
        estimate_cycleword_lengths(cipher_indices, cipher_len, cfg->max_cycleword_len,
            cfg->n_sigma_threshold, cfg->ioc_threshold,
            &n_cycleword_lengths, cycleword_lengths, cfg->verbose);
        if (n_cycleword_lengths == 0) {
            if (cfg->verbose) printf("No periodicities found above threshold. Not running the hillclimber and exiting.\n");
            return 0;
        }
    }

    // --- keyword length bounds (cipher-type specific; mutates cfg as before) ---
    int min_kw = cfg->min_keyword_len;
    int pt_max = cfg->plaintext_max_keyword_len;
    int ct_max = cfg->ciphertext_max_keyword_len;

    if (cfg->cipher_type == VIGENERE || cfg->cipher_type == GRONSFELD ||
        cfg->cipher_type == BEAUFORT ||
        cfg->cipher_type == PORTA ||
        (cfg->cipher_type >= AUTOKEY_0 && cfg->cipher_type <= AUTOKEY_2) ||
        cfg->cipher_type == QUAGMIRE_1 || cfg->cipher_type == QUAGMIRE_2 ||
        cfg->cipher_type == AUTOKEY_BEAU || cfg->cipher_type == AUTOKEY_PORTA) {
        min_kw = 1;
    }

    if (cfg->cipher_type == VIGENERE || cfg->cipher_type == GRONSFELD ||
        cfg->cipher_type == AUTOKEY_0 ||
        cfg->cipher_type == AUTOKEY_BEAU || cfg->cipher_type == AUTOKEY_PORTA) {
        pt_max = 2; ct_max = 2;
        cfg->plaintext_keyword_len = 0;
        cfg->ciphertext_keyword_len = 0;
    } else if (cfg->cipher_type == BEAUFORT) {
        pt_max = 2;
        cfg->plaintext_keyword_len = 1;
    } else if (cfg->cipher_type == PORTA) {
        pt_max = 2; ct_max = 2;
        cfg->plaintext_keyword_len = 1;
        cfg->ciphertext_keyword_len = 1;
    } else if (cfg->cipher_type == QUAGMIRE_1 || cfg->cipher_type == AUTOKEY_1) {
        ct_max = 2;
        cfg->ciphertext_keyword_len = 1;
    } else if (cfg->cipher_type == QUAGMIRE_2 || cfg->cipher_type == AUTOKEY_2) {
        pt_max = 2;
        cfg->plaintext_keyword_len = 1;
    }

    // --- enumerate (cycleword_len, j, k) with per-type pruning + crib gate ---
    int n = 0;
    for (int i = 0; i < n_cycleword_lengths; i++) {
        if (cfg->verbose) printf("\ncycleword length = %d\n", cycleword_lengths[i]);
        for (int j = min(min_kw, cfg->plaintext_keyword_len); j < pt_max; j++) {
            for (int k = min(min_kw, cfg->ciphertext_keyword_len); k < ct_max; k++) {
                if (cfg->verbose) printf("\npt/ct keyword len = %d, %d\n", j, k);
                if (cfg->plaintext_keyword_len_present && j != cfg->plaintext_keyword_len) continue;
                if (cfg->ciphertext_keyword_len_present && k != cfg->ciphertext_keyword_len) continue;

                if (cfg->cipher_type == QUAGMIRE_3 && j != k) continue;
                if (cfg->cipher_type == BEAUFORT && ! (j == 1 && k == 1)) continue;
                if (cfg->cipher_type == VIGENERE && ! (j == 1 && k == 1)) continue;
                if (cfg->cipher_type == GRONSFELD && ! (j == 1 && k == 1)) continue;
                if (cfg->cipher_type == PORTA && ! (j == 1 && k == 1)) continue;
                if (cfg->cipher_type == AUTOKEY_0 && ! (j == 1 && k == 1)) continue;
                if (cfg->cipher_type == AUTOKEY_1 && k != 1) continue;
                if (cfg->cipher_type == AUTOKEY_2 && j != 1) continue;
                if (cfg->cipher_type == AUTOKEY_3 && j != k) continue;
                if (cfg->cipher_type == AUTOKEY_BEAU && ! (j == 1 && k == 1)) continue;
                if (cfg->cipher_type == AUTOKEY_PORTA && ! (j == 1 && k == 1)) continue;

                if (cfg->cipher_type != AUTOKEY_0 && cfg->cipher_type != AUTOKEY_1 &&
                    cfg->cipher_type != AUTOKEY_2 && cfg->cipher_type != AUTOKEY_3 &&
                    cfg->cipher_type != AUTOKEY_4 && cfg->cipher_type != AUTOKEY_BEAU &&
                    cfg->cipher_type != AUTOKEY_PORTA) {
                    if (!cribs_satisfied_p(cfg, cipher_indices, cipher_len, crib_indices,
                        crib_positions, n_cribs, cycleword_lengths[i], cfg->verbose)) {
                        #if CRIB_CHECK
                        continue;
                        #endif
                    }
                }

                if (n < cap) {
                    out[n].period = cycleword_lengths[i];
                    out[n].j = j; out[n].k = k;
                    out[n].aux[0] = 0; out[n].aux[1] = 0;
                    n++;
                }
            }
        }
    }

    if (n == 0) {
        printf("\n\nERROR: No valid configuration found. Check your length constraints.\n");
        printf("Debug: Type=%d, PT_Len=%d, CT_Len=%d\n",
               cfg->cipher_type, cfg->plaintext_keyword_len, cfg->ciphertext_keyword_len);
    }
    return n;
}

static void polyalpha_seed(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st) {
    ColossusConfig *cfg = ctx->cfg;
    int *current_plaintext_keyword_state = st->pt_keyword;
    int *current_ciphertext_keyword_state = st->ct_keyword;
    int *current_cycleword_state = st->cycleword;
    int cycleword_len = cc->period;
    int plaintext_keyword_len = cc->j, ciphertext_keyword_len = cc->k;
    bool is_autokey = poly_is_autokey(cfg->cipher_type);
    int i;

    switch (cfg->cipher_type) {
        case VIGENERE:
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case GRONSFELD:
            // Gronsfeld = Vigenere with straight alphabets, but the shifts are
            // numeric key digits, so the cycleword domain is 0..9 (not 0..25).
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, GRONSFELD_DIGITS, cycleword_len);
            break ;
        case QUAGMIRE_1:
            if (cfg->user_plaintext_keyword_present) {
                make_keyed_alphabet(cfg->user_plaintext_keyword, current_plaintext_keyword_state);
            } else {
                random_keyword(current_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
            }
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case QUAGMIRE_2:
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            if (cfg->user_ciphertext_keyword_present) {
                make_keyed_alphabet(cfg->user_ciphertext_keyword, current_ciphertext_keyword_state);
            } else {
                random_keyword(current_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
            }
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case QUAGMIRE_3:
            if (cfg->user_plaintext_keyword_present) {
                make_keyed_alphabet(cfg->user_plaintext_keyword, current_plaintext_keyword_state);
            } else if (cfg->user_ciphertext_keyword_present) {
                make_keyed_alphabet(cfg->user_ciphertext_keyword, current_plaintext_keyword_state);
            } else {
                random_keyword(current_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
            }
            vec_copy(current_plaintext_keyword_state, current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case QUAGMIRE_4:
            if (cfg->user_plaintext_keyword_present) {
                make_keyed_alphabet(cfg->user_plaintext_keyword, current_plaintext_keyword_state);
            } else {
                random_keyword(current_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
            }
            if (cfg->user_ciphertext_keyword_present) {
                make_keyed_alphabet(cfg->user_ciphertext_keyword, current_ciphertext_keyword_state);
            } else {
                random_keyword(current_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
            }
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case BEAUFORT:
            plaintext_keyword_len = g_alpha;
            ciphertext_keyword_len = g_alpha;
            for (i = 0; i < g_alpha; i++) current_plaintext_keyword_state[i] = i;
            vec_copy(current_plaintext_keyword_state, current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case PORTA:
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break ;
        case AUTOKEY_0:
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break;
        case AUTOKEY_1:
            if (cfg->user_plaintext_keyword_present) {
                make_keyed_alphabet(cfg->user_plaintext_keyword, current_plaintext_keyword_state);
            } else {
                random_keyword(current_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
            }
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break;
        case AUTOKEY_2:
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            if (cfg->user_ciphertext_keyword_present) {
                make_keyed_alphabet(cfg->user_ciphertext_keyword, current_ciphertext_keyword_state);
            } else {
                random_keyword(current_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
            }
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break;
        case AUTOKEY_3:
            if (cfg->user_plaintext_keyword_present) {
                make_keyed_alphabet(cfg->user_plaintext_keyword, current_plaintext_keyword_state);
            } else if (cfg->user_ciphertext_keyword_present) {
                make_keyed_alphabet(cfg->user_ciphertext_keyword, current_plaintext_keyword_state);
            } else {
                random_keyword(current_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
            }
            vec_copy(current_plaintext_keyword_state, current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break;
        case AUTOKEY_4:
            if (cfg->user_plaintext_keyword_present) {
                make_keyed_alphabet(cfg->user_plaintext_keyword, current_plaintext_keyword_state);
            } else {
                random_keyword(current_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
            }
            if (cfg->user_ciphertext_keyword_present) {
                make_keyed_alphabet(cfg->user_ciphertext_keyword, current_ciphertext_keyword_state);
            } else {
                random_keyword(current_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
            }
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break;
        case AUTOKEY_BEAU:
        case AUTOKEY_PORTA:
            plaintext_keyword_len = g_alpha;
            ciphertext_keyword_len = g_alpha;
            straight_alphabet(current_plaintext_keyword_state, g_alpha);
            straight_alphabet(current_ciphertext_keyword_state, g_alpha);
            random_cycleword(current_cycleword_state, g_alpha, cycleword_len);
            break;
    }

    if (cfg->same_key_cycle) {
        vec_copy(current_plaintext_keyword_state, current_ciphertext_keyword_state, g_alpha);
        vec_copy(current_ciphertext_keyword_state, current_cycleword_state, g_alpha);
    }

    if (cfg->optimal_cycleword && ! is_autokey) {
        derive_optimal_cycleword(cfg, ctx->cipher, ctx->cipher_len,
            current_plaintext_keyword_state, current_ciphertext_keyword_state,
            current_cycleword_state, cycleword_len, ctx->hist_by_col);
    }
}

static void polyalpha_perturb(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                              bool *force_primary) {
    ColossusConfig *cfg = ctx->cfg;
    int *cipher_indices = ctx->cipher;
    int cipher_len = ctx->cipher_len;
    int *crib_indices = ctx->crib_indices, *crib_positions = ctx->crib_positions;
    int n_cribs = ctx->n_cribs;
    // st is the engine's `local` (already a copy of `current`); operate on it directly.
    int *local_plaintext_keyword_state = st->pt_keyword;
    int *local_ciphertext_keyword_state = st->ct_keyword;
    int *local_cycleword_state = st->cycleword;
    int cycleword_len = cc->period;
    int plaintext_keyword_len = cc->j, ciphertext_keyword_len = cc->k;
    bool is_autokey = poly_is_autokey(cfg->cipher_type);
    bool perturbate_keyword_p = *force_primary;
    bool did_perturb_keyword = false;
    bool contradiction;

    if (perturbate_keyword_p ||
            cfg->cipher_type == VIGENERE || cfg->cipher_type == GRONSFELD ||
            is_autokey || frand() < cfg->keyword_permutation_probability) {
        switch (cfg->cipher_type) {
            case VIGENERE:
            case GRONSFELD:
            case PORTA:
            case BEAUFORT:
            case AUTOKEY_0:
            case AUTOKEY_BEAU:
            case AUTOKEY_PORTA:
                did_perturb_keyword = false;
                break ;
            case QUAGMIRE_1:
                if (!cfg->user_plaintext_keyword_present) {
                    perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                    did_perturb_keyword = true;
                }
                break ;
            case QUAGMIRE_2:
                if (!cfg->user_ciphertext_keyword_present) {
                    perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                    did_perturb_keyword = true;
                }
                break ;
            case QUAGMIRE_3:
                if (!cfg->user_plaintext_keyword_present && !cfg->user_ciphertext_keyword_present) {
                    perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                    vec_copy(local_plaintext_keyword_state, local_ciphertext_keyword_state, g_alpha);
                    did_perturb_keyword = true;
                }
                break ;
            case QUAGMIRE_4:
                if (cfg->user_plaintext_keyword_present && cfg->user_ciphertext_keyword_present) {
                    did_perturb_keyword = false;
                } else if (cfg->user_plaintext_keyword_present) {
                    perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                    did_perturb_keyword = true;
                } else if (cfg->user_ciphertext_keyword_present) {
                    perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                    did_perturb_keyword = true;
                } else {
                    if (frand() < 0.5) {
                        perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                    } else {
                        perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                    }
                    did_perturb_keyword = true;
                }
                break ;
            case AUTOKEY_1:
                 if (!cfg->user_plaintext_keyword_present) {
                     perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                     did_perturb_keyword = true;
                 }
                 break;
            case AUTOKEY_2:
                 if (!cfg->user_ciphertext_keyword_present) {
                     perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                     did_perturb_keyword = true;
                 }
                 break;
            case AUTOKEY_3:
                 if (!cfg->user_plaintext_keyword_present && !cfg->user_ciphertext_keyword_present) {
                     perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                     vec_copy(local_plaintext_keyword_state, local_ciphertext_keyword_state, g_alpha);
                     did_perturb_keyword = true;
                 }
                 break;
            case AUTOKEY_4:
                 if (cfg->user_plaintext_keyword_present && cfg->user_ciphertext_keyword_present) {
                     did_perturb_keyword = false;
                 } else if (cfg->user_plaintext_keyword_present) {
                     perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                     did_perturb_keyword = true;
                 } else if (cfg->user_ciphertext_keyword_present) {
                     perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                     did_perturb_keyword = true;
                 } else {
                     if (frand() < 0.5) perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                     else perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                     did_perturb_keyword = true;
                 }
                 break;
        }
    } else {
        did_perturb_keyword = false;
    }

    if (cfg->optimal_cycleword && ! is_autokey) {
        if (!did_perturb_keyword && cfg->cipher_type != BEAUFORT && cfg->cipher_type != VIGENERE && cfg->cipher_type != GRONSFELD && cfg->cipher_type != PORTA) {
            if (cfg->cipher_type == QUAGMIRE_3 && !(cfg->user_plaintext_keyword_present || cfg->user_ciphertext_keyword_present)) {
                 perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                 vec_copy(local_plaintext_keyword_state, local_ciphertext_keyword_state, g_alpha);
                 did_perturb_keyword = true;
            }
            else if (cfg->cipher_type == QUAGMIRE_1 && !cfg->user_plaintext_keyword_present) {
                perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                did_perturb_keyword = true;
            }
            else if (cfg->cipher_type == QUAGMIRE_2 && !cfg->user_ciphertext_keyword_present) {
                perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                did_perturb_keyword = true;
            }
            else if (cfg->cipher_type == QUAGMIRE_4) {
                if (!cfg->user_plaintext_keyword_present && !cfg->user_ciphertext_keyword_present) {
                    if (frand() < 0.5) {
                        perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                    } else {
                        perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                    }
                    did_perturb_keyword = true;
                } else if (!cfg->user_plaintext_keyword_present) {
                     perturbate_keyword(local_plaintext_keyword_state, g_alpha, plaintext_keyword_len);
                     did_perturb_keyword = true;
                } else if (!cfg->user_ciphertext_keyword_present) {
                     perturbate_keyword(local_ciphertext_keyword_state, g_alpha, ciphertext_keyword_len);
                     did_perturb_keyword = true;
                }
            }
        }

        derive_optimal_cycleword(cfg, cipher_indices, cipher_len,
            local_plaintext_keyword_state, local_ciphertext_keyword_state,
            local_cycleword_state, cycleword_len, ctx->hist_by_col);

    } else {
        // Gronsfeld shifts are key digits, so its cycleword domain is 0..9 (not 0..25).
        int cw_max = (cfg->cipher_type == GRONSFELD) ? GRONSFELD_DIGITS : g_alpha;
        if (cfg->cipher_type == VIGENERE || cfg->cipher_type == GRONSFELD ||
            cfg->cipher_type == PORTA || is_autokey) {
             perturbate_cycleword(local_cycleword_state, cw_max, cycleword_len);
        }
        else if (!did_perturb_keyword) {
            perturbate_cycleword(local_cycleword_state, cw_max, cycleword_len);
        }

        if (cfg->cipher_type != VIGENERE && cfg->cipher_type != GRONSFELD && cfg->cipher_type != BEAUFORT && cfg->cipher_type != PORTA && ! is_autokey) {
            perturbate_keyword_p = false;

            if (did_perturb_keyword) {
                contradiction = constrain_cycleword(cfg, cipher_indices, cipher_len, crib_indices,
                    crib_positions, n_cribs,
                    local_plaintext_keyword_state, local_ciphertext_keyword_state,
                    local_cycleword_state, cycleword_len, cfg->variant, cfg->verbose);

                if (contradiction) {
                    if (ctx->model_scratch) (*(long *) ctx->model_scratch) += 1;
                    perturbate_keyword_p = true;
                }
            }
        }
    }

    if (cfg->same_key_cycle) {
        vec_copy(local_plaintext_keyword_state, local_ciphertext_keyword_state, g_alpha);
        vec_copy(local_ciphertext_keyword_state, local_cycleword_state, g_alpha);
    }

    *force_primary = perturbate_keyword_p;
}

static void polyalpha_copy(const SolverConfig *cc, const SolverState *src, SolverState *dst) {
    for (int i = 0; i < g_alpha; i++) {
        dst->pt_keyword[i] = src->pt_keyword[i];
        dst->ct_keyword[i] = src->ct_keyword[i];
    }
    for (int i = 0; i < cc->period; i++) dst->cycleword[i] = src->cycleword[i];
}

static void polyalpha_decrypt(const SolverCtx *ctx, const SolverConfig *cc, SolverState *st,
                              int *out, double *score_adjust) {
    (void) score_adjust;
    decrypt_state(ctx->cfg, ctx->cipher, ctx->cipher_len,
        st->pt_keyword, st->ct_keyword, st->cycleword, cc->period, out);
}

static void polyalpha_report_verbose(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted, const EngineStats *stats) {
    (void) decrypted;
    ColossusConfig *cfg = ctx->cfg;
    int cipher_len = ctx->cipher_len;
    int *cipher_indices = ctx->cipher;
    int cycleword_len = cc->period;
    bool is_autokey = poly_is_autokey(cfg->cipher_type);
    int buf[MAX_CIPHER_LENGTH];
    int j, k, indx, offset = 0;

    const int *best_plaintext_keyword_state = st->pt_keyword;
    const int *best_ciphertext_keyword_state = st->ct_keyword;
    const int *best_cycleword_state = st->cycleword;

    if (cfg->cipher_type == PORTA) {
        porta_decrypt(buf, cipher_indices, cipher_len, (int *) best_cycleword_state, cycleword_len);
    } else if (cfg->cipher_type == BEAUFORT) {
        beaufort_decrypt(buf, cipher_indices, cipher_len, (int *) best_cycleword_state, cycleword_len);
    } else if (is_autokey) {
        autokey_decrypt(cfg, buf, cipher_indices, cipher_len,
            (int *) best_plaintext_keyword_state, (int *) best_ciphertext_keyword_state,
            (int *) best_cycleword_state, cycleword_len);
    } else {
        quagmire_decrypt(buf, cipher_indices, cipher_len,
            (int *) best_plaintext_keyword_state, (int *) best_ciphertext_keyword_state,
            (int *) best_cycleword_state, cycleword_len, cfg->variant);
    }

    if (cfg->transperoffset_present)
        transperoffset(buf, cipher_len, cfg->trans_period, cfg->trans_offset);
    if (cfg->transmatrix_present)
        transmatrix(buf, cipher_len, cfg->trans_w1, cfg->trans_w2, cfg->trans_clockwise);

    double ioc = index_of_coincidence(buf, cipher_len);
    double chi = chi_squared(buf, cipher_len);
    double entropy_score = entropy(buf, cipher_len);
    double elapsed = engine_elapsed_sec(stats);
    double n_iter_per_sec = ((double) stats->n_iterations)/elapsed;

    printf("\n%.2f\t[sec]\n", elapsed);
    printf("%.0fK\t[it/sec/thread]\n", 1.e-3*n_iter_per_sec);
    printf("%d\t[backtracks]\n", stats->n_backtracks);
    printf("%d\t[restarts]\n", stats->n_restarts);
    printf("%d\t[slips]\n", stats->n_slips);
    long n_contradictions = ctx->model_scratch ? *(long *) ctx->model_scratch : 0;
    printf("%.2f\t[contradiction pct]\n", ((double) n_contradictions)/stats->n_iterations);
    printf("%.4f\t[IOC]\n", ioc);
    printf("%.4f\t[entropy]\n", entropy_score);
    printf("%.2f\t[chi-squared]\n", chi);
    printf("%.2f\t[score]\n", score);

    if (cfg->cipher_type != PORTA) {
        print_text((int *) best_plaintext_keyword_state, g_alpha); printf("\n");
        print_text((int *) best_ciphertext_keyword_state, g_alpha); printf("\n");
    }
    print_text((int *) best_cycleword_state, cycleword_len); printf("\n");

    printf("\n");
    if (cfg->cipher_type != PORTA) {
        for (k = 0; k < cycleword_len; k++) {
            for (j = 0; j < g_alpha; j++) {
                if (best_ciphertext_keyword_state[j] == best_cycleword_state[k]) offset = j;
            }
            for (j = 0; j < g_alpha; j++) {
                indx = (j + offset) % g_alpha;
                printf("%c", index_to_char(best_ciphertext_keyword_state[indx]));
            }
            printf("\n");
        }
    }
    printf("\n");
    print_solution_check(buf, cipher_len);
    print_text(buf, cipher_len); printf("\n");
    fflush(stdout);
}

static void polyalpha_report(const SolverCtx *ctx, const SolverConfig *cc,
        const SolverState *st, double score, int *decrypted) {
    ColossusConfig *cfg = ctx->cfg;
    SharedData *shared = ctx->shared;
    int cipher_len = ctx->cipher_len;
    int n_words_found = 0;

    char plaintext_string[MAX_CIPHER_LENGTH];
    for (int i = 0; i < cipher_len; i++) plaintext_string[i] = index_to_char(decrypted[i]);
    plaintext_string[cipher_len] = '\0';

    if (cfg->dictionary_present && shared->dict != NULL) {
        n_words_found = find_dictionary_words(plaintext_string, shared->dict,
            shared->n_dict_words, shared->max_dict_word_len);
    }

    SolveResult local_res;
    SolveResult *res = ctx->result ? ctx->result : &local_res;
    res->solved = true;
    res->cipher_type = cfg->cipher_type;
    res->score = score;
    res->n_words = n_words_found;
    res->cycleword_len = cc->period;
    vec_copy((int *) st->pt_keyword, res->plaintext_keyword, g_alpha);
    vec_copy((int *) st->ct_keyword, res->ciphertext_keyword, g_alpha);
    vec_copy((int *) st->cycleword, res->cycleword, cc->period);
    vec_copy(decrypted, res->decrypted, cipher_len);
    res->decrypted_len = cipher_len;

    report_solution(cfg, ctx->cribtext, ctx->cipher, res);
}

static const CipherModel POLYALPHA_MODEL = {
    .name = "polyalphabetic",
    .shape = SHAPE_SHOTGUN,          // overridden per-solve (DETERMINISTIC for optimal Vig/Beau/Porta)
    .needs_hist = false,             // set per-solve (optimal && !autokey)
    .enumerate_configs = polyalpha_enumerate,
    .key_len = NULL,
    .seed = polyalpha_seed,
    .perturb = polyalpha_perturb,
    .copy_state = polyalpha_copy,
    .decrypt = polyalpha_decrypt,
    .report = polyalpha_report,
    .report_verbose = polyalpha_report_verbose,
};


// Core Solver

// Prints the human-readable result block and the ">>> ..." one-line CSV summary
// from a populated SolveResult. Output is byte-for-byte the same as the inline
// reporting it replaced, so batch grep/sort over the ">>>" lines is unaffected.
void report_solution(ColossusConfig *cfg, char *cribtext_str,
    int cipher_indices[], SolveResult *res) {

    int cipher_len = res->decrypted_len;
    int n_cribs = 0;
    for (int i = 0; cribtext_str[i] != '\0'; i++) if (cribtext_str[i] != '_') n_cribs++;

    if (cfg->transperoffset_present) {
        printf("\ntransperiodoffset: period = %d, offset = %d\n", cfg->trans_period, cfg->trans_offset);
    }

    if (cfg->transmatrix_present) {
        printf("\ntransmatrix: w1 = %d, w2 = %d, direction = %s\n", cfg->trans_w1, cfg->trans_w2, cfg->trans_clockwise ? "cw" : "ccw");
    }

    // Results Output
    printf("\nResult Score: %.2f | Words: %d\n", res->score, res->n_words);

    print_text(cipher_indices, cipher_len);
    printf("\n");

    if (cfg->cipher_type != PORTA) {
        print_text(res->plaintext_keyword, g_alpha);
        printf("\n");
        print_text(res->ciphertext_keyword, g_alpha);
        printf("\n");
    }

    print_text(res->cycleword, res->cycleword_len);
    printf("\n");
    if (cfg->cipher_type == GRONSFELD) {
        // The Gronsfeld key is numeric: show the recovered per-column shifts as digits.
        printf("key (numeric): ");
        for (int i = 0; i < res->cycleword_len; i++) printf("%d", res->cycleword[i]);
        printf("\n");
    }
    print_text(res->decrypted, cipher_len);
    printf("\n");
    print_solution_check(res->decrypted, cipher_len);
    print_spaces_line(g_spaces_table, res->decrypted, cipher_len);
    printf("%s\n", cribtext_str);

    if (PARTIAL_CRIB_MATCH && n_cribs > 0) {
        // Index the crib by cipher position via cribtext_str ('_' = no crib),
        // not the positionally-packed crib_indices array. 0 = exact match,
        // a small digit = near miss, '*' = far (>= 10) from the crib char.
        for (int i = 0; i < cipher_len; i++) {
            if (cribtext_str[i] == '_') {
                printf("_"); // No crib defined for this position.
            } else {
                int diff = abs(res->decrypted[i] - (g_char_to_idx[toupper((unsigned char)cribtext_str[i]) & 127]));
                if (diff < 10) {
                    printf("%d", diff);
                } else {
                    printf("*");
                }
            }
        }
    }
    printf("\n\n");

    // One-liner summary
    if (cfg->transperoffset_present) {
        if (cfg->dictionary_present) {
            printf(">>> %.2f, %d, %d, %d, %d, %s, ", res->score, res->n_words, cfg->cipher_type, cfg->trans_period, cfg->trans_offset, cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        } else {
            printf(">>> %.2f, %d, %d, %d, %s, ", res->score, cfg->cipher_type, cfg->trans_period, cfg->trans_offset, cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        }
    } else if (cfg->transmatrix_present) {
        if (cfg->dictionary_present) {
            printf(">>> %.2f, %d, %d, %d, %d, %d, %s, ", res->score, res->n_words, cfg->cipher_type, cfg->trans_w1, cfg->trans_w2, cfg->trans_clockwise, cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        } else {
            printf(">>> %.2f, %d, %d, %d, %d, %s, ", res->score, cfg->cipher_type, cfg->trans_w1, cfg->trans_w2, cfg->trans_clockwise, cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        }
    } else {
        if (cfg->dictionary_present) {
            printf(">>> %.2f, %d, %d, %s, ", res->score, res->n_words, cfg->cipher_type, cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        } else {
            printf(">>> %.2f, %d, %s, ", res->score, cfg->cipher_type, cfg->batch_present ? "BATCH" : cfg->ciphertext_file);
        }
    }

    print_text(cipher_indices, cipher_len);
    printf(", ");

    if (cfg->cipher_type != PORTA) {
        print_text(res->plaintext_keyword, g_alpha);
        printf(", ");
        print_text(res->ciphertext_keyword, g_alpha);
        printf(", ");
    }

    print_text(res->cycleword, res->cycleword_len);
    printf(", ");

    print_text(res->decrypted, cipher_len);
    printf("\n");
}







// derive_optimal_cycleword() lives in optimal_cycleword.c (prototype in the
// shared header). It deterministically solves each period column for the key
// character that best matches English monogram frequencies.


int get_matrix_rotate_old_idx(int target_idx, int len, int width, int clockwise) {
    if (width <= 1 || width >= len) return target_idx;
    
    int R = (len + width - 1) / width;
    int W = width;
    int current_idx = 0;

    if (clockwise) {
        // Trace left-to-right, bottom-to-top
        for (int c = 0; c < W; c++) {
            for (int r = R - 1; r >= 0; r--) {
                int old_idx = r * W + c;
                if (old_idx < len) {
                    if (current_idx == target_idx) return old_idx;
                    current_idx++;
                }
            }
        }
    } else {
        // Trace right-to-left, top-to-bottom
        for (int c = W - 1; c >= 0; c--) {
            for (int r = 0; r < R; r++) {
                int old_idx = r * W + c;
                if (old_idx < len) {
                    if (current_idx == target_idx) return old_idx;
                    current_idx++;
                }
            }
        }
    }
    return target_idx;
}

int map_crib_to_cipher_pos(ColossusConfig *cfg, int crib_pos, int cipher_len) {
    int pos = crib_pos;

    // Un-map transmatrix (in reverse order of decryption: w2, then w1)
    if (cfg->transmatrix_present) {
        pos = get_matrix_rotate_old_idx(pos, cipher_len, cfg->trans_w2, cfg->trans_clockwise);
        pos = get_matrix_rotate_old_idx(pos, cipher_len, cfg->trans_w1, cfg->trans_clockwise);
    }

    // Un-map transperoffset.
    if (cfg->transperoffset_present) {
        pos = (pos + cfg->trans_offset) % cipher_len;
        if (pos < 0) pos += cipher_len;
        pos = (cfg->trans_period * pos) % cipher_len;
    }

    return pos;
}



bool cribs_satisfied_p(ColossusConfig *cfg, int cipher_indices[], int cipher_len, int crib_indices[], 
    int crib_positions[], int n_cribs, int cycleword_len, bool verbose) {

    int i, j, k, ii, jj, total, column_length, ciphertext_column_indices[MAX_CIPHER_LENGTH], 
        ciphertext_column[MAX_CIPHER_LENGTH], crib_frequencies[ALPHABET_SIZE][ALPHABET_SIZE];
    int mapped_crib_positions[MAX_CIPHER_LENGTH];

    if (n_cribs == 0) return true;

    // Map crib positions to their original ciphertext positions
    for (i = 0; i < n_cribs; i++) {
        mapped_crib_positions[i] = map_crib_to_cipher_pos(cfg, crib_positions[i], cipher_len);
    }

    for (j = 0; j < cycleword_len; j++) {
        if (verbose) {
            printf("\nCOLUMN = %d \n", j);
        }

        k = 0;
        while (cycleword_len*k + j < cipher_len) {
            ciphertext_column_indices[k] = cycleword_len*k + j;
            ciphertext_column[k] = cipher_indices[ciphertext_column_indices[k]];
            k++;
        }
        column_length = k;

        for (i = 0; i < g_alpha; i++) {
            for (k = 0; k < g_alpha; k++) {
                crib_frequencies[i][k] = 0;
            }
        }

        for (i = 0; i < n_cribs; i++) {
            for (k = 0; k < column_length; k++) {
                // Use mapped_crib_positions here instead of crib_positions
                if (mapped_crib_positions[i] == ciphertext_column_indices[k]) {
                    if (verbose) {
                        printf("CT = %c, PT = %c\n", index_to_char(ciphertext_column[k]), index_to_char(crib_indices[i]));
                    }
                    
                    crib_frequencies[crib_indices[i]][ciphertext_column[k]] = 1;

                    for (ii = 0; ii < g_alpha; ii++) {
                        total = 0;
                        for (jj = 0; jj < g_alpha; jj++) {
                            total += crib_frequencies[ii][jj];
                            if (total > 1) {
                                if (verbose) {
                                    printf("\n\nContradiction at col %d, crib char %c\n\n", j, index_to_char(crib_indices[i]));
                                }
                                return false;
                            }
                        }
                    }

                    for (jj = 0; jj < g_alpha; jj++) {
                        total = 0;
                        for (ii = 0; ii < g_alpha; ii++) {
                            total += crib_frequencies[ii][jj];
                            if (total > 1) {
                                if (verbose) {
                                    printf("\n\nContradiction at col %d, crib char %c\n\n", j, index_to_char(crib_indices[i]));
                                }
                                return false;
                            }
                        }
                    }                   
                }
            }
        }
    }
    return true;
}



bool constrain_cycleword(ColossusConfig *cfg, int cipher_indices[], int cipher_len, 
    int crib_indices[], int crib_positions[], int n_cribs, 
    int plaintext_keyword_indices[], int ciphertext_keyword_indices[], 
    int cycleword_indices[], int cycleword_len, 
    bool variant, bool verbose) {

    int i, j, k, crib_char, ciphertext_char, posn_keyword, posn_cycleword, 
        indx, crib_cyclewords[MAX_CYCLEWORD_LEN], mapped_pos;

    if (n_cribs == 0) return false; 

    for (i = 0; i < cycleword_len; i++) crib_cyclewords[i] = INACTIVE; 

    for (i = 0; i < cycleword_len; i++) {
        for (j = 0; j < n_cribs; j++) {
            
            mapped_pos = map_crib_to_cipher_pos(cfg, crib_positions[j], cipher_len);
            
            // Check against the mapped ciphertext position
            if (mapped_pos % cycleword_len == i) {
                crib_char = crib_indices[j];
                ciphertext_char = cipher_indices[mapped_pos];
                
                for (k = 0; k < g_alpha; k++) {
                    if (ciphertext_keyword_indices[k] == ciphertext_char) {
                        posn_keyword = k; 
                        break ;
                    }
                }
                for (k = 0; k < g_alpha; k++) {
                    if (plaintext_keyword_indices[k] == crib_char) {
                        posn_cycleword = k; 
                        break ;
                    }
                }

                if (variant) {
                    indx = (posn_cycleword - posn_keyword) % g_alpha;
                } else {
                    indx = (posn_keyword - posn_cycleword) % g_alpha;
                }
                
                if (indx < 0) indx += g_alpha;                    

                if (crib_cyclewords[i] == INACTIVE) {
                    crib_cyclewords[i] = ciphertext_keyword_indices[indx];
                    cycleword_indices[i] = ciphertext_keyword_indices[indx]; 
                } else if (crib_cyclewords[i] != ciphertext_keyword_indices[indx]) { 
                    /*
                    if (verbose) {
                        printf("\n\nContradiction at crib %c, posn %d; rejecting keyword ", 
                            index_to_char(crib_indices[j]), crib_positions[j]);
                    } */
                    return true; 
                }
            }
        }
    }
    return false;
}



void decrypt_state(ColossusConfig *cfg, int cipher_indices[], int cipher_len, 
                   int plaintext_keyword_state[], int ciphertext_keyword_state[], 
                   int cycleword_state[], int cycleword_len, 
                   int decrypted[]) {
                   
    bool is_autokey = ((cfg->cipher_type >= AUTOKEY_0 && cfg->cipher_type <= AUTOKEY_4) || 
        cfg->cipher_type == AUTOKEY_BEAU || cfg->cipher_type == AUTOKEY_PORTA);

    if (cfg->cipher_type == PORTA) { 
        porta_decrypt(decrypted, cipher_indices, cipher_len, 
                     cycleword_state, cycleword_len);
    } else if (cfg->cipher_type == BEAUFORT) { 
        beaufort_decrypt(decrypted, cipher_indices, cipher_len, 
                     cycleword_state, cycleword_len);
    } else if (is_autokey) {
        autokey_decrypt(cfg, decrypted, cipher_indices, cipher_len, 
            plaintext_keyword_state, ciphertext_keyword_state,
            cycleword_state, cycleword_len);
    } else if (cfg->cipher_type == VIGENERE) {
        vigenere_decrypt(decrypted, cipher_indices, cipher_len,
                         cycleword_state, cycleword_len, cfg->variant);
    } else if (cfg->cipher_type == GRONSFELD) {
        // Gronsfeld = Vigenere with numeric key digits carried in the cycleword lane.
        gronsfeld_decrypt(decrypted, cipher_indices, cipher_len,
                          cycleword_state, cycleword_len);
    } else {
        // Quagmire I-IV
        quagmire_decrypt(decrypted, cipher_indices, cipher_len, 
            plaintext_keyword_state, ciphertext_keyword_state, 
            cycleword_state, cycleword_len, cfg->variant);
    }

    // Apply transposition. 
    if (cfg->transperoffset_present) {
        transperoffset(decrypted, cipher_len, cfg->trans_period, cfg->trans_offset);
    }

    if (cfg->transmatrix_present) {
        transmatrix(decrypted, cipher_len, cfg->trans_w1, cfg->trans_w2, cfg->trans_clockwise);
    }
}





// Polyalphabetic solve entry, lifted out of solve_cipher so POLYALPHA_MODEL and the
// per-type seed/perturb ladders stay private to this module. Period estimation runs
// inside the engine via the model's enumerate hook.
void solve_polyalpha(char *ciphertext_str, char *cribtext_str, ColossusConfig *cfg,
    SharedData *shared, int cipher_indices[], int cipher_len,
    int crib_indices[], int crib_positions[], int n_cribs, SolveResult *result) {

    (void) ciphertext_str;
    // --- POLYALPHABETIC CIPHERS ---
    // Vigenere / Quagmire I-IV / Beaufort / Porta / Autokey* all share one model;
    // the cipher-agnostic engine enumerates (cycleword_len, j, k) and shotgun /
    // hill-climbs each via the model's seed/perturb/decrypt hooks.
    SolverCtx ctx = make_solver_ctx(cfg, shared, cribtext_str,
        cipher_indices, cipher_len, crib_indices, crib_positions, n_cribs);
    ctx.result = result;

    bool is_autokey = poly_is_autokey(cfg->cipher_type);
    long poly_contradictions = 0;            // verbose-only telemetry for the perturb hook
    ctx.model_scratch = &poly_contradictions;

    CipherModel model = POLYALPHA_MODEL;
    model.needs_hist = cfg->optimal_cycleword && !is_autokey;
    // Vigenere/Beaufort/Porta in optimal mode use fixed straight alphabets and a
    // deterministically-derived cycleword, so the climb is deterministic: stop at
    // the first recorded best. (-samekey breaks the determinism, so exclude it.)
    if (cfg->optimal_cycleword && !is_autokey && !cfg->same_key_cycle &&
        (cfg->cipher_type == VIGENERE || cfg->cipher_type == BEAUFORT ||
         cfg->cipher_type == PORTA)) {
        model.shape = SHAPE_DETERMINISTIC;
    }

    run_solver(&model, &ctx);
}

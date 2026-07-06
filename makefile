#
#   makefile
#

CC=gcc -Wall -O3 -funroll-loops

# CC=gcc -Wall -lm -g -O0

# Sources live under src/<cipher-class>/. All local #includes are flat
# (#include "foo.h"), so the compiler finds every header via these -I paths
# regardless of which subdirectory it lives in.
SRC=src
INCLUDES=-I$(SRC)/core -I$(SRC)/polyalphabetic -I$(SRC)/transposition -I$(SRC)/polygraphic -I$(SRC)/substitution

CORE=$(SRC)/core
POLY=$(SRC)/polyalphabetic
TRANS=$(SRC)/transposition
GRAPH=$(SRC)/polygraphic
SUBST=$(SRC)/substitution

# Cipher primitives (decrypt math, shared by the solvers and the unit tests).
PRIMITIVES=$(CORE)/utils.c $(CORE)/parse.c $(CORE)/dict.c $(TRANS)/transpositions.c $(TRANS)/period_column.c $(CORE)/perioc.c $(POLY)/quagmire.c $(POLY)/vigenere.c $(POLY)/gronsfeld.c $(POLY)/gromark.c $(POLY)/nicodemus.c $(POLY)/porta.c $(POLY)/beaufort.c $(POLY)/autokey.c $(CORE)/optimal_cycleword.c $(GRAPH)/playfair.c $(GRAPH)/bifid.c $(GRAPH)/trifid.c $(GRAPH)/hill.c $(GRAPH)/phillips.c $(GRAPH)/twosquare.c $(GRAPH)/foursquare.c $(GRAPH)/adfgvx.c $(GRAPH)/nihilist_sub.c $(GRAPH)/bazeries.c $(GRAPH)/portax.c $(POLY)/progkey.c $(GRAPH)/slidefair.c $(GRAPH)/seriated_playfair.c $(GRAPH)/digrafid.c $(GRAPH)/cm_bifid.c $(GRAPH)/trisquare.c $(POLY)/intkey.c $(POLY)/condi.c $(GRAPH)/fracmorse.c

# Cipher-agnostic core + per-cipher-type solver modules (split out of colossus.c).
SOLVERS=$(CORE)/engine.c $(CORE)/scoring.c $(TRANS)/trans_common.c $(POLY)/polyalpha_solver.c $(POLY)/gromark_solver.c $(POLY)/nicodemus_solver.c $(TRANS)/transmatrix_solver.c $(TRANS)/permutation_solver.c $(TRANS)/columnar_solver.c $(TRANS)/columnar_track_solver.c $(TRANS)/route_chain_solver.c $(TRANS)/tile_solver.c $(TRANS)/period_column_solver.c $(TRANS)/period_column_space_solver.c $(TRANS)/double_transposition_solver.c $(TRANS)/railfence_solver.c $(TRANS)/route_solver.c $(TRANS)/amsco_solver.c $(TRANS)/myszkowski_solver.c $(TRANS)/redefence_solver.c $(TRANS)/cadenus_solver.c $(TRANS)/nihilist_solver.c $(TRANS)/swagman_solver.c $(TRANS)/grille_solver.c $(SUBST)/indep_solver.c $(SUBST)/homophonic_solver.c $(GRAPH)/playfair_solver.c $(GRAPH)/bifid_solver.c $(GRAPH)/trifid_solver.c $(GRAPH)/hill_solver.c $(GRAPH)/phillips_solver.c $(GRAPH)/twosquare_solver.c $(GRAPH)/foursquare_solver.c $(GRAPH)/adfgvx_solver.c $(GRAPH)/nihilist_sub_solver.c $(GRAPH)/bazeries_solver.c $(GRAPH)/portax_solver.c $(POLY)/progkey_solver.c $(GRAPH)/slidefair_solver.c $(GRAPH)/seriated_playfair_solver.c $(GRAPH)/digrafid_solver.c $(GRAPH)/cm_bifid_solver.c $(GRAPH)/trisquare_solver.c $(POLY)/intkey_solver.c $(POLY)/condi_solver.c $(GRAPH)/fracmorse_solver.c

# The full solver translation-unit set (everything but the test harnesses).
SOLVER_SRC=$(PRIMITIVES) $(SOLVERS) $(CORE)/colossus.c

all:
	$(CC) $(INCLUDES) $(SOLVER_SRC) -o colossus
	cp colossus ..
	cp colossus ../quagmire

# Fast unit tests of the primitives (sub-second). Add -lm on Linux.
#   test_transpositions   : the transposition primitives (transpositions.c)
#   test_ciphers          : the cipher primitives (vigenere/beaufort/porta/quagmire/autokey)
#   test_optimal_cycleword : deterministic optimal-cycleword recovery
#   test_playfair         : the Playfair primitives (grid build / encrypt / decrypt / prepare)
test:
	$(CC) $(INCLUDES) tests/test_transpositions.c $(CORE)/utils.c $(TRANS)/transpositions.c -o tests/test_transpositions
	./tests/test_transpositions
	$(CC) $(INCLUDES) tests/test_ciphers.c $(CORE)/utils.c $(POLY)/quagmire.c $(POLY)/vigenere.c $(POLY)/porta.c $(POLY)/beaufort.c $(POLY)/autokey.c -o tests/test_ciphers
	./tests/test_ciphers
	$(CC) $(INCLUDES) tests/test_gronsfeld.c $(CORE)/utils.c $(POLY)/gronsfeld.c $(POLY)/vigenere.c $(POLY)/quagmire.c -o tests/test_gronsfeld
	./tests/test_gronsfeld
	$(CC) $(INCLUDES) tests/test_optimal_cycleword.c $(CORE)/utils.c $(POLY)/quagmire.c $(POLY)/vigenere.c $(POLY)/porta.c $(POLY)/beaufort.c $(CORE)/optimal_cycleword.c -o tests/test_optimal_cycleword
	./tests/test_optimal_cycleword
	$(CC) $(INCLUDES) tests/test_held_karp.c $(TRANS)/trans_common.c $(CORE)/scoring.c $(CORE)/utils.c $(CORE)/dict.c -o tests/test_held_karp
	./tests/test_held_karp
	$(CC) $(INCLUDES) tests/test_playfair.c $(CORE)/utils.c $(GRAPH)/playfair.c -o tests/test_playfair
	./tests/test_playfair
	$(CC) $(INCLUDES) tests/test_bifid.c $(CORE)/utils.c $(GRAPH)/bifid.c -o tests/test_bifid
	./tests/test_bifid
	$(CC) $(INCLUDES) tests/test_trifid.c $(CORE)/utils.c $(GRAPH)/trifid.c -o tests/test_trifid
	./tests/test_trifid
	$(CC) $(INCLUDES) tests/test_hill.c $(CORE)/utils.c $(GRAPH)/hill.c -o tests/test_hill
	./tests/test_hill
	$(CC) $(INCLUDES) tests/test_phillips.c $(CORE)/utils.c $(GRAPH)/phillips.c -o tests/test_phillips
	./tests/test_phillips
	$(CC) $(INCLUDES) tests/test_twosquare.c $(CORE)/utils.c $(GRAPH)/twosquare.c $(GRAPH)/playfair.c -o tests/test_twosquare
	./tests/test_twosquare
	$(CC) $(INCLUDES) tests/test_foursquare.c $(CORE)/utils.c $(GRAPH)/foursquare.c -o tests/test_foursquare
	./tests/test_foursquare
	$(CC) $(INCLUDES) tests/test_trisquare.c $(CORE)/utils.c $(GRAPH)/trisquare.c $(GRAPH)/bifid.c -o tests/test_trisquare
	./tests/test_trisquare
	$(CC) $(INCLUDES) tests/test_adfgvx.c $(CORE)/utils.c $(GRAPH)/adfgvx.c $(GRAPH)/bifid.c $(TRANS)/transpositions.c -o tests/test_adfgvx
	./tests/test_adfgvx
	$(CC) $(INCLUDES) tests/test_nihilist_sub.c $(CORE)/utils.c $(GRAPH)/nihilist_sub.c -o tests/test_nihilist_sub
	./tests/test_nihilist_sub
	$(CC) $(INCLUDES) tests/test_gromark.c $(CORE)/utils.c $(POLY)/gromark.c -o tests/test_gromark
	./tests/test_gromark
	$(CC) $(INCLUDES) tests/test_nicodemus.c $(CORE)/utils.c $(POLY)/nicodemus.c $(POLY)/quagmire.c $(POLY)/vigenere.c $(POLY)/beaufort.c -o tests/test_nicodemus
	./tests/test_nicodemus
	$(CC) $(INCLUDES) tests/test_bazeries.c $(CORE)/utils.c $(GRAPH)/bazeries.c $(GRAPH)/bifid.c -o tests/test_bazeries
	./tests/test_bazeries
	$(CC) $(INCLUDES) tests/test_portax.c $(CORE)/utils.c $(GRAPH)/portax.c -o tests/test_portax
	./tests/test_portax
	$(CC) $(INCLUDES) tests/test_progkey.c $(CORE)/utils.c $(POLY)/progkey.c $(POLY)/vigenere.c $(POLY)/beaufort.c $(POLY)/quagmire.c -o tests/test_progkey
	./tests/test_progkey
	$(CC) $(INCLUDES) tests/test_intkey.c $(CORE)/utils.c $(POLY)/intkey.c $(POLY)/vigenere.c $(POLY)/beaufort.c $(POLY)/quagmire.c -o tests/test_intkey
	./tests/test_intkey
	$(CC) $(INCLUDES) tests/test_condi.c $(CORE)/utils.c $(POLY)/condi.c -o tests/test_condi
	./tests/test_condi
	$(CC) $(INCLUDES) tests/test_slidefair.c $(CORE)/utils.c $(GRAPH)/slidefair.c -o tests/test_slidefair
	./tests/test_slidefair
	$(CC) $(INCLUDES) tests/test_seriated_playfair.c $(CORE)/utils.c $(GRAPH)/seriated_playfair.c $(GRAPH)/playfair.c -o tests/test_seriated_playfair
	./tests/test_seriated_playfair
	$(CC) $(INCLUDES) tests/test_digrafid.c $(CORE)/utils.c $(GRAPH)/digrafid.c $(GRAPH)/bifid.c -o tests/test_digrafid
	./tests/test_digrafid
	$(CC) $(INCLUDES) tests/test_cm_bifid.c $(CORE)/utils.c $(GRAPH)/cm_bifid.c $(GRAPH)/bifid.c -o tests/test_cm_bifid
	./tests/test_cm_bifid
	$(CC) $(INCLUDES) tests/test_fracmorse.c $(CORE)/utils.c $(GRAPH)/fracmorse.c -o tests/test_fracmorse
	./tests/test_fracmorse
	$(CC) $(INCLUDES) tests/test_period_column.c $(CORE)/utils.c $(TRANS)/period_column.c -o tests/test_period_column
	./tests/test_period_column
	$(CC) $(INCLUDES) tests/test_double_transposition.c $(TRANS)/double_transposition_solver.c $(TRANS)/transpositions.c $(CORE)/scoring.c $(CORE)/utils.c $(CORE)/dict.c -o tests/test_double_transposition
	./tests/test_double_transposition

# Slow optimizer regression suite (~30s): planted-cipher recovery through the
# full solve_cipher hill climber at fixed seeds and budgets. Kept separate from
# `make test` so the fast primitive checks stay in the quick edit/build loop.
testopt:
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_solver.c $(SOLVER_SRC) -o tests/test_solver
	./tests/test_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_gronsfeld_solver.c $(SOLVER_SRC) -o tests/test_gronsfeld_solver
	./tests/test_gronsfeld_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_playfair_solver.c $(SOLVER_SRC) -o tests/test_playfair_solver
	./tests/test_playfair_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_bifid_solver.c $(SOLVER_SRC) -o tests/test_bifid_solver
	./tests/test_bifid_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_trifid_solver.c $(SOLVER_SRC) -o tests/test_trifid_solver
	./tests/test_trifid_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_hill_solver.c $(SOLVER_SRC) -o tests/test_hill_solver
	./tests/test_hill_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_phillips_solver.c $(SOLVER_SRC) -o tests/test_phillips_solver
	./tests/test_phillips_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_twosquare_solver.c $(SOLVER_SRC) -o tests/test_twosquare_solver
	./tests/test_twosquare_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_foursquare_solver.c $(SOLVER_SRC) -o tests/test_foursquare_solver
	./tests/test_foursquare_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_trisquare_solver.c $(SOLVER_SRC) -o tests/test_trisquare_solver
	./tests/test_trisquare_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_adfgvx_solver.c $(SOLVER_SRC) -o tests/test_adfgvx_solver
	./tests/test_adfgvx_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_nihilist_sub_solver.c $(SOLVER_SRC) -o tests/test_nihilist_sub_solver
	./tests/test_nihilist_sub_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_gromark_solver.c $(SOLVER_SRC) -o tests/test_gromark_solver
	./tests/test_gromark_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_nicodemus_solver.c $(SOLVER_SRC) -o tests/test_nicodemus_solver
	./tests/test_nicodemus_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_bazeries_solver.c $(SOLVER_SRC) -o tests/test_bazeries_solver
	./tests/test_bazeries_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_portax_solver.c $(SOLVER_SRC) -o tests/test_portax_solver
	./tests/test_portax_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_progkey_solver.c $(SOLVER_SRC) -o tests/test_progkey_solver
	./tests/test_progkey_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_intkey_solver.c $(SOLVER_SRC) -o tests/test_intkey_solver
	./tests/test_intkey_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_condi_solver.c $(SOLVER_SRC) -o tests/test_condi_solver
	./tests/test_condi_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_slidefair_solver.c $(SOLVER_SRC) -o tests/test_slidefair_solver
	./tests/test_slidefair_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_seriated_playfair_solver.c $(SOLVER_SRC) -o tests/test_seriated_playfair_solver
	./tests/test_seriated_playfair_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_digrafid_solver.c $(SOLVER_SRC) -o tests/test_digrafid_solver
	./tests/test_digrafid_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_cm_bifid_solver.c $(SOLVER_SRC) -o tests/test_cm_bifid_solver
	./tests/test_cm_bifid_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_fracmorse_solver.c $(SOLVER_SRC) -o tests/test_fracmorse_solver
	./tests/test_fracmorse_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_period_column_solver.c $(SOLVER_SRC) -o tests/test_period_column_solver
	./tests/test_period_column_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_period_column_space_solver.c $(SOLVER_SRC) -o tests/test_period_column_space_solver
	./tests/test_period_column_space_solver
	$(CC) $(INCLUDES) -DCOLOSSUS_NO_MAIN tests/test_double_transposition_solver.c $(SOLVER_SRC) -o tests/test_double_transposition_solver
	./tests/test_double_transposition_solver

# Everything.
testall: test testopt

# Standalone test-data generator for homophonic ciphers (not part of the solver
# build). Emits a comma-separated homophonic ciphertext + its plaintext solution;
# used to mint ciphers/tests/homophonic_test.*.
homophonic_gen:
	$(CC) $(INCLUDES) tools/homophonic_gen.c -o tools/homophonic_gen

# Standalone test-data generator for Playfair ciphers. Reuses the real cipher code
# (playfair.c + utils.c) so the generator and solver can never drift in convention.
playfair_gen:
	$(CC) $(INCLUDES) tools/playfair_gen.c $(GRAPH)/playfair.c $(CORE)/utils.c -o tools/playfair_gen

# Standalone test-data generator for Bifid ciphers. Reuses the real cipher code
# (bifid.c + utils.c) so the generator and solver can never drift in convention.
bifid_gen:
	$(CC) $(INCLUDES) tools/bifid_gen.c $(GRAPH)/bifid.c $(CORE)/utils.c -o tools/bifid_gen

# Standalone test-data generator for Trifid ciphers. Reuses the real cipher code
# (trifid.c + utils.c) so the generator and solver can never drift in convention.
trifid_gen:
	$(CC) $(INCLUDES) tools/trifid_gen.c $(GRAPH)/trifid.c $(CORE)/utils.c -o tools/trifid_gen

# Standalone test-data generator for Hill ciphers. Reuses the real cipher code
# (hill.c + utils.c) so the generator and solver can never drift in convention.
hill_gen:
	$(CC) $(INCLUDES) tools/hill_gen.c $(GRAPH)/hill.c $(CORE)/utils.c -o tools/hill_gen

# Standalone test-data generator for Gronsfeld ciphers. Reuses the real cipher code
# (gronsfeld.c + utils.c) so the generator and solver can never drift in convention.
gronsfeld_gen:
	$(CC) $(INCLUDES) tools/gronsfeld_gen.c $(POLY)/gronsfeld.c $(CORE)/utils.c -o tools/gronsfeld_gen

# Standalone test-data generator for Phillips ciphers. Reuses the real cipher code
# (phillips.c + utils.c) so the generator and solver can never drift in convention.
phillips_gen:
	$(CC) $(INCLUDES) tools/phillips_gen.c $(GRAPH)/phillips.c $(CORE)/utils.c -o tools/phillips_gen

# Standalone test-data generators for Two-Square / Four-Square ciphers. Reuse the real
# cipher code (twosquare.c / foursquare.c + playfair.c keyword build + utils.c) so the
# generators and the solver can never drift in convention.
twosquare_gen:
	$(CC) $(INCLUDES) tools/twosquare_gen.c $(GRAPH)/twosquare.c $(GRAPH)/playfair.c $(CORE)/utils.c -o tools/twosquare_gen

foursquare_gen:
	$(CC) $(INCLUDES) tools/foursquare_gen.c $(GRAPH)/foursquare.c $(GRAPH)/playfair.c $(CORE)/utils.c -o tools/foursquare_gen

trisquare_gen:
	$(CC) $(INCLUDES) tools/trisquare_gen.c $(GRAPH)/trisquare.c $(GRAPH)/bifid.c $(GRAPH)/playfair.c $(CORE)/utils.c -o tools/trisquare_gen

# Standalone test-data generator for ADFGVX / ADFGX ciphers. Reuses the real cipher code
# (adfgvx.c + bifid.c keyed-square build/inverse + transpositions.c columnar + utils.c)
# so the generator and the solver can never drift in convention.
adfgvx_gen:
	$(CC) $(INCLUDES) tools/adfgvx_gen.c $(GRAPH)/adfgvx.c $(GRAPH)/bifid.c $(TRANS)/transpositions.c $(CORE)/utils.c -o tools/adfgvx_gen

# Standalone test-data generator for Nihilist Substitution ciphers. Reuses the real cipher
# code (nihilist_sub.c + bifid.c keyed-square build/inverse + utils.c) so the generator and
# the solver can never drift in convention.
nihilist_sub_gen:
	$(CC) $(INCLUDES) tools/nihilist_sub_gen.c $(GRAPH)/nihilist_sub.c $(GRAPH)/bifid.c $(CORE)/utils.c -o tools/nihilist_sub_gen

# Standalone test-data generator for Gromark / Periodic Gromark ciphers. Reuses the real
# cipher code (gromark.c + utils.c) so the generator and the solver can never drift.
gromark_gen:
	$(CC) $(INCLUDES) tools/gromark_gen.c $(POLY)/gromark.c $(CORE)/utils.c -o tools/gromark_gen

# Standalone test-data generator for Nicodemus ciphers. Reuses the real cipher code
# (nicodemus.c + utils.c) so the generator and the solver can never drift.
nicodemus_gen:
	$(CC) $(INCLUDES) tools/nicodemus_gen.c $(POLY)/nicodemus.c $(CORE)/utils.c -o tools/nicodemus_gen

# Standalone test-data generator for Bazeries ciphers. Reuses the real cipher code
# (bazeries.c + bifid.c keyed-square build/inverse + utils.c) so the generator and the
# solver can never drift in convention.
bazeries_gen:
	$(CC) $(INCLUDES) tools/bazeries_gen.c $(GRAPH)/bazeries.c $(GRAPH)/bifid.c $(CORE)/utils.c -o tools/bazeries_gen

# Standalone test-data generator for Portax ciphers. Reuses the real cipher code
# (portax.c + utils.c) so the generator and the solver can never drift in convention.
portax_gen:
	$(CC) $(INCLUDES) tools/portax_gen.c $(GRAPH)/portax.c $(CORE)/utils.c -o tools/portax_gen

# Standalone test-data generator for Progressive Key ciphers. Reuses the real cipher code
# (progkey.c + utils.c) so the generator and the solver can never drift in convention.
progkey_gen:
	$(CC) $(INCLUDES) tools/progkey_gen.c $(POLY)/progkey.c $(CORE)/utils.c -o tools/progkey_gen

intkey_gen:
	$(CC) $(INCLUDES) tools/intkey_gen.c $(POLY)/intkey.c $(CORE)/utils.c -o tools/intkey_gen

condi_gen:
	$(CC) $(INCLUDES) tools/condi_gen.c $(POLY)/condi.c $(CORE)/utils.c -o tools/condi_gen

# Standalone test-data generator for Slidefair ciphers. Reuses the real cipher code
# (slidefair.c + utils.c) so the generator and the solver can never drift in convention.
slidefair_gen:
	$(CC) $(INCLUDES) tools/slidefair_gen.c $(GRAPH)/slidefair.c $(CORE)/utils.c -o tools/slidefair_gen

# Standalone test-data generator for Seriated Playfair ciphers. Reuses the real cipher
# code (seriated_playfair.c + playfair.c keyword/grid build + utils.c) so the generator
# and the solver can never drift in convention.
seriated_playfair_gen:
	$(CC) $(INCLUDES) tools/seriated_playfair_gen.c $(GRAPH)/seriated_playfair.c $(GRAPH)/playfair.c $(CORE)/utils.c -o tools/seriated_playfair_gen

# Standalone test-data generator for Digrafid ciphers. Reuses the real cipher code
# (digrafid.c + bifid.c build-inverse + utils.c) so the generator and the solver can never
# drift in convention.
digrafid_gen:
	$(CC) $(INCLUDES) tools/digrafid_gen.c $(GRAPH)/digrafid.c $(GRAPH)/bifid.c $(CORE)/utils.c -o tools/digrafid_gen

# Standalone test-data generator for CM Bifid ciphers. Reuses the real cipher code
# (cm_bifid.c + bifid.c build-inverse/keyword + utils.c) so the generator and the solver
# can never drift in convention.
cm_bifid_gen:
	$(CC) $(INCLUDES) tools/cm_bifid_gen.c $(GRAPH)/cm_bifid.c $(GRAPH)/bifid.c $(CORE)/utils.c -o tools/cm_bifid_gen

# Standalone test-data generator for Fractionated Morse ciphers. Reuses the real cipher code
# (fracmorse.c + utils.c) so the generator and the solver can never drift in convention.
fracmorse_gen:
	$(CC) $(INCLUDES) tools/fracmorse_gen.c $(GRAPH)/fracmorse.c $(CORE)/utils.c -o tools/fracmorse_gen

# Period column order generator (links the real primitive; carries spaces as grid cells).
period_column_gen:
	$(CC) $(INCLUDES) tools/period_column_gen.c $(TRANS)/period_column.c -o tools/period_column_gen

# Standalone test-data generator for double columnar transposition ciphers. Reuses the
# real columnar primitive (transpositions.c) so the generator and the solver can never
# drift in convention.
double_transposition_gen:
	$(CC) $(INCLUDES) tools/double_transposition_gen.c $(TRANS)/transpositions.c $(CORE)/utils.c -o tools/double_transposition_gen

clean:
	rm -f colossus tests/test_transpositions tests/test_ciphers tests/test_optimal_cycleword tests/test_held_karp tests/test_solver tests/test_playfair tests/test_playfair_solver tests/test_bifid tests/test_bifid_solver tests/test_trifid tests/test_trifid_solver tests/test_hill tests/test_hill_solver tests/test_gronsfeld tests/test_gronsfeld_solver tests/test_phillips tests/test_phillips_solver tests/test_twosquare tests/test_twosquare_solver tests/test_foursquare tests/test_foursquare_solver tests/test_adfgvx tests/test_adfgvx_solver tests/test_nihilist_sub tests/test_nihilist_sub_solver tests/test_gromark tests/test_gromark_solver tests/test_nicodemus tests/test_nicodemus_solver tests/test_bazeries tests/test_bazeries_solver tests/test_portax tests/test_portax_solver tests/test_progkey tests/test_progkey_solver tests/test_slidefair tests/test_slidefair_solver tests/test_seriated_playfair tests/test_seriated_playfair_solver tests/test_digrafid tests/test_digrafid_solver tests/test_cm_bifid tests/test_cm_bifid_solver tests/test_trisquare tests/test_trisquare_solver tests/test_intkey tests/test_intkey_solver tests/test_condi tests/test_condi_solver tests/test_fracmorse tests/test_fracmorse_solver tests/test_period_column tests/test_period_column_solver tests/test_period_column_space_solver tests/test_double_transposition tests/test_double_transposition_solver tools/period_column_gen tools/double_transposition_gen tools/fracmorse_gen tools/trisquare_gen tools/intkey_gen tools/condi_gen tools/homophonic_gen tools/digrafid_gen tools/cm_bifid_gen tools/seriated_playfair_gen tools/slidefair_gen tools/progkey_gen tools/portax_gen tools/bazeries_gen tools/nihilist_sub_gen tools/playfair_gen tools/bifid_gen tools/trifid_gen tools/hill_gen tools/gronsfeld_gen tools/phillips_gen tools/twosquare_gen tools/foursquare_gen tools/adfgvx_gen tools/gromark_gen tools/nicodemus_gen

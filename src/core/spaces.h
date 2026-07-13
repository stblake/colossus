// -spaces: post-decrypt readability pass. A recovered plaintext is normally one unbroken
// run of letters (or, for the space-significant cipher types, a mix of letters and the
// cipher's own carried-through spaces). This module inserts ADDITIONAL spaces into the
// letters-only runs to make the plaintext readable, scored by an independently-loaded
// space-inclusive character n-gram table (see spaces.c for the file format and the exact
// Viterbi segmentation this runs).

#ifndef SPACES_H
#define SPACES_H
#include "colossus.h"     // SpacesNgramTable is forward-declared there (shared with every solver)

// Loads a space-inclusive character n-gram table from `filename`: each line is
// "<order-char run of A-Z and ' '> <frequency>" (frequency separated by the RIGHTMOST
// whitespace, since the n-gram itself may contain embedded spaces). Builds a dense
// log10-probability table over the 27-symbol {A..Z, ' '} alphabet with an unseen-window
// floor (the same convention as load_ngrams()'s -logprob mode). Returns NULL on failure
// (bad file, order out of range, or the table would be too large to allocate).
SpacesNgramTable *load_spaces_ngrams(const char *filename, int order, bool verbose);

void free_spaces_ngrams(SpacesNgramTable *tbl);

// Runs an exact Viterbi word segmentation over indices[0..len-1] (the standard colossus
// convention: indices[i] >= 0 is a live plaintext alphabet symbol, negative is a passthrough
// sentinel already carrying a literal character -- e.g. a cipher-native space/punctuation
// mark). Existing passthrough characters are treated as fixed hard boundaries; only maximal
// runs of live symbols are segmented. Returns a malloc'd NUL-terminated string (caller frees)
// that reproduces what print_text() would print, with additional single spaces inserted at
// the positions the n-gram table favours. NULL if tbl is NULL.
char *spaces_insert(const SpacesNgramTable *tbl, int indices[], int len);

// If tbl is non-NULL, runs spaces_insert() and prints the result as a labelled report line
// (no-op, prints nothing, if tbl is NULL). The single call site every solver's final report
// should use.
void print_spaces_line(const SpacesNgramTable *tbl, int indices[], int len);

#endif

#include "colossus.h"
#include <limits.h>

// =====================================================================
//  Period column order transposition primitive (PERIOD_COLUMN)
// =====================================================================
//
//  A faithful port of AZdecrypt's "Period column order" transposition (its
//  case 13 / the "Period column order" optimize block). It is a COLUMNAR
//  transposition whose column permutation is fixed not by a keyword but by a
//  PERIOD p: the dx grid columns are visited in the periodic order
//
//      0, p, 2p, ...   then   1, 1+p, ...   then   2, 2+p, ...   (< dx)
//
//  and the k-th visited column receives (utp==0) / supplies (utp==1) the
//  natural k-th column. The message is laid row-major into a dx-wide grid of
//  dy = ceil(len/dx) rows, the columns are permuted, and the grid is read back
//  row-major. Padding cells (the dx*dy - len trailing cells of the last row)
//  are marked empty and skipped on readout, so the transform is exactly
//  length-preserving. Composing several stages reproduces AZdecrypt's stacked
//  "Period column order" solutions (the motivating 168-char cipher was a
//  two-stage composition, 4x42 then 56x3).

// Padding marker. Real cipher values are letters 0..g_alpha-1 (>= 0) or negative
// space/punctuation sentinels from ord() (small negatives, ~[-256,-1]); INT_MIN is
// far outside both ranges, so a padding cell is unambiguous and a real space rides
// through the transposition like any letter.
#define PC_EMPTY INT_MIN

// Grid scratch. dx*dy < len + dx <= 2*len for dx <= len, so 2*MAX_CIPHER_LENGTH
// bounds any grid this solver builds. File-static (the whole engine is single-threaded).
static int pc_gr0[2 * MAX_CIPHER_LENGTH];
static int pc_gr1[2 * MAX_CIPHER_LENGTH];

void period_column_transform(const int *in, int *out, int len, int dx, int p, int utp) {
    if (dx < 1) dx = 1;
    if (p < 1) p = 1;
    int dy = (len + dx - 1) / dx;             // ceil(len / dx)
    int cells = dx * dy;

    // Lay the input row-major; cell (x, y) lives at gr0[y*dx + x] = in[y*dx + x].
    // The trailing cells.tail (indices >= len) are padding.
    for (int i = 0; i < cells; i++) {
        pc_gr0[i] = (i < len) ? in[i] : PC_EMPTY;
        pc_gr1[i] = PC_EMPTY;
    }

    // Permute columns by the periodic visiting order. k counts visited columns
    // 0,1,2,...; the k-th visited column index is j. utp==0 places source column
    // k at destination column j (transpose); utp==1 is the inverse (dest k <- src j).
    int k = 0;
    for (int i = 0; i < p; i++) {
        for (int j = i; j < dx; j += p) {
            for (int y = 0; y < dy; y++) {
                if (utp == 0) pc_gr1[y * dx + j] = pc_gr0[y * dx + k];
                else          pc_gr1[y * dx + k] = pc_gr0[y * dx + j];
            }
            k++;
        }
    }

    // Read back row-major, dropping the padding cells; exactly len values emerge.
    int a = 0;
    for (int i = 0; i < cells; i++)
        if (pc_gr1[i] != PC_EMPTY) out[a++] = pc_gr1[i];
}

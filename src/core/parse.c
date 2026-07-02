
// Parser for cipher types


#include "colossus.h"


// Helper for case-insensitive comparison
int str_eq(const char *a, const char *b) {
    return strcasecmp(a, b) == 0;
}



int parse_cipher_type(const char *arg) {
    // Check if the argument is a pure integer.
    char *endptr;
    long val = strtol(arg, &endptr, 10);
    
    // If endptr points to the null terminator, the whole string was a number.
    if (*arg != '\0' && *endptr == '\0') {
        return (int)val;
    }

    // Check string aliases (case insensitive.)
    
    // Vigenere
    if (str_eq(arg, "vig") || str_eq(arg, "vigenere")) return VIGENERE;

    // Quagmire I
    if (str_eq(arg, "q1") || str_eq(arg, "quag1") || str_eq(arg, "quagmire1")) return QUAGMIRE_1;

    // Quagmire II
    if (str_eq(arg, "q2") || str_eq(arg, "quag2") || str_eq(arg, "quagmire2")) return QUAGMIRE_2;

    // Quagmire III
    if (str_eq(arg, "q3") || str_eq(arg, "quag3") || str_eq(arg, "quagmire3")) return QUAGMIRE_3;

    // Quagmire IV
    if (str_eq(arg, "q4") || str_eq(arg, "quag4") || str_eq(arg, "quagmire4")) return QUAGMIRE_4;

    // Beaufort
    if (str_eq(arg, "beau") || str_eq(arg, "beaufort")) return BEAUFORT;

    // Porta
    if (str_eq(arg, "porta")) return PORTA;

    // Autokey (Vigenere Tableau)
    if (str_eq(arg, "auto") || str_eq(arg, "autokey") || str_eq(arg, "auto0") || str_eq(arg, "autovig")) return AUTOKEY_0;

    // Autokey Variants
    if (str_eq(arg, "auto1") || str_eq(arg, "autokey1") || str_eq(arg, "autoquag1")) return AUTOKEY_1;
    if (str_eq(arg, "auto2") || str_eq(arg, "autokey2") || str_eq(arg, "autoquag2")) return AUTOKEY_2;
    if (str_eq(arg, "auto3") || str_eq(arg, "autokey3") || str_eq(arg, "autoquag3")) return AUTOKEY_3;
    if (str_eq(arg, "auto4") || str_eq(arg, "autokey4") || str_eq(arg, "autoquag4")) return AUTOKEY_4;

    // Autokey with Beaufort and Porta tableau
	if (str_eq(arg, "auto5") || str_eq(arg, "abeau") || str_eq(arg, "autobeau") || str_eq(arg, "autobeaufort")) return AUTOKEY_BEAU;
	if (str_eq(arg, "auto6") || str_eq(arg, "aporta") || str_eq(arg, "autoporta")) return AUTOKEY_PORTA;

    // Transposition ciphers (solved by optimization over the transform parameters)
    if (str_eq(arg, "transmatrix") || str_eq(arg, "tmatrix") || str_eq(arg, "matrix")) return TRANSMATRIX;
    if (str_eq(arg, "transperoffset") || str_eq(arg, "transperiodoffset") ||
        str_eq(arg, "transperoff") || str_eq(arg, "tpo")) return TRANSPEROFFSET;
    if (str_eq(arg, "trans") || str_eq(arg, "transpo") || str_eq(arg, "transposition")) return TRANSPOSITION;

    // Columnar transposition (dedicated solver over the column-order permutation)
    if (str_eq(arg, "transcol") || str_eq(arg, "transpocolumnar") ||
        str_eq(arg, "columnar") || str_eq(arg, "col")) return TRANSCOL;
    if (str_eq(arg, "transcol2") || str_eq(arg, "doublecolumnar") ||
        str_eq(arg, "doublecol") || str_eq(arg, "dcol")) return TRANSCOL2;

    // Columnar with within-column track permutation L (exact seam best-L)
    if (str_eq(arg, "transcol-l") || str_eq(arg, "transcoll") || str_eq(arg, "columnar-track") ||
        str_eq(arg, "coltrack") || str_eq(arg, "transcoltrack")) return TRANSCOL_L;

    // Route + column-key two-stage chain (seam best-L reading)
    if (str_eq(arg, "transroutecol") || str_eq(arg, "routecol") || str_eq(arg, "routecolumnar") ||
        str_eq(arg, "chain")) return TRANSROUTECOL;

    // Sub-grid / tile transposition (uniform h x w tile cell permutation)
    if (str_eq(arg, "transtile") || str_eq(arg, "tile") || str_eq(arg, "subgrid")) return TRANSTILE;

    // Rail fence (covers the variant rail fence via the phase-offset sweep)
    if (str_eq(arg, "railfence") || str_eq(arg, "rail") || str_eq(arg, "rails") ||
        str_eq(arg, "varrailfence")) return RAILFENCE;

    // Route transposition (snake / spiral over an R x C grid)
    if (str_eq(arg, "route") || str_eq(arg, "routetransposition") ||
        str_eq(arg, "routetramp")) return ROUTE;

    // Amsco (alternating 1/2-letter columnar)
    if (str_eq(arg, "amsco")) return AMSCO;

    // Myszkowski (columnar with tied keyword ranks)
    if (str_eq(arg, "myszkowski") || str_eq(arg, "mysz")) return MYSZKOWSKI;

    // Redefence (keyed rail fence)
    if (str_eq(arg, "redefence") || str_eq(arg, "rede")) return REDEFENCE;

    // Cadenus (rotated-column transposition, 25 rows)
    if (str_eq(arg, "cadenus")) return CADENUS;

    // Nihilist transposition (single permutation on rows + columns)
    if (str_eq(arg, "nihilist") || str_eq(arg, "nihilisttransposition") ||
        str_eq(arg, "nihilisttramp")) return NIHILIST;

    // Swagman (N x N Latin-square column transposition)
    if (str_eq(arg, "swagman")) return SWAGMAN;

    // Turning grille
    if (str_eq(arg, "grille") || str_eq(arg, "turninggrille")) return GRILLE;

    // Independent periodic substitution (P independent mixed alphabets)
    if (str_eq(arg, "indep") || str_eq(arg, "indperiodic") || str_eq(arg, "periodicsub")
        || str_eq(arg, "indepperiodic")) return INDEP_PERIODIC;

    // Homophonic substitution (ciphertext alphabet larger than the plaintext alphabet)
    if (str_eq(arg, "homophonic") || str_eq(arg, "homophone") || str_eq(arg, "homo"))
        return HOMOPHONIC;

    // Playfair (digraphic substitution over a 5x5 keyed grid)
    if (str_eq(arg, "playfair") || str_eq(arg, "pf")) return PLAYFAIR;

    // Bifid (Delastelle fractionation over a keyed Polybius square)
    if (str_eq(arg, "bifid") || str_eq(arg, "bf")) return BIFID;

    // Trifid (Delastelle fractionation over a 3x3x3 keyed cube)
    if (str_eq(arg, "trifid") || str_eq(arg, "tf") || str_eq(arg, "tri")) return TRIFID;

    // Hill (polygraphic substitution by a k x k matrix mod 26)
    if (str_eq(arg, "hill")) return HILL;

    // Gronsfeld (Vigenere with a numeric key: per-column shifts 0..9)
    if (str_eq(arg, "gronsfeld") || str_eq(arg, "gron")) return GRONSFELD;

    // Phillips (8-square keyed-Polybius substitution) and its column / row-column variants.
    // Test the variants before the bare "phillips" so the longer aliases win.
    if (str_eq(arg, "phillips-c") || str_eq(arg, "phillipsc") || str_eq(arg, "phillips_c") ||
        str_eq(arg, "phillipscolumn")) return PHILLIPS_C;
    if (str_eq(arg, "phillips-rc") || str_eq(arg, "phillipsrc") || str_eq(arg, "phillips_rc") ||
        str_eq(arg, "phillipsrowcolumn")) return PHILLIPS_RC;
    if (str_eq(arg, "phillips") || str_eq(arg, "phil")) return PHILLIPS;

    // Two-Square (two keyed 5x5 squares). Test the vertical variant before the bare
    // "twosquare"/"ts" so the longer aliases win; the bare name is the ACA horizontal type.
    if (str_eq(arg, "twosquare-v") || str_eq(arg, "twosquarev") || str_eq(arg, "two-square-v") ||
        str_eq(arg, "2square-v") || str_eq(arg, "2sq-v") || str_eq(arg, "tsv") ||
        str_eq(arg, "twosquarevertical")) return TWO_SQUARE_V;
    if (str_eq(arg, "twosquare") || str_eq(arg, "two-square") || str_eq(arg, "2square") ||
        str_eq(arg, "2sq") || str_eq(arg, "ts")) return TWO_SQUARE;

    // Four-Square (two keyed ciphertext squares + two standard plaintext squares)
    if (str_eq(arg, "foursquare") || str_eq(arg, "four-square") || str_eq(arg, "4square") ||
        str_eq(arg, "4sq") || str_eq(arg, "fs")) return FOUR_SQUARE;

    // Tri-Square (three keyed 5x5 squares; digraph -> trigraph). 'ts' is taken by Two-Square,
    // so use 'trisq'/'3sq'.
    if (str_eq(arg, "trisquare") || str_eq(arg, "tri-square") || str_eq(arg, "3square") ||
        str_eq(arg, "3sq") || str_eq(arg, "trisq")) return TRI_SQUARE;

    // ADFGVX / ADFGX (keyed-square fractionation + keyed columnar transposition)
    if (str_eq(arg, "adfgvx") || str_eq(arg, "adfg")) return ADFGVX;
    if (str_eq(arg, "adfgx")) return ADFGX;

    // Nihilist Substitution (periodic additive over a keyed Polybius square). Distinct from
    // the Nihilist TRANSPOSITION above (alias "nihilist"). Test the convention variants before
    // the bare name so the longer aliases win; the bare name is the ACA carry convention.
    if (str_eq(arg, "nihilist-sub-nc") || str_eq(arg, "nihilistsubnc") ||
        str_eq(arg, "nihilist-sub-nocarry") || str_eq(arg, "nihsub-nc")) return NIHILIST_SUB_NC;
    if (str_eq(arg, "nihilist-sub-m100") || str_eq(arg, "nihilistsubm100") ||
        str_eq(arg, "nihilist-sub-mod100") || str_eq(arg, "nihsub-m100")) return NIHILIST_SUB_M100;
    if (str_eq(arg, "nihilist-sub") || str_eq(arg, "nihilistsub") ||
        str_eq(arg, "nihilistsubstitution") || str_eq(arg, "nihsub")) return NIHILIST_SUB;

    // Gromark (keyed-alphabet substitution + chain-addition running key) and its Periodic
    // variant (+ per-group offset). Test the periodic alias before the bare name.
    if (str_eq(arg, "gromark-periodic") || str_eq(arg, "periodicgromark") ||
        str_eq(arg, "gromark-p") || str_eq(arg, "pgromark")) return GROMARK_PERIODIC;
    if (str_eq(arg, "gromark") || str_eq(arg, "gromark-basic") || str_eq(arg, "gm")) return GROMARK;

    // Nicodemus (periodic substitution + per-block columnar). The substitution variant is
    // part of the type: bare = Vigenere, -v = Variant, -b = Beaufort.
    if (str_eq(arg, "nicodemus-variant") || str_eq(arg, "nicodemus-v") || str_eq(arg, "nicov"))
        return NICODEMUS_VARIANT;
    if (str_eq(arg, "nicodemus-beaufort") || str_eq(arg, "nicodemus-b") || str_eq(arg, "nicob"))
        return NICODEMUS_BEAUFORT;
    if (str_eq(arg, "nicodemus") || str_eq(arg, "nico")) return NICODEMUS;

    // Bazeries (keyed-square substitution + digit-grouped reversal, one number key).
    if (str_eq(arg, "bazeries") || str_eq(arg, "baz")) return BAZERIES;

    // Portax (periodic digraphic Porta; vertical pairs over a Porta slide).
    if (str_eq(arg, "portax") || str_eq(arg, "ptx")) return PORTAX;

    // Slidefair (periodic digraphic Vigenere/Variant/Beaufort). Check the variant/beaufort
    // aliases before the bare slidefair so a substring never shadows them.
    if (str_eq(arg, "slidefair-var") || str_eq(arg, "slidefair-v") || str_eq(arg, "sfv"))
        return SLIDEFAIR_VAR;
    if (str_eq(arg, "slidefair-beau") || str_eq(arg, "slidefair-b") || str_eq(arg, "sfb"))
        return SLIDEFAIR_BEAU;
    if (str_eq(arg, "slidefair") || str_eq(arg, "sf")) return SLIDEFAIR;

    // Seriated Playfair (digraphic Playfair over vertical pairs of a two-row seriated
    // layout; the seriation period is swept). 'sf' is taken by Slidefair, so use 'spf'.
    if (str_eq(arg, "seriated-playfair") || str_eq(arg, "seriatedplayfair") ||
        str_eq(arg, "serpf") || str_eq(arg, "spf")) return SERIATED_PLAYFAIR;

    // Digrafid (digraphic fractionation over two keyed 27-symbol alphabets; period swept).
    if (str_eq(arg, "digrafid") || str_eq(arg, "df") || str_eq(arg, "dgf")) return DIGRAFID;

    // CM Bifid (Conjugated Matrix Bifid: Bifid fractionation over two keyed squares; period swept).
    if (str_eq(arg, "cm-bifid") || str_eq(arg, "cmbifid") || str_eq(arg, "cmb")) return CM_BIFID;

    // Progressive Key (periodic base cipher + per-group constant key drift). Check the
    // variant/beaufort aliases before the bare progkey so a substring never shadows them.
    if (str_eq(arg, "progkey-var") || str_eq(arg, "progkey-v") || str_eq(arg, "pkv"))
        return PROGKEY_VAR;
    if (str_eq(arg, "progkey-beau") || str_eq(arg, "progkey-b") || str_eq(arg, "pkb"))
        return PROGKEY_BEAU;
    if (str_eq(arg, "progkey") || str_eq(arg, "pk") || str_eq(arg, "progressivekey"))
        return PROGKEY;

    // Interrupted Key (periodic keyword that resets to the first key letter at break points).
    // Check the variant/beaufort aliases before the bare intkey so a substring never shadows them.
    if (str_eq(arg, "interrupted-key-var") || str_eq(arg, "intkey-var") || str_eq(arg, "ikv"))
        return INTERRUPTED_KEY_VAR;
    if (str_eq(arg, "interrupted-key-beau") || str_eq(arg, "intkey-beau") || str_eq(arg, "ikb"))
        return INTERRUPTED_KEY_BEAU;
    if (str_eq(arg, "interrupted-key") || str_eq(arg, "interruptedkey") ||
        str_eq(arg, "intkey") || str_eq(arg, "ik"))
        return INTERRUPTED_KEY;

    // Return -1 to indicate invalid/unknown type.
    return -1;
}





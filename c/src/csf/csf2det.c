/* csf/csf2det.c -- CSF -> Slater determinant expansion */
/**
 * @file csf2det.c
 * @brief Expand one CSF into the Slater determinants it is built from.
 *
 * A CSF produced by the DRT walk is a spin eigenfunction; a determinant is
 * not. The expansion
 *
 *     |CSF_I> = sum_k c_k |det_k>
 *
 * depends only on the step vector and on the spin-coupling convention, not on
 * integrals, not on how the CI space was built. It therefore lives here, next
 * to csf_lookup.c, and not in a consumer.
 *
 * ==================================================================
 *  CONVENTIONS -- all fixed here, none inferred from a passing test
 * ==================================================================
 *
 * (1) SPIN COUPLING ORDER. The genealogical coupling of Eq. (5) of the hexpr
 *     paper, which is what the DRT walk and the expression path both mean by
 *     a CSF:
 *
 *         |Phi_k> = ( f_k  Phi_{k-1} )^{Phi_k}
 *
 *     The NEW one-orbital state f_k stands on the LEFT and the accumulated
 *     partial CSF Phi_{k-1} on the RIGHT. The coupling at each orbital is
 *     therefore
 *
 *         (1/2) (x) S_{k-1}  ->  S_k,
 *
 *     not S_{k-1} (x) (1/2).
 *
 *     THIS ORDER IS NOT COSMETIC. Clebsch-Gordan coefficients are not
 *     symmetric under exchange of the two coupled angular momenta:
 *
 *         <j1 m1, j2 m2 | J M> = (-1)^(j1 + j2 - J) <j2 m2, j1 m1 | J M>.
 *
 *     With j1 = 1/2 and j2 = S_{k-1}, the exponent is 0 for a U step
 *     (S_k = S + 1/2) and 1 for a D step (S_k = S - 1/2). Writing the chain
 *     in the other order therefore multiplies the whole CSF by
 *     (-1)^(number of D steps) -- a per-CSF phase that is invisible to norm,
 *     to orthogonality and to Sz, and that shows up only when the expansion
 *     is compared against something outside the CSF basis.
 *
 *     An earlier version of this file used the S (x) 1/2 order. Every
 *     stage-one self-check passed; the matrix comparison against PySCF
 *     (validation/csf2det_vs_pyscf.py) then failed on exactly the CSFs with
 *     an odd number of D steps. The four closed forms below are written in
 *     the Eq. (5) order.
 *
 * (2) CLEBSCH-GORDAN PHASE. Condon-Shortley. Only one family of CG
 *     coefficients is ever needed -- coupling a single spin-1/2 with an
 *     intermediate spin S -- so the four closed forms are written out here
 *     rather than obtained from a general routine. c/wig/ does NOT supply
 *     them: it exposes wigner_delta, wigner_f3l, wigner_6j and wigner_9j,
 *     none of which take magnetic quantum numbers, so no CG can be built
 *     from it.
 *
 *     In the Eq. (5) order, with mu = +/-1/2 the spin of the new orbital and
 *     m the accumulated projection BEFORE this orbital:
 *
 *       <1/2,+1/2; S,m | S+1/2, m+1/2> = +sqrt( (S+m+1) / (2S+1) )
 *       <1/2,-1/2; S,m | S+1/2, m-1/2> = +sqrt( (S-m+1) / (2S+1) )
 *       <1/2,+1/2; S,m | S-1/2, m+1/2> = +sqrt( (S-m)   / (2S+1) )
 *       <1/2,-1/2; S,m | S-1/2, m-1/2> = -sqrt( (S+m)   / (2S+1) )
 *
 *     The minus sign on the last line is the Condon-Shortley phase carried
 *     through the exchange of Eq. (1)'s two factors. It is the only place in
 *     this file where the phase convention enters.
 *
 * (3) DETERMINANT OPERATOR ORDER. Two orderings appear here. The
 *     genealogical product of (2) gives the coefficient in the
 *     ORBITAL-ORDERED ("interleaved") convention
 *
 *         |det>_orb = prod over p ascending of
 *                       [a^dag_{p,alpha}] [a^dag_{p,beta}]
 *
 *     (those present, alpha first), which is the ordering the DRT walk
 *     induces. The PUBLIC API returns coefficients in the BLOCK
 *     ("all alpha then all beta") convention
 *
 *         |det>_blk = (prod_{p asc, alpha occ} a^dag_{p,alpha})
 *                     (prod_{q asc, beta  occ} a^dag_{q,beta }) |vac>
 *
 *     because that is what most determinant CI codes and most
 *     second-quantized qubit mappings use. The two differ by a
 *     determinant-dependent permutation sign, computed in det_block_sign():
 *
 *         sign = (-1)^N,
 *         N = sum over beta-occupied q of #{ alpha-occupied p : p > q }
 *
 *     A consumer that needs the interleaved convention must undo this sign
 *     itself; the rule above is the whole of it.
 *
 *     Note that Eq. (5) builds the CSF as a descending product
 *     a^dag_{Norb} ... a^dag_1 acting on the vacuum, whereas both orderings
 *     above ascend. For a fixed electron count that difference is one fixed
 *     permutation, the same for every CSF of a DRT, so it contributes a
 *     GLOBAL phase and cancels out of every relative sign, of <I|J>, and of
 *     U H_det U^T. It is stated here because a global phase is a convention
 *     even when nothing measurable depends on it.
 *
 * (4) Ms. Passed as 2*Ms in ms_x2. Any |Ms| <= S with the right parity is
 *     accepted, including negative Ms. Ms is NOT assumed equal to S: a
 *     comparison against a determinant program run at nelec=(na,nb) fixes
 *     Ms = (na-nb)/2 on the other side, and that side's choice is not ours.
 *
 * (5) NORMALIZATION. sum_k c_k^2 == 1 for every CSF, at every Ms. Note that
 *     the generator-form CSFs of Eq. (5) are NOT unit-normalized -- Eq. (21)
 *     gives their norm as sqrt([Phi_Norb]) -- and the expression path strips
 *     that factor separately (Eq. (22)). What this file returns is the
 *     unit-normalized CSF.
 *
 * (6) ZERO COEFFICIENTS ARE EMITTED. The determinant list is every
 *     alpha/beta assignment of the open shells compatible with ms_x2, so the
 *     count is exactly C(n_open, n_alpha_open) and follows from the step
 *     vector alone. Some assignments have a vanishing CG product; they are
 *     still emitted, with coefficient 0.0. Dropping them would require a
 *     threshold, and a threshold is one more convention to get wrong.
 *
 * (7) DETERMINANT ORDER. Lexicographically ascending in the set of
 *     open-shell POSITIONS (0-based, in walk order) that carry alpha. This
 *     is a convention of this API and nothing else depends on it; it is
 *     stated so that two calls can be compared entry by entry.
 */
#include "hexpr.h"
#include "drt.h"
#include "hexpr_error.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Orbital index of walk position k.                                   */
/*                                                                     */
/* step[0] is the FIRST spatial orbital: orbitals ascend core ->        */
/* valence -> external along the step vector. Measured, not assumed --  */
/* see test_orbital_direction() in c/tests/test_csf2det.c, which pins   */
/* it through SOCI admitting a doubly occupied last orbital where FOCI  */
/* does not.                                                           */
/* ------------------------------------------------------------------ */
static int orb_of_walk_pos(int k, int norb)
{
    (void)norb;
    return k;
}

/* ------------------------------------------------------------------ */
/* Clebsch-Gordan: couple ONE spin-1/2 (left) with the accumulated      */
/* intermediate spin S (right), in the Eq. (5) order of convention (1). */
/*                                                                      */
/* Spins are doubled throughout: b = 2S, mm = 2m, both taken BEFORE     */
/* this orbital. Substituting the doubling into the four closed forms   */
/* of convention (2):                                                   */
/*                                                                      */
/*   up,   mu=+1/2:  +sqrt( (b+mm+2) / (2*(b+1)) )                      */
/*   up,   mu=-1/2:  +sqrt( (b-mm+2) / (2*(b+1)) )                      */
/*   down, mu=+1/2:  +sqrt( (b-mm)   / (2*(b+1)) )                      */
/*   down, mu=-1/2:  -sqrt( (b+mm)   / (2*(b+1)) )                      */
/*                                                                      */
/* Every numerator is an integer, so nothing rounds before the square   */
/* root. A non-positive numerator means the coupling is impossible      */
/* (|m| > S on one side) and the coefficient is exactly 0.              */
/* ------------------------------------------------------------------ */
/**
 * @brief One genealogical coupling coefficient, in the Eq. (5) order.
 * @param b        2*S before this orbital.
 * @param mm       2*m before this orbital.
 * @param is_up    nonzero for step U (S -> S+1/2), zero for step D.
 * @param is_alpha nonzero if this orbital carries alpha (mu = +1/2).
 * @return the Condon-Shortley CG coefficient; 0.0 for an impossible coupling.
 */
static double cg_half(int b, int mm, int is_up, int is_alpha)
{
    int    num;
    double sign = 1.0;

    if (is_up) {
        num = is_alpha ? (b + mm + 2) : (b - mm + 2);
    } else if (is_alpha) {
        num = b - mm;
    } else {
        num  = b + mm;
        sign = -1.0;              /* <- the Condon-Shortley phase */
    }

    if (num <= 0 || b < 0) return 0.0;
    return sign * sqrt((double)num / (2.0 * (double)(b + 1)));
}

/* ------------------------------------------------------------------ */
/* Permutation sign from orbital order to all-alpha-then-all-beta.      */
/* See convention (3). occ_a/occ_b are norb-length 0/1 arrays.          */
/* ------------------------------------------------------------------ */
static double det_block_sign(const unsigned char *occ_a,
                             const unsigned char *occ_b,
                             int norb)
{
    long n = 0;
    for (int q = 0; q < norb; q++) {
        if (!occ_b[q]) continue;
        for (int p = q + 1; p < norb; p++)
            if (occ_a[p]) n++;
    }
    return (n & 1L) ? -1.0 : 1.0;
}

/* C(n, k) in 64-bit, with overflow reported rather than wrapped. */
static hexpr_status_t binom(int n, int k, long long *out)
{
    if (k < 0 || k > n) { *out = 0; return HEXPR_OK; }
    if (k > n - k) k = n - k;

    long long r = 1;
    for (int i = 1; i <= k; i++) {
        long long num = (long long)(n - k + i);
        if (r > LLONG_MAX / num) return HEXPR_ERR_RANGE;
        r = r * num / (long long)i;   /* exact at every step */
    }
    *out = r;
    return HEXPR_OK;
}

/* ------------------------------------------------------------------ */
/* Shared front end: validate, fetch the step vector, classify the      */
/* orbitals. `fn` is the public entry point's name, used so that        */
/* hexpr_last_error() names the function the caller actually called.    */
/* ------------------------------------------------------------------ */
typedef struct {
    int  norb;
    int  spin_x2;
    int  n_open;
    int  n_alpha_open;
    int *open_pos;        /* walk positions of the open shells, ascending */
    unsigned char *step;  /* norb step codes */
} csf2det_ctx_t;

static void ctx_free(csf2det_ctx_t *c)
{
    free(c->open_pos);
    free(c->step);
    c->open_pos = NULL;
    c->step     = NULL;
}

static hexpr_status_t ctx_build(const char *fn, const hexpr_drt_t *drt,
                                int index, int ms_x2, csf2det_ctx_t *c)
{
    memset(c, 0, sizeof(*c));

    if (!drt) {
        hexpr_set_error("%s: drt is NULL", fn);
        return HEXPR_ERR_INVALID_ARG;
    }

    c->norb    = hexpr_drt_norb(drt);
    c->spin_x2 = hexpr_drt_spin_x2(drt);
    if (c->norb <= 0) {
        hexpr_set_error("%s: drt reports norb=%d", fn, c->norb);
        return HEXPR_ERR_NOT_BUILT;
    }

    /* Ms must have the same parity as 2S and lie within |Ms| <= S. */
    if (ms_x2 > c->spin_x2 || ms_x2 < -c->spin_x2) {
        hexpr_set_error("%s: ms_x2=%d outside |Ms| <= S (spin_x2=%d)",
                        fn, ms_x2, c->spin_x2);
        return HEXPR_ERR_INVALID_ARG;
    }
    if (((ms_x2 - c->spin_x2) & 1) != 0) {
        hexpr_set_error("%s: ms_x2=%d has the wrong parity for spin_x2=%d",
                        fn, ms_x2, c->spin_x2);
        return HEXPR_ERR_INVALID_ARG;
    }

    c->step = (unsigned char *)malloc((size_t)c->norb);
    if (!c->step) {
        hexpr_set_error("%s: out of memory for a %d-byte step vector",
                        fn, c->norb);
        return HEXPR_ERR_OOM;
    }

    /* hexpr_csf_steps sets its own message on failure; do not overwrite it,
     * since it names the actual cause (index out of range, undecodable
     * walk). Only note which entry point was in progress. */
    hexpr_status_t st = hexpr_csf_steps(drt, index, c->step);
    if (st != HEXPR_OK) { ctx_free(c); return st; }

    c->open_pos = (int *)malloc(sizeof(int) * (size_t)c->norb);
    if (!c->open_pos) {
        ctx_free(c);
        hexpr_set_error("%s: out of memory for the open-shell list", fn);
        return HEXPR_ERR_OOM;
    }

    int b = 0;
    for (int k = 0; k < c->norb; k++) {
        switch (c->step[k]) {
        case HEXPR_STEP_E: break;
        case HEXPR_STEP_F: break;
        case HEXPR_STEP_U: c->open_pos[c->n_open++] = k; b += 1; break;
        case HEXPR_STEP_D: c->open_pos[c->n_open++] = k; b -= 1; break;
        default:
            hexpr_set_error("%s: csf %d has step byte %d at orbital %d",
                            fn, index, (int)c->step[k], k);
            ctx_free(c);
            return HEXPR_ERR_INTERNAL;
        }
        if (b < 0) {
            hexpr_set_error("%s: csf %d drives the intermediate spin negative "
                            "at orbital %d", fn, index, k);
            ctx_free(c);
            return HEXPR_ERR_INTERNAL;
        }
    }

    /* The walk must land on the DRT's own spin. This is the same class of
     * check as the construction-time ncsfs guard: the invariant should hold
     * by construction, and the defect-2 work is why it is checked anyway. */
    if (b != c->spin_x2) {
        hexpr_set_error("%s: csf %d walks to spin_x2=%d but the DRT is %d",
                        fn, index, b, c->spin_x2);
        ctx_free(c);
        return HEXPR_ERR_INTERNAL;
    }

    /* n_alpha_open - n_beta_open = ms_x2, n_alpha_open + n_beta_open = n_open */
    if (((c->n_open + ms_x2) & 1) != 0) {
        hexpr_set_error("%s: csf %d has %d open shells, which cannot carry "
                        "ms_x2=%d", fn, index, c->n_open, ms_x2);
        ctx_free(c);
        return HEXPR_ERR_INTERNAL;
    }
    c->n_alpha_open = (c->n_open + ms_x2) / 2;
    if (c->n_alpha_open < 0 || c->n_alpha_open > c->n_open) {
        hexpr_set_error("%s: csf %d would need %d alpha among %d open shells",
                        fn, index, c->n_alpha_open, c->n_open);
        ctx_free(c);
        return HEXPR_ERR_INTERNAL;
    }
    return HEXPR_OK;
}

/* ================================================================== */
/* Public API                                                          */
/* ================================================================== */

/**
 * @brief Number of determinants in the expansion of CSF @p index at @p ms_x2.
 * @param[out] out_ndets receives C(n_open, n_alpha_open) (untouched on error).
 * @return @c HEXPR_OK, @c HEXPR_ERR_INVALID_ARG for a NULL drt or an ms_x2
 *         outside |Ms| <= S or of the wrong parity, @c HEXPR_ERR_RANGE for an
 *         index out of range or a count that exceeds int, or the status of the
 *         underlying hexpr_csf_steps(). @c hexpr_last_error() describes the
 *         failure in every case.
 */
hexpr_status_t hexpr_csf_ndets(const hexpr_drt_t *drt, int index, int ms_x2,
                               int *out_ndets)
{
    static const char FN[] = "hexpr_csf_ndets";

    if (!out_ndets) {
        hexpr_set_error("%s: out_ndets is NULL", FN);
        return HEXPR_ERR_INVALID_ARG;
    }

    csf2det_ctx_t c;
    hexpr_status_t st = ctx_build(FN, drt, index, ms_x2, &c);
    if (st != HEXPR_OK) return st;

    long long nd = 0;
    st = binom(c.n_open, c.n_alpha_open, &nd);
    const int n_open = c.n_open, n_alpha = c.n_alpha_open;
    ctx_free(&c);

    if (st != HEXPR_OK) {
        hexpr_set_error("%s: C(%d, %d) overflows", FN, n_open, n_alpha);
        return st;
    }
    if (nd > (long long)INT_MAX) {
        hexpr_set_error("%s: C(%d, %d) = %lld exceeds int",
                        FN, n_open, n_alpha, nd);
        return HEXPR_ERR_RANGE;
    }

    *out_ndets = (int)nd;
    return HEXPR_OK;
}

/**
 * @brief Expand CSF @p index into determinants at @p ms_x2.
 *
 * Buffer sizes, with ndets from hexpr_csf_ndets() and
 * nwords = (hexpr_drt_norb(drt) + 63) / 64:
 *   out_alpha, out_beta : ndets * nwords uint64_t each
 *   out_coef            : ndets double
 * Determinant k occupies out_alpha[k*nwords .. k*nwords+nwords-1]; bit p of
 * word w is spatial orbital 64*w + p.
 *
 * @param[out] out_count receives the number written (== hexpr_csf_ndets()).
 * @return as hexpr_csf_ndets().
 */
hexpr_status_t hexpr_csf_to_dets(const hexpr_drt_t *drt, int index, int ms_x2,
                                 int *out_count,
                                 uint64_t *out_alpha,
                                 uint64_t *out_beta,
                                 double   *out_coef)
{
    static const char FN[] = "hexpr_csf_to_dets";

    if (!out_count || !out_alpha || !out_beta || !out_coef) {
        hexpr_set_error("%s: out_count/out_alpha/out_beta/out_coef must all "
                        "be non-NULL", FN);
        return HEXPR_ERR_INVALID_ARG;
    }

    csf2det_ctx_t c;
    hexpr_status_t st = ctx_build(FN, drt, index, ms_x2, &c);
    if (st != HEXPR_OK) return st;

    long long nd = 0;
    st = binom(c.n_open, c.n_alpha_open, &nd);
    if (st != HEXPR_OK) {
        hexpr_set_error("%s: C(%d, %d) overflows", FN, c.n_open, c.n_alpha_open);
        ctx_free(&c);
        return st;
    }
    if (nd > (long long)INT_MAX) {
        hexpr_set_error("%s: C(%d, %d) = %lld exceeds int",
                        FN, c.n_open, c.n_alpha_open, nd);
        ctx_free(&c);
        return HEXPR_ERR_RANGE;
    }

    const int ndets  = (int)nd;
    const int norb   = c.norb;
    const int nwords = (norb + 63) / 64;
    const int nopen  = c.n_open;
    const int nalpha = c.n_alpha_open;

    /* comb[] holds the open-shell POSITIONS (indices into open_pos) that
     * carry alpha, ascending. Lexicographic successor gives convention (7). */
    int *comb  = (int *)malloc(sizeof(int) * (size_t)(nalpha > 0 ? nalpha : 1));
    unsigned char *is_a  = (unsigned char *)malloc((size_t)(nopen > 0 ? nopen : 1));
    unsigned char *occ_a = (unsigned char *)calloc((size_t)norb, 1);
    unsigned char *occ_b = (unsigned char *)calloc((size_t)norb, 1);
    if (!comb || !is_a || !occ_a || !occ_b) {
        free(comb); free(is_a); free(occ_a); free(occ_b);
        ctx_free(&c);
        hexpr_set_error("%s: out of memory for %d determinants over %d orbitals",
                        FN, ndets, norb);
        return HEXPR_ERR_OOM;
    }

    for (int i = 0; i < nalpha; i++) comb[i] = i;

    for (int d = 0; d < ndets; d++) {
        memset(is_a, 0, (size_t)(nopen > 0 ? nopen : 1));
        for (int i = 0; i < nalpha; i++) is_a[comb[i]] = 1;

        /* --- occupations, in orbital order --- */
        memset(occ_a, 0, (size_t)norb);
        memset(occ_b, 0, (size_t)norb);
        int open_seen = 0;
        for (int k = 0; k < norb; k++) {
            const int p = orb_of_walk_pos(k, norb);
            if (c.step[k] == HEXPR_STEP_F) {
                occ_a[p] = 1;
                occ_b[p] = 1;
            } else if (c.step[k] == HEXPR_STEP_U || c.step[k] == HEXPR_STEP_D) {
                if (is_a[open_seen]) occ_a[p] = 1; else occ_b[p] = 1;
                open_seen++;
            }
        }

        /* --- genealogical coefficient, in the interleaved convention --- */
        double coef = 1.0;
        int b = 0, mm = 0;
        open_seen = 0;
        for (int k = 0; k < norb && coef != 0.0; k++) {
            if (c.step[k] != HEXPR_STEP_U && c.step[k] != HEXPR_STEP_D) continue;
            const int up    = (c.step[k] == HEXPR_STEP_U);
            const int alpha = is_a[open_seen++];
            coef *= cg_half(b, mm, up, alpha);
            b  += up ? 1 : -1;
            mm += alpha ? 1 : -1;
        }

        /* --- to the block convention --- */
        if (coef != 0.0) coef *= det_block_sign(occ_a, occ_b, norb);

        /* --- emit --- */
        uint64_t *da = out_alpha + (size_t)d * (size_t)nwords;
        uint64_t *db = out_beta  + (size_t)d * (size_t)nwords;
        for (int w = 0; w < nwords; w++) { da[w] = 0; db[w] = 0; }
        for (int p = 0; p < norb; p++) {
            if (occ_a[p]) da[p >> 6] |= (uint64_t)1 << (p & 63);
            if (occ_b[p]) db[p >> 6] |= (uint64_t)1 << (p & 63);
        }
        out_coef[d] = coef;

        /* --- next combination, lexicographic --- */
        if (d + 1 < ndets) {
            int i = nalpha - 1;
            while (i >= 0 && comb[i] == nopen - nalpha + i) i--;
            if (i < 0) {           /* cannot happen: d+1 < ndets */
                hexpr_set_error("%s: combination enumeration ran out after %d "
                                "of %d determinants", FN, d + 1, ndets);
                free(comb); free(is_a); free(occ_a); free(occ_b);
                ctx_free(&c);
                return HEXPR_ERR_INTERNAL;
            }
            comb[i]++;
            for (int j = i + 1; j < nalpha; j++) comb[j] = comb[j - 1] + 1;
        }
    }

    free(comb); free(is_a); free(occ_a); free(occ_b);
    ctx_free(&c);

    *out_count = ndets;
    return HEXPR_OK;
}

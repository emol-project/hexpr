/* test_drt_from_csfset.c -- unit tests for hexpr_drt_from_csfset */
#include "hexpr.h"
#include "loopgen/drt.h"
#include "csf/csf_lookup.h"   /* internal walk_csf_index / walk_csf_steps */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ------------------------------------------------------------------ */
/* T1: round-trip on full s_full_4o4e CSF set                         */
/* ------------------------------------------------------------------ */

static int t1_round_trip(void)
{
    int ok = 1;

    hexpr_drt_t *soci = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    if (!soci) { printf("  T1: FAIL (soci NULL)\n"); return 1; }

    int norb  = hexpr_drt_norb(soci);
    int ncsfs = hexpr_drt_ncsfs(soci);
    int nn    = hexpr_drt_nnodes(soci);

    /* copy packed CSF buffer before building sub-DRT */
    unsigned char *csf_buf = malloc((size_t)ncsfs * (size_t)norb);
    memcpy(csf_buf, hexpr_drt_csfs(soci), (size_t)ncsfs * (size_t)norb);

    hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, csf_buf, ncsfs);
    free(csf_buf);

    if (!sub) {
        printf("  T1: FAIL (sub NULL: %s)\n", hexpr_last_error());
        hexpr_drt_destroy(soci);
        return 1;
    }

    /* nnodes and norb must match */
    if (hexpr_drt_nnodes(sub) != nn) {
        printf("  T1: nnodes mismatch: sub=%d soci=%d\n",
               hexpr_drt_nnodes(sub), nn);
        ok = 0;
    }
    if (hexpr_drt_norb(sub) != norb) {
        printf("  T1: norb mismatch\n"); ok = 0;
    }

    /* per-node accessor comparison */
    int nn_cmp = (hexpr_drt_nnodes(sub) < nn) ? hexpr_drt_nnodes(sub) : nn;
    for (int i = 0; i < nn_cmp; i++) {
        if (hexpr_drt_orb(sub, i)  != hexpr_drt_orb(soci, i))  { ok = 0; }
        if (hexpr_drt_nve(sub, i)  != hexpr_drt_nve(soci, i))  { ok = 0; }
        if (hexpr_drt_spin(sub, i) != hexpr_drt_spin(soci, i)) { ok = 0; }
        if (hexpr_drt_number_upper_arcs(sub, i) !=
            hexpr_drt_number_upper_arcs(soci, i)) { ok = 0; }
        if (hexpr_drt_number_lower_arcs(sub, i) !=
            hexpr_drt_number_lower_arcs(soci, i)) { ok = 0; }
        for (int w = 0; w < 4; w++) {
            if (hexpr_drt_upper_arc(sub, i, w) !=
                hexpr_drt_upper_arc(soci, i, w))  { ok = 0; }
            if (hexpr_drt_lower_arc(sub, i, w) !=
                hexpr_drt_lower_arc(soci, i, w))  { ok = 0; }
            if (hexpr_drt_arc_weight(sub, i, w) !=
                hexpr_drt_arc_weight(soci, i, w)) { ok = 0; }
        }

        int n1, n2;
        const int *lp1 = hexpr_drt_lower_path_weights(soci, i, &n1);
        const int *lp2 = hexpr_drt_lower_path_weights(sub,  i, &n2);
        if (n1 != n2) { ok = 0; }
        else for (int k = 0; k < n1; k++) if (lp1[k] != lp2[k]) { ok = 0; }

        const int *up1 = hexpr_drt_upper_path_weights(soci, i, &n1);
        const int *up2 = hexpr_drt_upper_path_weights(sub,  i, &n2);
        if (n1 != n2) { ok = 0; }
        else for (int k = 0; k < n1; k++) if (up1[k] != up2[k]) { ok = 0; }
    }
    if (!ok) printf("  T1: per-node accessor mismatch detected\n");

    /* packed CSF buffers must be identical */
    int ncsfs2 = hexpr_drt_ncsfs(sub);
    if (ncsfs2 != ncsfs) {
        printf("  T1: ncsfs mismatch: sub=%d soci=%d\n", ncsfs2, ncsfs);
        ok = 0;
    } else if (memcmp(hexpr_drt_csfs(sub), hexpr_drt_csfs(soci),
                      (size_t)ncsfs * (size_t)norb) != 0) {
        printf("  T1: CSF buffer mismatch\n");
        ok = 0;
    }

    hexpr_drt_destroy(soci);
    hexpr_drt_destroy(sub);

    printf("  T1 round-trip s_full_4o4e: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* T2: structural invariants on a 2-CSF subset                        */
/* ------------------------------------------------------------------ */

static int t2_two_csf_subset(void)
{
    int ok = 1;

    hexpr_drt_t *soci = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    if (!soci) { printf("  T2: FAIL (soci NULL)\n"); return 1; }

    int norb    = hexpr_drt_norb(soci);
    int nn_soci = hexpr_drt_nnodes(soci);

    unsigned char subset[8];   /* first 2 CSFs, norb=4 */
    memcpy(subset, hexpr_drt_csfs(soci), 2 * (size_t)norb);
    hexpr_drt_destroy(soci);

    hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, subset, 2);
    if (!sub) {
        printf("  T2: FAIL (sub NULL: %s)\n", hexpr_last_error());
        return 1;
    }

    int nn_sub     = hexpr_drt_nnodes(sub);
    int n_sub_csfs = hexpr_drt_ncsfs(sub);        /* user-facing: input count N */
    int n_walk     = hexpr_drt_walk_ncsfs(sub);   /* canonical count M */

    /* sub-DRT nnode must be strictly smaller than full DRT */
    if (nn_sub >= nn_soci) {
        printf("  T2: nnode %d >= soci nnode %d\n", nn_sub, nn_soci);
        ok = 0;
    }

    /* user-facing csfs()/ncsfs() speak the input set: exactly the 2 rows
     * supplied, in input order; canonical is a superset (M >= N). */
    if (n_sub_csfs != 2) {
        printf("  T2: ncsfs (input count) = %d, want 2\n", n_sub_csfs);
        ok = 0;
    }
    if (n_walk < n_sub_csfs) {
        printf("  T2: walk_ncsfs %d < ncsfs %d\n", n_walk, n_sub_csfs);
        ok = 0;
    }
    const unsigned char *uc = hexpr_drt_csfs(sub);
    if (memcmp(uc, subset, (size_t)2 * (size_t)norb) != 0) {
        printf("  T2: csfs() not equal to input rows in input order\n");
        ok = 0;
    }

    /* both input CSFs must appear in the canonical (walk) CSF list */
    const unsigned char *sc = hexpr_drt_walk_csfs(sub);
    int found0 = 0, found1 = 0;
    for (int i = 0; i < n_walk; i++) {
        if (memcmp(sc + (size_t)i * (size_t)norb, subset,              norb) == 0) found0 = 1;
        if (memcmp(sc + (size_t)i * (size_t)norb, subset + norb, norb) == 0) found1 = 1;
    }
    if (!found0 || !found1) {
        printf("  T2: input CSFs not found in sub-DRT (found0=%d found1=%d)\n",
               found0, found1);
        ok = 0;
    }

    /* arc consistency: upper_arc[inode][w]=j  =>  lower_arc[j][w]=inode */
    for (int inode = 0; inode < nn_sub; inode++) {
        for (int w = 0; w < 4; w++) {
            int j = hexpr_drt_upper_arc(sub, inode, w);
            if (j >= 0 && hexpr_drt_lower_arc(sub, j, w) != inode) {
                printf("  T2: arc inconsistency at inode=%d w=%d j=%d\n", inode, w, j);
                ok = 0;
            }
        }
    }

    /* lower_path_weights at top node has size == canonical (walk) ncsfs */
    int top_idx = nn_sub - 1;
    int n_lpw;
    hexpr_drt_lower_path_weights(sub, top_idx, &n_lpw);
    if (n_lpw != n_walk) {
        printf("  T2: lower_path_weights at top: size=%d walk_ncsfs=%d\n",
               n_lpw, n_walk);
        ok = 0;
    }

    hexpr_drt_destroy(sub);
    printf("  T2 2-CSF subset: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* T3: input contract enforcement                                      */
/* ------------------------------------------------------------------ */

#define T3_CHECK_NULL(label, norb_, buf_, n_)                           \
    do {                                                                \
        hexpr_drt_t *_d = hexpr_drt_from_csfset((norb_), (buf_), (n_)); \
        if (_d) {                                                       \
            printf("  T3: %s should return NULL\n", (label)); ok = 0; \
            hexpr_drt_destroy(_d);                                     \
        } else if (!hexpr_last_error() || !hexpr_last_error()[0]) {    \
            printf("  T3: %s should set last_error\n", (label)); ok = 0; \
        }                                                               \
    } while (0)

static int t3_contracts(void)
{
    int ok = 1;

    /* NULL csfs pointer */
    {
        hexpr_drt_t *d = hexpr_drt_from_csfset(4, NULL, 1);
        if (d) { printf("  T3: NULL csfs should return NULL\n"); ok = 0; hexpr_drt_destroy(d); }
    }

    /* n_csfs = 0 */
    {
        unsigned char buf[4] = {0,0,3,3};
        hexpr_drt_t *d = hexpr_drt_from_csfset(4, buf, 0);
        if (d) { printf("  T3: n_csfs=0 should return NULL\n"); ok = 0; hexpr_drt_destroy(d); }
    }

    /* step out of range (4) */
    {
        unsigned char buf[4] = {4,0,0,0};
        T3_CHECK_NULL("step=4 out of range", 4, buf, 1);
    }

    /* is < 0: step D before any U -> is=-1 at orb=1 */
    {
        unsigned char buf[4] = {2,2,0,0};   /* D,D,E,E */
        T3_CHECK_NULL("is<0 from D first", 4, buf, 1);
    }

    /* different top nodes:
     * CSF0 {3,0,0,0} = feee -> top (4,2,0)
     * CSF1 {1,1,0,0} = uuee -> top (4,2,2) */
    {
        unsigned char buf[8] = {3,0,0,0, 1,1,0,0};
        T3_CHECK_NULL("different top nodes", 4, buf, 2);
    }

    printf("  T3 input contracts: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* T4: corner cases                                                    */
/* ------------------------------------------------------------------ */

static int t4_corner_cases(void)
{
    int ok = 1;

    /* T4a: single CSF "eeff" -> linear chain of 5 nodes */
    /* eeff={0,0,3,3}: (0,0,0)->E->(1,0,0)->E->(2,0,0)->F->(3,2,0)->F->(4,4,0) */
    {
        unsigned char csf[4] = {0,0,3,3};
        hexpr_drt_t *d = hexpr_drt_from_csfset(4, csf, 1);
        if (!d) {
            printf("  T4a: 1-CSF returned NULL\n"); ok = 0;
        } else {
            if (hexpr_drt_nnodes(d) != 5) {
                printf("  T4a: expected nnode=5, got %d\n", hexpr_drt_nnodes(d)); ok = 0;
            }
            if (hexpr_drt_ncsfs(d) != 1) {
                printf("  T4a: expected ncsfs=1, got %d\n", hexpr_drt_ncsfs(d)); ok = 0;
            }
            hexpr_drt_destroy(d);
        }
    }

    /* T4b: duplicate CSF -> canonical walk dedups to 1 path (5 nodes), but
     * the user-facing input view retains both supplied rows (ncsfs == N == 2).
     * Demonstrates walk_ncsfs (canonical) can be < ncsfs (input) under dedup. */
    {
        unsigned char csfs[8] = {0,0,3,3, 0,0,3,3};
        hexpr_drt_t *d = hexpr_drt_from_csfset(4, csfs, 2);
        if (!d) {
            printf("  T4b: dup CSF returned NULL\n"); ok = 0;
        } else {
            if (hexpr_drt_nnodes(d) != 5) {
                printf("  T4b: expected nnode=5, got %d\n", hexpr_drt_nnodes(d)); ok = 0;
            }
            if (hexpr_drt_ncsfs(d) != 2) {   /* input count, dups retained */
                printf("  T4b: expected ncsfs(input)=2, got %d\n", hexpr_drt_ncsfs(d)); ok = 0;
            }
            if (hexpr_drt_walk_ncsfs(d) != 1) {  /* canonical dedup */
                printf("  T4b: expected walk_ncsfs=1, got %d\n", hexpr_drt_walk_ncsfs(d)); ok = 0;
            }
            hexpr_drt_destroy(d);
        }
    }

    /* T4c: no shared intermediate nodes
     * eeff={0,0,3,3}: (0,0,0)(1,0,0)(2,0,0)(3,2,0)(4,4,0) -- 5 nodes
     * ffee={3,3,0,0}: (0,0,0)(1,2,0)(2,4,0)(3,4,0)(4,4,0) -- 5 nodes
     * shared: (0,0,0) and (4,4,0) -> 8 unique nodes */
    {
        unsigned char csfs[8] = {0,0,3,3, 3,3,0,0};
        hexpr_drt_t *d = hexpr_drt_from_csfset(4, csfs, 2);
        if (!d) {
            printf("  T4c: no-share returned NULL\n"); ok = 0;
        } else {
            if (hexpr_drt_nnodes(d) != 8) {
                printf("  T4c: expected nnode=8, got %d\n", hexpr_drt_nnodes(d)); ok = 0;
            }
            hexpr_drt_destroy(d);
        }
    }

    /* T4d: CSFs that re-merge at an intermediate node
     * eudf={0,1,2,3}: (0,0,0)->E->(1,0,0)->U->(2,1,1)->D->(3,2,0)->F->(4,4,0)
     * uefd={1,0,3,2}: (0,0,0)->U->(1,1,1)->E->(2,1,1)->F->(3,3,1)->D->(4,4,0)
     * shared: (0,0,0),(2,1,1),(4,4,0) -> 7 unique nodes
     * shared node (2,1,1) has 2 lower paths and 2 upper paths */
    {
        unsigned char csfs[8] = {0,1,2,3, 1,0,3,2};
        hexpr_drt_t *d = hexpr_drt_from_csfset(4, csfs, 2);
        if (!d) {
            printf("  T4d: shared-intermediate returned NULL\n"); ok = 0;
        } else {
            if (hexpr_drt_nnodes(d) != 7) {
                printf("  T4d: expected nnode=7, got %d\n", hexpr_drt_nnodes(d)); ok = 0;
            }
            /* find node (2,1,1) and verify nlower=2, nupper=2 */
            int nn = hexpr_drt_nnodes(d);
            int found = 0;
            for (int i = 0; i < nn; i++) {
                int o, nv, is;
                hexpr_drt_node_info(d, i, &o, &nv, &is);
                if (o == 2 && nv == 1 && is == 1) {
                    int nl = hexpr_drt_number_lower_arcs(d, i);
                    int nu = hexpr_drt_number_upper_arcs(d, i);
                    if (nl != 2 || nu != 2) {
                        printf("  T4d: node(2,1,1) nlower=%d nupper=%d (want 2,2)\n",
                               nl, nu);
                        ok = 0;
                    }
                    found = 1;
                    break;
                }
            }
            if (!found) {
                printf("  T4d: node (2,1,1) not found\n"); ok = 0;
            }
            hexpr_drt_destroy(d);
        }
    }

    printf("  T4 corner cases: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* T5: non-canonical input order + superset -> input-order accessors    */
/* ------------------------------------------------------------------ */

/* Input (NON-canonical order): row0 = uefd {1,0,3,2}, row1 = eudf {0,1,2,3}.
 * These share node (2,1,1), so the canonical GUGA walk expands to a SUPERSET
 * of M=4: {eudf, eufd, uedf, uefd} (canonical indices 0..3). The two
 * supplied rows are canonical 3 (uefd) and 0 (eudf) -- deliberately reversed
 * from canonical order so input order != canonical order.
 *
 * This is the case that separates the two index spaces, so it is also where
 * the eval API's input-order contract is pinned: indices 0..1 evaluate the
 * two supplied rows, and the canonical indices (3 and 0) are NOT accepted as
 * a way to reach them. */
static int t5_input_order_superset(void)
{
    int ok = 1;
    const int norb = 4;
    unsigned char in[8] = {1,0,3,2,  0,1,2,3};   /* uefd, eudf */
    unsigned char uefd[4] = {1,0,3,2};
    unsigned char eudf[4] = {0,1,2,3};
    unsigned char eufd[4] = {0,1,3,2};   /* superset-only (canonical 1) */
    unsigned char uedf[4] = {1,0,2,3};   /* superset-only (canonical 2) */

    hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, in, 2);
    if (!sub) { printf("  T5: FAIL (sub NULL: %s)\n", hexpr_last_error()); return 1; }

    /* counts: user-facing N=2, canonical superset M=4 */
    if (hexpr_drt_ncsfs(sub) != 2) {
        printf("  T5: ncsfs=%d, want 2 (input count)\n", hexpr_drt_ncsfs(sub)); ok = 0; }
    if (hexpr_drt_walk_ncsfs(sub) != 4) {
        printf("  T5: walk_ncsfs=%d, want 4 (canonical superset)\n",
               hexpr_drt_walk_ncsfs(sub)); ok = 0; }

    /* csfs() == the 2 input rows, in input order */
    const unsigned char *uc = hexpr_drt_csfs(sub);
    if (memcmp(uc, in, (size_t)2 * norb) != 0) {
        printf("  T5: csfs() != input rows in input order\n"); ok = 0; }

    /* csf_steps(i) == input row i */
    unsigned char buf[4];
    if (hexpr_csf_steps(sub, 0, buf) != HEXPR_OK || memcmp(buf, uefd, norb) != 0) {
        printf("  T5: csf_steps(0) != uefd\n"); ok = 0; }
    if (hexpr_csf_steps(sub, 1, buf) != HEXPR_OK || memcmp(buf, eudf, norb) != 0) {
        printf("  T5: csf_steps(1) != eudf\n"); ok = 0; }
    if (hexpr_csf_steps(sub, 2, buf) != HEXPR_ERR_RANGE) {
        printf("  T5: csf_steps(2) should be RANGE (only 2 input rows)\n"); ok = 0; }

    /* csf_index(input row i) == i (input position) */
    int idx = -1;
    if (hexpr_csf_index(sub, uefd, &idx) != HEXPR_OK || idx != 0) {
        printf("  T5: csf_index(uefd)=%d, want 0\n", idx); ok = 0; }
    if (hexpr_csf_index(sub, eudf, &idx) != HEXPR_OK || idx != 1) {
        printf("  T5: csf_index(eudf)=%d, want 1\n", idx); ok = 0; }

    /* csf_index of a superset-only (valid but not input) CSF -> RANGE */
    if (hexpr_csf_index(sub, eufd, &idx) != HEXPR_ERR_RANGE) {
        printf("  T5: csf_index(eufd) should be RANGE (superset-only)\n"); ok = 0; }
    if (hexpr_csf_index(sub, uedf, &idx) != HEXPR_ERR_RANGE) {
        printf("  T5: csf_index(uedf) should be RANGE (superset-only)\n"); ok = 0; }

    /* walk_csfs() is the canonical superset and contains the superset-only ones */
    const unsigned char *wc = hexpr_drt_walk_csfs(sub);
    int found_eufd = 0, found_uedf = 0;
    for (int i = 0; i < 4; i++) {
        if (memcmp(wc + (size_t)i * norb, eufd, norb) == 0) found_eufd = 1;
        if (memcmp(wc + (size_t)i * norb, uedf, norb) == 0) found_uedf = 1;
    }
    if (!found_eufd || !found_uedf) {
        printf("  T5: walk_csfs() missing superset members (eufd=%d uedf=%d)\n",
               found_eufd, found_uedf); ok = 0; }

    /* Eval (requirement 3): evaluate the 2 input rows. The eval API takes
     * INPUT-ORDER indices, so this is simply 0..ncsfs-1 -- no mapping. The
     * canonical indices of these rows are 3 and 0, and passing THOSE is now
     * an error: 3 is outside [0, ncsfs). That is checked below.
     * bra_pos/ket_pos must index the array we pass (values in 0..1),
     * agreeing with csfs() row order. */
    int idx2[2] = { 0, 1 };
    int mt = 0;
    hexpr_subspace_max_terms(sub, 2, HEXPR_EXPR_FULL, &mt);
    if (mt < 1) mt = 1;
    int32_t *bp = malloc((size_t)mt * sizeof(int32_t));
    int32_t *kp = malloc((size_t)mt * sizeof(int32_t));
    int32_t *ij = malloc((size_t)mt * sizeof(int32_t));
    int32_t *pq = malloc((size_t)mt * sizeof(int32_t));
    double  *co = malloc((size_t)mt * sizeof(double));
    int cnt = 0;
    hexpr_status_t est = hexpr_subspace_eval(sub, idx2, 2, HEXPR_EXPR_FULL,
                                             &cnt, bp, kp, ij, pq, co);
    if (est != HEXPR_OK || cnt <= 0) {
        printf("  T5: subspace_eval failed (st=%d cnt=%d)\n", (int)est, cnt); ok = 0;
    } else {
        for (int k = 0; k < cnt; k++) {
            if (bp[k] < 0 || bp[k] >= 2 || kp[k] < 0 || kp[k] >= 2) {
                printf("  T5: bra_pos/ket_pos out of input range: %d,%d\n",
                       bp[k], kp[k]); ok = 0; break;
            }
        }
    }
    free(bp); free(kp); free(ij); free(pq); free(co);

    /* The eval API is closed over the INPUT set. A canonical index that is
     * not an input position -- here 3, the canonical index of uefd, and 2,
     * one past the last input row -- must be rejected rather than silently
     * evaluating some other CSF. This is the trap the input-order contract
     * exists to remove, so it is checked, not assumed. */
    int cbra = -1;
    if (hexpr_drt_walk_csf_index(sub, uefd, &cbra) != HEXPR_OK || cbra != 3) {
        printf("  T5: walk_csf_index(uefd)=%d, want 3\n", cbra); ok = 0; }
    int bad[2] = { cbra, 0 };
    int bcnt = 0;
    if (hexpr_subspace_count(sub, bad, 2, HEXPR_EXPR_FULL, &bcnt)
        != HEXPR_ERR_RANGE) {
        printf("  T5: canonical index %d should be RANGE (ncsfs=2)\n", cbra);
        ok = 0; }
    int one_past[1] = { 2 };
    if (hexpr_subspace_count(sub, one_past, 1, HEXPR_EXPR_FULL, &bcnt)
        != HEXPR_ERR_RANGE) {
        printf("  T5: index 2 should be RANGE (ncsfs=2)\n"); ok = 0; }
    if (hexpr_pair_count(sub, 0, 2, HEXPR_EXPR_FULL, &bcnt)
        != HEXPR_ERR_RANGE) {
        printf("  T5: pair_count ket=2 should be RANGE (ncsfs=2)\n"); ok = 0; }

    hexpr_drt_destroy(sub);
    printf("  T5 input-order + superset: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* T6: internal canonical lookups (walk_csf_index / walk_csf_steps)    */
/* ------------------------------------------------------------------ */

static int t6_internal_walk_lookups(void)
{
    int ok = 1;
    const int norb = 4;
    unsigned char in[8] = {1,0,3,2,  0,1,2,3};   /* uefd, eudf (same as T5) */
    unsigned char eudf[4] = {0,1,2,3};
    unsigned char eufd[4] = {0,1,3,2};
    unsigned char uedf[4] = {1,0,2,3};
    unsigned char uefd[4] = {1,0,3,2};

    hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, in, 2);
    if (!sub) { printf("  T6: FAIL (sub NULL)\n"); return 1; }

    /* canonical indices (GUGA DFS e<u<d<f): eudf0 eufd1 uedf2 uefd3 */
    struct { unsigned char *s; int ci; } exp[4] = {
        {eudf,0}, {eufd,1}, {uedf,2}, {uefd,3} };
    for (int t = 0; t < 4; t++) {
        int ci = -1;
        if (hexpr_drt_walk_csf_index(sub, exp[t].s, &ci) != HEXPR_OK || ci != exp[t].ci) {
            printf("  T6: walk_csf_index mismatch (got %d want %d)\n", ci, exp[t].ci);
            ok = 0;
        }
    }

    /* round-trip over the full canonical superset (M=4) */
    int m = hexpr_drt_walk_ncsfs(sub);
    for (int i = 0; i < m; i++) {
        unsigned char buf[4];
        int ci = -1;
        if (hexpr_drt_walk_csf_steps(sub, i, buf) != HEXPR_OK) { ok = 0; continue; }
        if (hexpr_drt_walk_csf_index(sub, buf, &ci) != HEXPR_OK || ci != i) {
            printf("  T6: canonical round-trip failed at i=%d (got %d)\n", i, ci);
            ok = 0;
        }
    }

    /* internal (canonical) != public (input) for an input row: eudf is
     * canonical 0 but input position 1. */
    int wc = -1, pub = -1;
    hexpr_drt_walk_csf_index(sub, eudf, &wc);   /* canonical 0 */
    hexpr_csf_index(sub, eudf, &pub);           /* input    1 */
    if (wc != 0 || pub != 1) {
        printf("  T6: expected walk=0/public=1 for eudf, got walk=%d public=%d\n",
               wc, pub); ok = 0;
    }

    hexpr_drt_destroy(sub);
    printf("  T6 internal walk lookups: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
    printf("=== test_drt_from_csfset ===\n");

    int failures = 0;
    failures += t1_round_trip();
    failures += t2_two_csf_subset();
    failures += t3_contracts();
    failures += t4_corner_cases();
    failures += t5_input_order_superset();
    failures += t6_internal_walk_lookups();

    printf("%s (%d failure(s))\n",
           failures == 0 ? "PASS" : "FAIL", failures);
    return failures ? 1 : 0;
}

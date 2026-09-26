/* test_ncsfs_consistency.c -- regression test for the SOCI/FOCI ncsfs
 * over-count (nodes[0].number_upper_arcs can exceed the true number of
 * top-reaching walks when a dead node survives prune_unreachable).
 *
 * Two independent, DRT-count-agnostic assertions per case:
 *
 *   (a) every canonical index in [0, hexpr_drt_ncsfs()) must decode via
 *       the public hexpr_csf_steps(). A dead node makes some indices
 *       undecodable (HEXPR_ERR_INTERNAL) even though hexpr_drt_ncsfs()
 *       claims they exist.
 *
 *   (b) hexpr_pair_eval() on a real CSF's diagonal must match the
 *       hexpr_drt_from_csfset() ground truth for that same CSF built in
 *       isolation. A dead node can make the canonical-index decode used
 *       internally by hexpr_pair_eval() fail on a REAL CSF's index too,
 *       silently returning zero terms instead of the true diagonal
 *       energy expression.
 *
 * Neither assertion hardcodes what the correct ncsfs value is -- both
 * simply check internal consistency, so they remain valid regardless of
 * which fix (or none) is in place.
 */
#include "hexpr.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>

/* Run both assertions for one (nve, spin_x2, ncore, nval, next) SOCI case.
 * Returns the number of failed assertions (0 = both pass). */
static int run_case(const char *label, int nve, int spin_x2,
                     int ncore, int nval, int next)
{
    int fail = 0;
    printf("-- %s: SOCI nve=%d spin_x2=%d ncore=%d nval=%d next=%d --\n",
           label, nve, spin_x2, ncore, nval, next);

    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, nve, spin_x2, ncore, nval, next);
    if (!drt) {
        printf("  hexpr_drt_create returned NULL (%s)\n", hexpr_last_error());
        return 1;
    }

    int ncsfs = hexpr_drt_ncsfs(drt);
    int norb  = hexpr_drt_norb(drt);
    printf("  ncsfs=%d walk_ncsfs=%d norb=%d\n", ncsfs, hexpr_drt_walk_ncsfs(drt), norb);

    /* (a) every canonical index in [0, ncsfs) must decode */
    unsigned char step[64];
    int decodable = 0;
    for (int i = 0; i < ncsfs; i++)
        if (hexpr_csf_steps(drt, i, step) == HEXPR_OK) decodable++;
    if (decodable != ncsfs) {
        printf("  (a) FAIL: only %d of ncsfs=%d canonical indices decode\n",
               decodable, ncsfs);
        fail++;
    } else {
        printf("  (a) PASS: all %d canonical indices decode\n", ncsfs);
    }

    /* (b) pair_eval on a real CSF's diagonal vs. from_csfset ground truth.
     * hexpr_drt_csfs() row 0 is always a real, valid CSF (it is the first
     * entry the packed buffer was actually filled with), regardless of
     * whether ncsfs itself is inflated. */
    unsigned char csf0[64];
    memcpy(csf0, hexpr_drt_csfs(drt), (size_t)norb);

    int idx = -1;
    hexpr_status_t st = hexpr_csf_index(drt, csf0, &idx);
    if (st != HEXPR_OK) {
        printf("  (b) FAIL: hexpr_csf_index on a real CSF returned status=%d\n", st);
        fail++;
    } else {
        int32_t ij[256], pqrs[256]; double coef[256];
        int n = -1;
        hexpr_pair_eval(drt, idx, idx, HEXPR_EXPR_FULL, &n, ij, pqrs, coef);

        hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, csf0, 1);
        int n_gt = -1;
        if (sub) hexpr_pair_eval(sub, 0, 0, HEXPR_EXPR_FULL, &n_gt, ij, pqrs, coef);

        if (!sub || n != n_gt) {
            printf("  (b) FAIL: hexpr_pair_eval(drt,%d,%d) gave %d terms, "
                   "from_csfset ground truth gives %d\n", idx, idx, n, n_gt);
            fail++;
        } else {
            printf("  (b) PASS: hexpr_pair_eval(drt,%d,%d) matches ground truth (%d terms)\n",
                   idx, idx, n);
        }
        if (sub) hexpr_drt_destroy(sub);
    }

    hexpr_drt_destroy(drt);
    return fail;
}

/* Run a case where NO valid top-reaching walk exists at all: the requested
 * nve exceeds what SOCI's <=2-particle (FOCI's <=1-particle) external
 * admission allows relative to an all-empty core+valence reference (with
 * ncore=nval=0 the reference has 0 electrons, so nve itself IS the particle
 * count being requested). hexpr_drt_create must return NULL.
 *
 * This is a different bug from case1/case2's: node_new correctly marks the
 * bottom (0,0,0) inadmissible here (it genuinely cannot reach the target
 * within budget), but nothing downstream was honoring that -- the bottom's
 * own invalidity was never stopping prune_unreachable's forward BFS from
 * radiating "reachable from bottom" out through it anyway, so a spurious
 * top-reaching path survived and hexpr_drt_create built a DRT around it,
 * one orbital short of the real bottom (nodes[0] ends up being the second
 * orbital's node, not (0,0,0)) -- explaining why the old assertions (a)/(b)
 * also catch this shape (see case3/case4 below), but the correct fix makes
 * construction itself fail, not just the decode. */
static int run_reject_case(const char *label, hexpr_drt_kind_t kind, const char *kind_name,
                            int nve, int spin_x2, int ncore, int nval, int next)
{
    printf("-- %s: %s nve=%d spin_x2=%d ncore=%d nval=%d next=%d "
           "(expect: no valid walk -> NULL) --\n",
           label, kind_name, nve, spin_x2, ncore, nval, next);
    hexpr_drt_t *drt = hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next);
    if (drt) {
        printf("  FAIL: hexpr_drt_create succeeded (ncsfs=%d, norb=%d) -- "
               "should have found no valid top-reaching walk\n",
               hexpr_drt_ncsfs(drt), hexpr_drt_norb(drt));
        hexpr_drt_destroy(drt);
        return 1;
    }
    printf("  PASS: correctly rejected (%s)\n", hexpr_last_error());
    return 0;
}

int main(void)
{
    printf("=== test_ncsfs_consistency ===\n");

    int failures = 0;
    /* smallest known combination that triggers the divergence */
    failures += run_case("case1", 6, 0, 1, 1, 2);
    /* smaller sibling: fewer decodable indices lost relative to ncsfs,
     * same (ncore, nval, next) partition */
    failures += run_case("case2", 5, 1, 1, 1, 2);

    /* ncore=0 AND nval=0 ("next-only" DRT, no core, no valence at all): the
     * four combinations found by the norb/nve/spin sweep where nve exceeds
     * SOCI/FOCI's own external particle-admission budget from an all-empty
     * reference. All four must be rejected outright. */
    failures += run_reject_case("case3a", HEXPR_DRT_SOCI, "SOCI", 3, 1, 0, 0, 2);
    failures += run_reject_case("case3b", HEXPR_DRT_FOCI, "FOCI", 3, 1, 0, 0, 2);
    failures += run_reject_case("case4a", HEXPR_DRT_SOCI, "SOCI", 4, 0, 0, 0, 2);
    failures += run_reject_case("case4b", HEXPR_DRT_FOCI, "FOCI", 4, 0, 0, 0, 2);

    printf("%s (%d failure(s))\n", failures == 0 ? "PASS" : "FAIL", failures);
    return failures ? 1 : 0;
}

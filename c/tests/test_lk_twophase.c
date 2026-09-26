/* test_lk_twophase.c -- acceptance test for the two-phase L(k) table.
 *
 * Builds every whole-matrix oracle case via hexpr_expr_build and asserts the
 * three SEARCH-phase invariants after each build:
 *     search misses == 0, wigner_9j calls == 0, matrix_el calls == 0.
 * i.e. warm_lk_table fills exactly the key set the search pass reads, so the
 * search pass is a pure table lookup (no hmlkva / wigner / matrix_el). The pair
 * path (OFF, fill-on-miss) is intentionally NOT asserted here. */
#include <stdio.h>
#include "../include/hexpr.h"

/* instrumentation getters (defined in expr_tree.c) */
long hexpr_lk_search_misses(void);
long hexpr_lk_search_wig_calls(void);
long hexpr_lk_search_mel_calls(void);
long hexpr_lk_table_size(void);

static int failures = 0;
static void check(const char *label, int ok)
{
    printf("  %-52s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

typedef struct { const char *name; int nve, spin, ncore, nval, next; hexpr_drt_kind_t type; } case_def_t;

/* mirrors run_reference.c ALL_CASES (the whole-matrix oracle set) */
static const case_def_t CASES[] = {
    {"s_full_4o4e",   4, 0, 0, 4, 0, HEXPR_DRT_SOCI},
    {"t_full_4o4e",   4, 2, 0, 4, 0, HEXPR_DRT_SOCI},
    {"s_full_6o6e",   6, 0, 0, 6, 0, HEXPR_DRT_SOCI},
    {"t_full_6o6e",   6, 2, 0, 6, 0, HEXPR_DRT_SOCI},
    {"s_soci_small",  6, 0, 1, 4, 1, HEXPR_DRT_SOCI},
    {"s_foci_small",  6, 0, 1, 4, 1, HEXPR_DRT_FOCI},
    {"t_soci_small",  6, 2, 1, 4, 1, HEXPR_DRT_SOCI},
    {"t_foci_small",  6, 2, 1, 4, 1, HEXPR_DRT_FOCI},
    {"s_foci_mid",    6, 0, 2, 2, 2, HEXPR_DRT_FOCI},
    {"s_soci_mid",    6, 0, 2, 2, 2, HEXPR_DRT_SOCI},
    {"s_soci_med",   10, 0, 2, 4, 4, HEXPR_DRT_SOCI},
    {"s_foci_med",   10, 0, 2, 4, 4, HEXPR_DRT_FOCI},
    {"t_soci_med",   10, 2, 2, 4, 4, HEXPR_DRT_SOCI},
    {"t_foci_med",   10, 2, 2, 4, 4, HEXPR_DRT_FOCI},
};

int main(void)
{
    printf("=== two-phase L(k): whole-path search must be lookup-only ===\n");
    hexpr_init();
    for (unsigned i = 0; i < sizeof(CASES)/sizeof(CASES[0]); i++) {
        const case_def_t *c = &CASES[i];
        hexpr_drt_t *drt = hexpr_drt_create(c->type, c->nve, c->spin,
                                            c->ncore, c->nval, c->next);
        if (!drt) { check(c->name, 0); continue; }
        hexpr_expr_t *e = hexpr_expr_build(drt);   /* runs generate/warm -> search */
        long miss = hexpr_lk_search_misses();
        long wig  = hexpr_lk_search_wig_calls();
        long mel  = hexpr_lk_search_mel_calls();
        char label[80];
        snprintf(label, sizeof(label), "%-14s miss=%ld wig=%ld mel=%ld",
                 c->name, miss, wig, mel);
        check(label, miss == 0 && wig == 0 && mel == 0);
        hexpr_expr_destroy(e);
        hexpr_drt_destroy(drt);
    }
    printf("shared table size after all builds: %ld keys\n", hexpr_lk_table_size());
    printf(failures ? "%d TWO-PHASE CHECK(S) FAILED\n" : "ALL TWO-PHASE TESTS PASSED\n",
           failures);
    return failures ? 1 : 0;
}

/* test_drt.c -- DRT regression test: compare DRT output against reference */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "loopgen/drt.h"
#include "fmt.h"
#include "util/compare.h"

#define REF_DIR "../../reference"
#define TMP_DIR "tmp/test_drt"

/* ------------------------------------------------------------------ */

typedef struct {
    const char *name;
    int  nve, spin, ncore, nval, next;
    hexpr_drt_kind_t kind;
} tcase_t;

static const tcase_t CASES[] = {
    {"s_full_4o4e",  4, 0, 0, 4, 0, HEXPR_DRT_SOCI},
    {"t_full_4o4e",  4, 2, 0, 4, 0, HEXPR_DRT_SOCI},
    {"s_full_6o6e",  6, 0, 0, 6, 0, HEXPR_DRT_SOCI},
    {"t_full_6o6e",  6, 2, 0, 6, 0, HEXPR_DRT_SOCI},
    {"s_soci_small", 6, 0, 1, 4, 1, HEXPR_DRT_SOCI},
    {"s_foci_small", 6, 0, 1, 4, 1, HEXPR_DRT_FOCI},
    {"t_soci_small", 6, 2, 1, 4, 1, HEXPR_DRT_SOCI},
    {"t_foci_small", 6, 2, 1, 4, 1, HEXPR_DRT_FOCI},
    {"s_soci_med",  10, 0, 2, 4, 4, HEXPR_DRT_SOCI},
    {"s_foci_med",  10, 0, 2, 4, 4, HEXPR_DRT_FOCI},
    {"t_soci_med",  10, 2, 2, 4, 4, HEXPR_DRT_SOCI},
    {"t_foci_med",  10, 2, 2, 4, 4, HEXPR_DRT_FOCI},
};
static const int NCASES = (int)(sizeof(CASES) / sizeof(CASES[0]));

/* ------------------------------------------------------------------ */

static int run_case(const tcase_t *tc)
{
    hexpr_drt_t *drt = hexpr_drt_create(tc->kind,
                                        tc->nve, tc->spin,
                                        tc->ncore, tc->nval, tc->next);
    if (!drt) {
        fprintf(stderr, "  [%s] hexpr_drt_create returned NULL\n", tc->name);
        return 1;
    }

    /* write to tmp file */
    char out_path[512];
    snprintf(out_path, sizeof(out_path), "%s/%s.drt.txt", TMP_DIR, tc->name);

    FILE *f = fopen(out_path, "w");
    if (!f) {
        perror(out_path);
        hexpr_drt_free(drt);
        return 1;
    }

    hexpr_drt_show(drt, f);

    /* CSF list */
    char **csfs = NULL;
    int ncsf = hexpr_drt_mk_csf_list(drt, &csfs);
    fprintf(f, "\n CSF List \n");
    int norb = hexpr_drt_norb(drt);
    for (int i = 0; i < ncsf; i++) {
        hexpr_print_csf(f, csfs[i], norb);
        free(csfs[i]);
    }
    free(csfs);

    fclose(f);
    hexpr_drt_free(drt);

    /* compare against reference */
    char ref_path[512];
    snprintf(ref_path, sizeof(ref_path), "%s/%s.drt.txt", REF_DIR, tc->name);

    int mismatches = compare_files(ref_path, out_path, 1e-10, 1e-10, 5);
    if (mismatches == 0) {
        printf("  %-20s  PASS\n", tc->name);
        return 0;
    } else {
        printf("  %-20s  FAIL (%d mismatch(es))\n", tc->name, mismatches);
        return 1;
    }
}

/* ------------------------------------------------------------------ */

int main(void)
{
    /* create tmp dir */
    system("mkdir -p " TMP_DIR);

    printf("=== DRT regression test ===\n");
    int failures = 0;
    for (int i = 0; i < NCASES; i++)
        failures += run_case(&CASES[i]);

    printf("%s (%d failure(s))\n",
           failures == 0 ? "PASS" : "FAIL", failures);
    return failures ? 1 : 0;
}

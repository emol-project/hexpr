#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/hexpr.h"

int main(int argc, char **argv) {
    int nve = atoi(argv[1]);
    int spin_x2 = atoi(argv[2]);
    int nval = atoi(argv[4]);

    hexpr_init();
    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, nve, spin_x2, 0, nval, 0);
    int ncsfs = hexpr_drt_ncsfs(drt);
    int nnodes = hexpr_drt_nnodes(drt);
    printf("ncsfs=%d, nnodes=%d\n", ncsfs, nnodes);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    hexpr_expr_t *expr = hexpr_expr_build(drt);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    printf("expr_build (DRT-wide brooks_all): %.3f s\n", elapsed);
    printf("nterms_full: %ld\n", hexpr_expr_nterms_full(expr));

    hexpr_expr_destroy(expr);
    hexpr_drt_destroy(drt);
    return 0;
}

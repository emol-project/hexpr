/* test_fortran_symbols.c -- smoke test that F77 underscore symbols exist */
#include <stdio.h>
#include <dlfcn.h>

#ifdef __APPLE__
#define LIBPATH "../libhexpr.dylib"
#else
#define LIBPATH "../libhexpr.so"
#endif

static int failures = 0;

static void check_sym(void *handle, const char *name)
{
    void *sym = dlsym(handle, name);
    const char *err = dlerror();
    int ok = (sym != NULL && err == NULL);
    printf("  %-50s %s\n", name, ok ? "PASS" : "FAIL");
    if (!ok) {
        printf("    error: %s\n", err ? err : "(sym is NULL)");
        failures++;
    }
}

int main(void)
{
    printf("=== F77 symbol smoke test ===\n");

    dlerror(); /* clear */
    void *h = dlopen(LIBPATH, RTLD_LAZY);
    if (!h) {
        printf("FAIL: cannot dlopen %s: %s\n", LIBPATH, dlerror());
        return 1;
    }

    /* Public API symbols (no underscore) */
    check_sym(h, "hexpr_init");
    check_sym(h, "hexpr_shutdown");
    check_sym(h, "hexpr_version_string");
    check_sym(h, "hexpr_last_error");
    check_sym(h, "hexpr_drt_create");
    check_sym(h, "hexpr_drt_destroy");
    check_sym(h, "hexpr_drt_kind");
    check_sym(h, "hexpr_drt_csfs");
    check_sym(h, "hexpr_drt_node_info");
    check_sym(h, "hexpr_drt_show");
    check_sym(h, "hexpr_expr_build");
    check_sym(h, "hexpr_expr_destroy");
    check_sym(h, "hexpr_expr_nterms_full");
    check_sym(h, "hexpr_expr_nterms_one");
    check_sym(h, "hexpr_expr_read");
    check_sym(h, "hexpr_expr_show_full");
    check_sym(h, "hexpr_expr_show_one");

    /* F77 underscore variants */
    check_sym(h, "hexpr_init_");
    check_sym(h, "hexpr_shutdown_");
    check_sym(h, "hexpr_version_string_");
    check_sym(h, "hexpr_version_major_");
    check_sym(h, "hexpr_version_minor_");
    check_sym(h, "hexpr_version_patch_");
    check_sym(h, "hexpr_last_error_");
    check_sym(h, "hexpr_drt_create_");
    check_sym(h, "hexpr_drt_destroy_");
    check_sym(h, "hexpr_drt_norb_");
    check_sym(h, "hexpr_drt_nnodes_");
    check_sym(h, "hexpr_drt_ncsfs_");
    check_sym(h, "hexpr_drt_kind_");
    check_sym(h, "hexpr_drt_csfs_");
    check_sym(h, "hexpr_drt_node_info_");
    check_sym(h, "hexpr_drt_show_");
    check_sym(h, "hexpr_expr_build_");
    check_sym(h, "hexpr_expr_destroy_");
    check_sym(h, "hexpr_expr_nterms_full_");
    check_sym(h, "hexpr_expr_nterms_one_");
    check_sym(h, "hexpr_expr_read_");
    check_sym(h, "hexpr_expr_show_full_");
    check_sym(h, "hexpr_expr_show_one_");

    /* Binding-completeness additions (v1.0.0) */
    check_sym(h, "hexpr_drt_from_csfset_");
    check_sym(h, "hexpr_drt_ncore_");
    check_sym(h, "hexpr_drt_nval_");
    check_sym(h, "hexpr_drt_next_");
    check_sym(h, "hexpr_drt_nve_state_");
    check_sym(h, "hexpr_drt_spin_x2_");
    check_sym(h, "hexpr_csf_index_");
    check_sym(h, "hexpr_csf_steps_");
    check_sym(h, "hexpr_step_parse_");

    /* Canonical GUGA-walk view. Declared in hexpr_fortran.h since the
     * CSFSET work, but never listed here until now. */
    check_sym(h, "hexpr_drt_walk_ncsfs_");
    check_sym(h, "hexpr_drt_walk_csfs_");

    /* Per-pair / block / subspace evaluation. Likewise declared all
     * along and never listed: this file asserted that the DRT and
     * expression layers were reachable from Fortran while the entire
     * evaluation API was absent from the check. */
    check_sym(h, "hexpr_pair_count_");
    check_sym(h, "hexpr_pair_max_terms_");
    check_sym(h, "hexpr_pair_eval_");
    check_sym(h, "hexpr_block_count_");
    check_sym(h, "hexpr_block_max_terms_");
    check_sym(h, "hexpr_block_eval_");
    check_sym(h, "hexpr_subspace_count_");
    check_sym(h, "hexpr_subspace_max_terms_");
    check_sym(h, "hexpr_subspace_eval_");

    /* CSF -> determinant expansion */
    check_sym(h, "hexpr_csf_ndets_");
    check_sym(h, "hexpr_csf_to_dets_");

    dlclose(h);

    printf("%s (%d failure(s))\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}

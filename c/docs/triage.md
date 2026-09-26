# SPEC-DOX-1 — whole-tree function triage

Analysis only. No `.c`/`.h` file was edited. Classification is judgment over
line-count:

- **API**  — behavior clear from the signature alone: trivial wrappers,
  accessors, simple arithmetic, thin glue, `free()`/getters.
- **RICH** — carries non-obvious algorithmic intent: DRT construction, expr
  generation core, CSF lookup/pairing logic, non-trivial loops/recursion;
  "what it's doing and why" needs prose later.
- **DOMAIN? = YES** — intent cannot be fully explained from the code alone;
  needs Hiroaki / Noro-san GUGA / quantum-chemistry-physics input. Flagging a
  genuine gap is more useful than bluffing.

`~body_LOC` is approximate. Scope: `c/src/common/*.c`, `c/src/loopgen/*.c`,
`c/src/csf/*.c`, `c/src/hexpr_api.c`; `hexpr_fortran.c` + `run_reference.c`
listed as low-priority glue/driver (default API); `c/wig` excluded from prose
scope entirely (Sasaki-sensei's code — never inject there).

---

## common/ — generic data-structure utilities

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| common/darr.c | darr_init | 3 | API | Initializes array to empty state | NO |
| common/darr.c | darr_grow | 3 | API | Doubles capacity or initializes to 4 | NO |
| common/darr.c | darr_push | 2 | API | Appends item, grows if needed | NO |
| common/darr.c | darr_get | 1 | API | Returns item at index | NO |
| common/darr.c | darr_size | 1 | API | Returns current size | NO |
| common/darr.c | darr_clear | 1 | API | Resets size to zero | NO |
| common/darr.c | darr_free | 4 | API | Frees memory and resets state | NO |
| common/dhash.c | hash_key | 7 | RICH | FNV-1a hash over two int64 fields | NO |
| common/dhash.c | key_eq | 1 | API | Compares two keys for equality | NO |
| common/dhash.c | probe | 15 | RICH | Linear probing with deleted-slot tracking | NO |
| common/dhash.c | dhash_create | 4 | API | Allocates and initializes hash table | NO |
| common/dhash.c | dhash_free | 9 | API | Frees table, optional per-value callback | NO |
| common/dhash.c | dhash_rehash | 16 | RICH | Rebuilds table with doubled capacity | NO |
| common/dhash.c | dhash_get | 2 | API | Returns value for key or NULL | NO |
| common/dhash.c | dhash_has | 1 | API | Checks if key exists | NO |
| common/dhash.c | dhash_set | 13 | RICH | Insert/update, load-factor-driven rehashing | NO |
| common/dhash.c | dhash_delete | 5 | API | Marks entry as deleted | NO |
| common/dhash.c | dhash_size | 1 | API | Returns number of occupied entries | NO |
| common/dhash.c | dhash_to_array | 9 | API | Collects live entries into arrays | NO |
| common/fmt.c | hexpr_print_double | 1 | API | Prints double with 24.16e format | NO |
| common/fmt.c | hexpr_print_csf | 6 | API | Prints char array as JSON array | NO |
| common/sort.c | merge_run | 11 | RICH | Merges sorted runs using temp buffer | NO |
| common/sort.c | msort | 5 | RICH | Recursive merge-sort implementation | NO |
| common/sort.c | hexpr_stable_sort | 8 | API | Stable-sort wrapper with qsort fallback | NO |

**common subtotal:** 24 functions — RICH 6, API 18, DOMAIN 0.

---

## loopgen/drt.c — Distinct Row Table construction

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| drt.c | htab_hash | 1 | API | Simple hash computation on integers | NO |
| drt.c | htab_new | 6 | API | Allocate and initialize hash table | NO |
| drt.c | htab_free | 13 | API | Free hash table and chain nodes | NO |
| drt.c | htab_lookup | 9 | API | Standard hash table lookup | NO |
| drt.c | htab_insert | 16 | API | Hash insertion with array growth | NO |
| drt.c | mo_init | 13 | API | Initialize MO orbital classification struct | NO |
| drt.c | node_new | 20 | RICH | State-space validation logic (orb,nve,is) | YES |
| drt.c | node_if_connect_area | 11 | RICH | Check connectivity to specific GUGA targets | YES |
| drt.c | h_store | 6 | API | Wrapper around htab lookup/insert with dedupe | NO |
| drt.c | soci_generate_core_node | 15 | RICH | Generate core orbital nodes per SOCI rules | YES |
| drt.c | foci_generate_core_node | 12 | RICH | Generate core orbital nodes per FOCI rules | YES |
| drt.c | soci_generate_external_node | 18 | RICH | Generate external orbital nodes per SOCI | YES |
| drt.c | foci_generate_external_node | 12 | RICH | Generate external orbital nodes per FOCI | YES |
| drt.c | generate_valence_node | 33 | RICH | Recursive valence node generation with connectivity | YES |
| drt.c | prune_unreachable | 93 | RICH | Bidirectional BFS reachability pruning | NO |
| drt.c | node_cmp | 6 | API | Lexicographic node comparison for qsort | NO |
| drt.c | imax | 1 | API | Return maximum of two integers | NO |
| drt.c | store_number_upper_arcs | 21 | RICH | Recursive upper arc count + weight computation | YES |
| drt.c | store_number_lower_arcs | 17 | RICH | Recursive lower arc count computation | YES |
| drt.c | ivec_push | 6 | API | Dynamic int-vector append with growth | NO |
| drt.c | tree_search_lower | 12 | RICH | Recursive lower path weight enumeration | YES |
| drt.c | tree_search_upper | 12 | RICH | Recursive upper path weight enumeration | YES |
| drt.c | csf_push | 9 | API | Allocate and append CSF string to list | NO |
| drt.c | tree_search_csf | 14 | RICH | Recursive CSF walk enumeration with encoding | YES |
| drt.c | hexpr_drt_create | 180 | RICH | Multi-step DRT construction pipeline (SOCI/FOCI) | YES |
| drt.c | hexpr_drt_from_csfset | 210 | RICH | Build DRT from external CSF set with validation | YES |
| drt.c | hexpr_drt_free | 14 | API | Free all DRT-owned memory | NO |
| drt.c | hexpr_drt_destroy | 1 | API | Thin wrapper to hexpr_drt_free | NO |
| drt.c | hexpr_drt_norb | 1 | API | Return orbital count accessor | NO |
| drt.c | hexpr_drt_nnodes | 1 | API | Return node count accessor | NO |
| drt.c | hexpr_drt_ncsfs | 3 | API | Return CSF count with CSFSET logic | NO |
| drt.c | hexpr_drt_kind | 1 | API | Return DRT kind accessor | NO |
| drt.c | hexpr_drt_ncore | 1 | API | Return core orbital count | NO |
| drt.c | hexpr_drt_nval | 1 | API | Return valence orbital count | NO |
| drt.c | hexpr_drt_next | 1 | API | Return external orbital count | NO |
| drt.c | hexpr_drt_nve_state | 1 | API | Return NVE state accessor | NO |
| drt.c | hexpr_drt_spin_x2 | 1 | API | Return spin*2 state accessor | NO |
| drt.c | hexpr_drt_node_info | 7 | API | Get node (orb,nve,is) with bounds check | NO |
| drt.c | hexpr_drt_csfs | 3 | API | Return CSF buffer (user-facing or canonical) | NO |
| drt.c | hexpr_drt_walk_csfs | 2 | API | Return canonical CSF buffer | NO |
| drt.c | hexpr_drt_walk_ncsfs | 2 | API | Return canonical CSF count | NO |
| drt.c | hexpr_drt_orb | 1 | API | Return orbital of node i | NO |
| drt.c | hexpr_drt_nve | 1 | API | Return NVE of node i | NO |
| drt.c | hexpr_drt_spin | 1 | API | Return spin of node i with offset | NO |
| drt.c | hexpr_drt_upper_arc | 1 | API | Return upper arc accessor | NO |
| drt.c | hexpr_drt_lower_arc | 1 | API | Return lower arc accessor | NO |
| drt.c | hexpr_drt_arc_weight | 1 | API | Return arc weight accessor | NO |
| drt.c | hexpr_drt_number_upper_arcs | 1 | API | Return upper arc count accessor | NO |
| drt.c | hexpr_drt_number_lower_arcs | 1 | API | Return lower arc count accessor | NO |
| drt.c | hexpr_drt_lower_path_weights | 3 | API | Return lower path weights with count | NO |
| drt.c | hexpr_drt_upper_path_weights | 3 | API | Return upper path weights with count | NO |
| drt.c | hexpr_drt_mk_csf_list | 13 | RICH | Orchestrate CSF enumeration from DRT | YES |
| drt.c | hexpr_drt_show | 70 | RICH | Formatted tabular DRT output | NO |

**drt.c subtotal:** 53 functions — RICH 17, API 36, DOMAIN 15.

## loopgen/expr.c — coupling-coefficient expression generation

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| loopgen/expr.c | imax2 | 1 | API | trivial comparator | NO |
| loopgen/expr.c | imin2 | 1 | API | trivial comparator | NO |
| loopgen/expr.c | triangle | 5 | API | triangular number formula | NO |
| loopgen/expr.c | tri | 5 | API | triangular number formula | NO |
| loopgen/expr.c | ct_cmp_key | 5 | API | qsort comparator for keys | NO |
| loopgen/expr.c | ct_cmp_ij_only | 3 | API | qsort comparator by ij only | NO |
| loopgen/expr.c | free_ct | 1 | API | free wrapper for hash values | NO |
| loopgen/expr.c | add_expr | 85 | RICH | accumulates segment coefficients via hash table | YES |
| loopgen/expr.c | ot_cmp_ij | 3 | API | qsort comparator for ij | NO |
| loopgen/expr.c | expand_expr | 23 | RICH | expands intermediate terms to output form | YES |
| loopgen/expr.c | append_ct | 10 | API | vector append with doubling growth | NO |
| loopgen/expr.c | TP | 8 | API | tensor product struct initializer | NO |
| loopgen/expr.c | bc_init_constants | 4 | API | initializes sqrt(2) and sqrt(3) | NO |
| loopgen/expr.c | case01 | 10 | RICH | segment factor case 01 coefficient generation | YES |
| loopgen/expr.c | case02 | 10 | RICH | segment factor case 02 coefficient generation | YES |
| loopgen/expr.c | case1 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case2 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case3 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case4 | 10 | RICH | segment factor case 4 coefficient generation | YES |
| loopgen/expr.c | case5 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case6 | 10 | RICH | segment factor case 6 coefficient generation | YES |
| loopgen/expr.c | case7 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case8 | 10 | RICH | segment factor case 8 coefficient generation | YES |
| loopgen/expr.c | case10 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case11 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case12 | 10 | RICH | segment factor case 12 coefficient generation | YES |
| loopgen/expr.c | case13 | 18 | RICH | dual-result segment factor case | YES |
| loopgen/expr.c | case14 | 10 | RICH | segment factor case 14 coefficient generation | YES |
| loopgen/expr.c | brooks_all | 58 | RICH | orchestrates all Brooks segment cases + final sort | YES |
| loopgen/expr.c | hexpr_expr_build | 26 | RICH | builds full coupled expression from DRT | YES |
| loopgen/expr.c | hexpr_expr_free | 6 | API | deallocates expression struct memory | NO |
| loopgen/expr.c | hexpr_expr_destroy | 1 | API | destroy alias for free function | NO |
| loopgen/expr.c | hexpr_expr_nterms_full | 3 | API | returns total term count accessor | NO |
| loopgen/expr.c | hexpr_expr_nterms_one | 3 | API | returns one-electron term count | NO |
| loopgen/expr.c | hexpr_expr_read | 14 | API | copies filtered expression terms to output arrays | NO |
| loopgen/expr.c | write_header | 9 | API | prints header and CSF list | NO |
| loopgen/expr.c | hexpr_expr_show_full | 10 | API | prints full expression with all terms | NO |
| loopgen/expr.c | hexpr_expr_show_one | 12 | API | prints one-electron expression terms only | NO |

**expr.c subtotal:** 38 functions — RICH 19, API 19, DOMAIN 19.

## loopgen/expr_tree.c — tensor-operator DRT traversal

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| loopgen/expr_tree.c | matrix_el | 37 | RICH | Wigner matrix element for operator cases | YES |
| loopgen/expr_tree.c | hmlkva | 28 | RICH | HML angular-momentum coupling coefficient computation | YES |
| loopgen/expr_tree.c | rhash_init | 4 | API | record hash bucket initialization | NO |
| loopgen/expr_tree.c | fnv64 | 8 | API | FNV-1a 64-bit hash function | NO |
| loopgen/expr_tree.c | rhash_has_or_add | 23 | API | deduplication hash lookup and insert | NO |
| loopgen/expr_tree.c | rhash_reset | 13 | API | cleanup record hash with pool management | NO |
| loopgen/expr_tree.c | le_push | 7 | API | loop element vector push append | NO |
| loopgen/expr_tree.c | le_cmp_top | 3 | API | qsort comparator for top field | NO |
| loopgen/expr_tree.c | tree_search | 83 | RICH | recursive DRT graph walk applying segment rules | YES |
| loopgen/expr_tree.c | calc_TensorOp_average | 60 | RICH | main tensor product DRT traversal entry point | YES |

**expr_tree.c subtotal:** 10 functions — RICH 4, API 6, DOMAIN 4.

**loopgen total:** 101 functions — RICH 40, API 61, DOMAIN 38.

---

## csf/csf_lookup.c — CSF ↔ canonical-index mapping

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| csf/csf_lookup.c | hexpr_drt_walk_csf_index | 14 | RICH | DRT walk top-down, summing arc weights for index | YES |
| csf/csf_lookup.c | hexpr_drt_walk_csf_steps | 25 | RICH | Inverse DRT walk; recover step codes from index | YES |
| csf/csf_lookup.c | hexpr_csf_index | 20 | RICH | CSFSET linear search vs SOCI/FOCI canonical dispatch | YES |
| csf/csf_lookup.c | hexpr_csf_steps | 12 | RICH | CSFSET direct copy vs SOCI/FOCI canonical dispatch | YES |
| csf/csf_lookup.c | hexpr_step_parse | 14 | API | Parse eudf string to step code bytes | NO |

**csf_lookup.c subtotal:** 5 functions — RICH 4, API 1, DOMAIN 4.

## csf/csf_pair.c — pair / block / subspace orchestration

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| csf/csf_pair.c | imax2 | 1 | API | Inline max function | NO |
| csf/csf_pair.c | imin2 | 1 | API | Inline min function | NO |
| csf/csf_pair.c | triangle | 1 | RICH | Symmetric pair indexing formula | NO |
| csf/csf_pair.c | pair_via_subdrt | 33 | RICH | Create per-pair sub-DRT, compute coupling, patch ij | YES |
| csf/csf_pair.c | hexpr_pair_count | 16 | RICH | Count pair terms, filter one-electron interactions | YES |
| csf/csf_pair.c | hexpr_pair_eval | 21 | RICH | Evaluate pair, extract ij/pqrs/coef with filter | YES |
| csf/csf_pair.c | block_eval_impl | 42 | RICH | Nested bra/ket loop, dual count-only or emit mode | YES |
| csf/csf_pair.c | hexpr_block_count | 7 | API | Count-only wrapper to block_eval_impl | NO |
| csf/csf_pair.c | hexpr_block_eval | 10 | API | Full eval wrapper to block_eval_impl | NO |
| csf/csf_pair.c | hexpr_subspace_count | 4 | API | Subspace count via block with same indices | NO |
| csf/csf_pair.c | hexpr_subspace_eval | 7 | API | Subspace eval via block with same indices | NO |
| csf/csf_pair.c | per_pair_max_terms | 4 | API | Max terms formula per pair type | NO |
| csf/csf_pair.c | hexpr_pair_max_terms | 5 | API | Pair max terms wrapper | NO |
| csf/csf_pair.c | hexpr_block_max_terms | 10 | API | Block max terms with int64 overflow check | NO |
| csf/csf_pair.c | hexpr_subspace_max_terms | 10 | API | Subspace max terms with overflow check | NO |

**csf_pair.c subtotal:** 15 functions — RICH 5, API 10, DOMAIN 4.

## csf/csf_pair_core.c — per-pair 14-case coupling core

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| csf/csf_pair_core.c | init_constants | 1 | API | Initialize sqrt(2) and sqrt(3) constants | NO |
| csf/csf_pair_core.c | imax2 | 1 | API | Inline max function | NO |
| csf/csf_pair_core.c | imin2 | 1 | API | Inline min function | NO |
| csf/csf_pair_core.c | triangle | 1 | RICH | Symmetric pair indexing formula | NO |
| csf/csf_pair_core.c | tri | 1 | RICH | Symmetric orbital pair indexing | NO |
| csf/csf_pair_core.c | TP | 8 | API | Construct tensor_product_t from operator/theta arrays | NO |
| csf/csf_pair_core.c | free_pt | 1 | API | Free function wrapper | NO |
| csf/csf_pair_core.c | add_expr_pair | 84 | RICH | D3 filter, dedup by ij_key/seq_num, expand paths | YES |
| csf/csf_pair_core.c | append_ot | 3 | API | Append output term vector | NO |
| csf/csf_pair_core.c | case01_pair | 8 | RICH | CAS coupling case; tensor product, filter, coeff | YES |
| csf/csf_pair_core.c | case02_pair | 8 | RICH | CA coupling case; tensor product, filter, coeff | YES |
| csf/csf_pair_core.c | case03_pair | 8 | RICH | AC coupling case; tensor product, filter, coeff | YES |
| csf/csf_pair_core.c | case1_pair | 16 | RICH | ACCA coupling; dual tensor products, index reorder | YES |
| csf/csf_pair_core.c | case2_pair | 16 | RICH | CACA coupling; dual tensor products, index reorder | YES |
| csf/csf_pair_core.c | case3_pair | 16 | RICH | CCAA coupling; dual tensor products, index reorder | YES |
| csf/csf_pair_core.c | case4_pair | 8 | RICH | A-CCS-A coupling case | YES |
| csf/csf_pair_core.c | case5_pair | 16 | RICH | C-CAS/CAT-A coupling; dual tensor products | YES |
| csf/csf_pair_core.c | case6_pair | 8 | RICH | CCS-A-A coupling case | YES |
| csf/csf_pair_core.c | case7_pair | 16 | RICH | CAS/CAT-C-A coupling; dual tensor products | YES |
| csf/csf_pair_core.c | case8_pair | 8 | RICH | C-C-AAS coupling case | YES |
| csf/csf_pair_core.c | case10_pair | 16 | RICH | C-A-CAS/CAT coupling; dual tensor products | YES |
| csf/csf_pair_core.c | case11_pair | 16 | RICH | C-CAA / CCA-A coupling; dual tensor products | YES |
| csf/csf_pair_core.c | case12_pair | 8 | RICH | CCS-AAS coupling case | YES |
| csf/csf_pair_core.c | case13_pair | 16 | RICH | CAS-CAS / CAT-CAT coupling; dual tensor products | YES |
| csf/csf_pair_core.c | case14_pair | 8 | RICH | CCAA four-body coupling tensor product | YES |
| csf/csf_pair_core.c | pair_eval_internal | 22 | RICH | Orchestrate all 14 coupling cases, append results | YES |

**csf_pair_core.c subtotal:** 26 functions — RICH 20, API 6, DOMAIN 18.

**csf total:** 46 functions — RICH 29, API 17, DOMAIN 26.

---

## hexpr_api.c — public library lifecycle / version / error

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| hexpr_api.c | hexpr_version_string | 1 | API | Version string accessor | NO |
| hexpr_api.c | hexpr_version_major | 1 | API | Version major accessor | NO |
| hexpr_api.c | hexpr_version_minor | 1 | API | Version minor accessor | NO |
| hexpr_api.c | hexpr_version_patch | 1 | API | Version patch accessor | NO |
| hexpr_api.c | hexpr_last_error | 1 | API | Thread-local error buffer accessor | NO |
| hexpr_api.c | hexpr_set_error | 4 | API | Format error message via vsnprintf | NO |
| hexpr_api.c | hexpr_init | 4 | API | Library initialization, calls wigner_init | NO |
| hexpr_api.c | hexpr_shutdown | 1 | API | Lifecycle shutdown, reset flag | NO |

**hexpr_api.c subtotal:** 8 functions — RICH 0, API 8, DOMAIN 0.

---

## hexpr_fortran.c — Fortran-binding glue (low priority, default API)

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| hexpr_fortran.c | hexpr_fortran_stdout | 1 | API | Return stdout pointer for Fortran | NO |
| hexpr_fortran.c | hexpr_version_string_ | 1 | API | Fortran wrapper for version string | NO |
| hexpr_fortran.c | hexpr_version_major_ | 1 | API | Fortran wrapper for version major | NO |
| hexpr_fortran.c | hexpr_version_minor_ | 1 | API | Fortran wrapper for version minor | NO |
| hexpr_fortran.c | hexpr_version_patch_ | 1 | API | Fortran wrapper for version patch | NO |
| hexpr_fortran.c | hexpr_last_error_ | 1 | API | Fortran wrapper for last error | NO |
| hexpr_fortran.c | hexpr_init_ | 2 | API | Fortran wrapper for init | NO |
| hexpr_fortran.c | hexpr_shutdown_ | 1 | API | Fortran wrapper for shutdown | NO |
| hexpr_fortran.c | hexpr_drt_create_ | 4 | API | Fortran wrapper, arg casting for DRT creation | NO |
| hexpr_fortran.c | hexpr_drt_from_csfset_ | 2 | API | Fortran wrapper for DRT from CSF set | NO |
| hexpr_fortran.c | hexpr_drt_destroy_ | 1 | API | Fortran wrapper for DRT destruction | NO |
| hexpr_fortran.c | hexpr_drt_norb_ | 1 | API | Fortran wrapper for DRT norb accessor | NO |
| hexpr_fortran.c | hexpr_drt_nnodes_ | 1 | API | Fortran wrapper for DRT nnodes accessor | NO |
| hexpr_fortran.c | hexpr_drt_ncsfs_ | 1 | API | Fortran wrapper for DRT ncsfs accessor | NO |
| hexpr_fortran.c | hexpr_drt_walk_ncsfs_ | 1 | API | Fortran wrapper for walk ncsfs accessor | NO |
| hexpr_fortran.c | hexpr_drt_kind_ | 1 | API | Fortran wrapper for DRT kind accessor | NO |
| hexpr_fortran.c | hexpr_drt_ncore_ | 1 | API | Fortran wrapper for ncore accessor | NO |
| hexpr_fortran.c | hexpr_drt_nval_ | 1 | API | Fortran wrapper for nval accessor | NO |
| hexpr_fortran.c | hexpr_drt_next_ | 1 | API | Fortran wrapper for next accessor | NO |
| hexpr_fortran.c | hexpr_drt_nve_state_ | 1 | API | Fortran wrapper for nve state accessor | NO |
| hexpr_fortran.c | hexpr_drt_spin_x2_ | 1 | API | Fortran wrapper for spin x2 accessor | NO |
| hexpr_fortran.c | hexpr_drt_csfs_ | 6 | API | Copy CSF buffer for Fortran caller | NO |
| hexpr_fortran.c | hexpr_drt_walk_csfs_ | 6 | API | Copy canonical GUGA walk CSF buffer | NO |
| hexpr_fortran.c | hexpr_drt_node_info_ | 4 | API | Fortran wrapper for DRT node info lookup | NO |
| hexpr_fortran.c | hexpr_drt_show_ | 1 | API | Fortran wrapper for DRT show | NO |
| hexpr_fortran.c | hexpr_expr_build_ | 1 | API | Fortran wrapper for expression builder | NO |
| hexpr_fortran.c | hexpr_expr_destroy_ | 1 | API | Fortran wrapper for expression destroy | NO |
| hexpr_fortran.c | hexpr_expr_nterms_full_ | 1 | API | Fortran wrapper for nterms full accessor | NO |
| hexpr_fortran.c | hexpr_expr_nterms_one_ | 1 | API | Fortran wrapper for nterms one accessor | NO |
| hexpr_fortran.c | hexpr_expr_read_ | 3 | API | Fortran wrapper for expression term read | NO |
| hexpr_fortran.c | hexpr_expr_show_full_ | 1 | API | Fortran wrapper for show full | NO |
| hexpr_fortran.c | hexpr_expr_show_one_ | 1 | API | Fortran wrapper for show one | NO |
| hexpr_fortran.c | hexpr_csf_index_ | 3 | API | Fortran wrapper for CSF index lookup | NO |
| hexpr_fortran.c | hexpr_csf_steps_ | 3 | API | Fortran wrapper for CSF steps | NO |
| hexpr_fortran.c | hexpr_step_parse_ | 10 | API | Parse Fortran string to step codes | NO |
| hexpr_fortran.c | hexpr_pair_count_ | 2 | API | Fortran wrapper for pair term count | NO |
| hexpr_fortran.c | hexpr_pair_max_terms_ | 2 | API | Fortran wrapper for pair max terms | NO |
| hexpr_fortran.c | hexpr_pair_eval_ | 3 | API | Fortran wrapper for pair evaluation | NO |
| hexpr_fortran.c | hexpr_block_count_ | 3 | API | Fortran wrapper for block count | NO |
| hexpr_fortran.c | hexpr_block_max_terms_ | 3 | API | Fortran wrapper for block max terms | NO |
| hexpr_fortran.c | hexpr_block_eval_ | 6 | API | Fortran wrapper for block evaluation | NO |
| hexpr_fortran.c | hexpr_subspace_count_ | 2 | API | Fortran wrapper for subspace count | NO |
| hexpr_fortran.c | hexpr_subspace_max_terms_ | 2 | API | Fortran wrapper for subspace max terms | NO |
| hexpr_fortran.c | hexpr_subspace_eval_ | 6 | API | Fortran wrapper for subspace evaluation | NO |

**hexpr_fortran.c subtotal:** 44 functions — RICH 0, API 44, DOMAIN 0.

---

## run_reference.c — reference CLI driver (low priority, default API)

| file | function | ~body_LOC | class | reason | DOMAIN? |
|------|----------|-----------|-------|--------|---------|
| run_reference.c | path_size | 3 | API | Get file size via stat syscall | NO |
| run_reference.c | wall_elapsed | 5 | API | Compute elapsed wall time, CLOCK_MONOTONIC | NO |
| run_reference.c | generate_case | 63 | API | Generate DRT, CSF list, full/one-electron expr | NO |
| run_reference.c | generate_pair_outputs | 124 | API | Generate per-pair expr output, dynamic buffers | NO |
| run_reference.c | main | 86 | API | CLI parser, case filter, output dispatch | NO |

**run_reference.c subtotal:** 5 functions — RICH 0, API 5, DOMAIN 0.

---

## Counts

| group | total | RICH | API | DOMAIN? |
|-------|-------|------|-----|---------|
| common/darr.c | 7 | 0 | 7 | 0 |
| common/dhash.c | 12 | 4 | 8 | 0 |
| common/fmt.c | 2 | 0 | 2 | 0 |
| common/sort.c | 3 | 2 | 1 | 0 |
| loopgen/drt.c | 53 | 17 | 36 | 15 |
| loopgen/expr.c | 38 | 19 | 19 | 19 |
| loopgen/expr_tree.c | 10 | 4 | 6 | 4 |
| csf/csf_lookup.c | 5 | 4 | 1 | 4 |
| csf/csf_pair.c | 15 | 5 | 10 | 4 |
| csf/csf_pair_core.c | 26 | 20 | 6 | 18 |
| hexpr_api.c | 8 | 0 | 8 | 0 |
| hexpr_fortran.c | 44 | 0 | 44 | 0 |
| run_reference.c | 5 | 0 | 5 | 0 |
| **TOTAL** | **228** | **75** | **153** | **64** |

### Where the prose work concentrates (RICH)

- **loopgen/** — 40 RICH (DRT construction in `drt.c`; the `case01..case14`
  segment-factor rules and `add_expr`/`brooks_all` generation core in `expr.c`;
  `tree_search`/`calc_TensorOp_average`/`matrix_el`/`hmlkva` in `expr_tree.c`).
- **csf/** — 29 RICH (canonical CSF↔index walks in `csf_lookup.c`;
  `pair_via_subdrt`/`block_eval_impl` in `csf_pair.c`; the per-pair
  `case*_pair` 14-case core + `add_expr_pair` in `csf_pair_core.c`).
- **common/** — 6 RICH (hash internals `hash_key`/`probe`/`rehash`/`set`;
  merge-sort `merge_run`/`msort`). Generic, no domain input needed.
- **hexpr_api.c / hexpr_fortran.c / run_reference.c** — 0 RICH; all
  accessors, Fortran wrappers, and CLI plumbing → API.

### DOMAIN-flagged (64) — needs Hiroaki / Noro-san GUGA / physics input

All 64 are in `loopgen/` (38) and `csf/` (26). They cluster on:

1. **Segment-factor / coupling-coefficient rules** — `case01..case14`
   (`expr.c`) and `case*_pair` (`csf_pair_core.c`): why each step-vector case
   yields its specific coefficients (Brooks–Moshinsky / Shavitt conventions).
2. **DRT node generation & connectivity** — `node_new`,
   `node_if_connect_area`, `{soci,foci}_generate_*_node`,
   `generate_valence_node`, `hexpr_drt_create`, `hexpr_drt_from_csfset`:
   SOCI/FOCI state-space and step-vector validity rules.
3. **Arc-weight / walk indexing** — `store_number_{upper,lower}_arcs`,
   `tree_search_*`, `hexpr_drt_walk_csf_*`, `hexpr_csf_index/steps`: the GUGA
   lexical-index (Paldus/Shavitt) weighting convention.
4. **Tensor-operator matrix elements** — `matrix_el`, `hmlkva`, `tree_search`,
   `calc_TensorOp_average`: angular-momentum coupling / Wigner factors.

Reachability pruning (`prune_unreachable`), formatted output (`hexpr_drt_show`),
and the symmetric index helpers (`triangle`/`tri`) are RICH but **not**
DOMAIN-flagged — their intent is explainable from the code alone.

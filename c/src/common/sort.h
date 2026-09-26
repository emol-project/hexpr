/* common/sort.h -- stable sort for hexpr
 *
 * DECISION: use hexpr_stable_sort (merge sort) everywhere output order must
 * match Ruby Array#sort, which has been stable (merge sort) since Ruby 1.9.
 * C's qsort is NOT guaranteed stable -- do not use it for expression-term or
 * node arrays where equal primary keys must preserve insertion order.
 *
 * RULE: a comparator may return 0 only when two elements are truly identical
 * for all purposes, OR callers must use hexpr_stable_sort (not qsort).
 *
 * API mirrors qsort(3): same signature, stable semantics.
 * Falls back to qsort if malloc fails (should not happen in practice).
 */
#ifndef HEXPR_SORT_H_
#define HEXPR_SORT_H_

#include <stddef.h>

void hexpr_stable_sort(void *base, size_t nmemb, size_t size,
                       int (*compar)(const void *, const void *));

#endif /* HEXPR_SORT_H_ */

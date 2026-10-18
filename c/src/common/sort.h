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
 * API mirrors qsort(3), with stable semantics, but returns 0 on success and
 * -1 if the scratch buffer cannot be allocated; the array is then left
 * unsorted (since 1.2: an unstable qsort fallback would silently change the
 * output order, so the caller reports HEXPR_ERR_OOM instead).
 */
#ifndef HEXPR_SORT_H_
#define HEXPR_SORT_H_

#include <stddef.h>

int hexpr_stable_sort(void *base, size_t nmemb, size_t size,
                      int (*compar)(const void *, const void *));

#endif /* HEXPR_SORT_H_ */

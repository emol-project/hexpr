/* common/dhash.h -- open-addressing hash table keyed on (int64, int64).
 *
 * Used by BrooksCases::add_expr to accumulate expression terms keyed on
 * (ij, pqrs).  Output is always sorted before use, so iteration order
 * (bucket order) does not need to match Ruby's insertion order.
 */
#ifndef HEXPR_DHASH_H_
#define HEXPR_DHASH_H_

#include <stdint.h>

typedef struct {
    int64_t k0;   /* first key component  (triangle index / ij) */
    int64_t k1;   /* second key component (orbital index / pqrs) */
} dhash_key_t;

typedef struct dhash dhash_t;

dhash_t *dhash_create(void);

/* Free the table.  If free_val != NULL, call it on every stored value. */
void     dhash_free  (dhash_t *h, void (*free_val)(void *));

/* Returns stored value, or NULL if not present. */
void    *dhash_get   (const dhash_t *h, dhash_key_t key);

/* Returns 1 if key is present, 0 otherwise. */
int      dhash_has   (const dhash_t *h, dhash_key_t key);

/* Insert or update.  val must not be NULL (NULL means absent). */
void     dhash_set   (dhash_t *h, dhash_key_t key, void *val);

/* Remove key; no-op if not present. */
void     dhash_delete(dhash_t *h, dhash_key_t key);

/* Number of live entries. */
int      dhash_size  (const dhash_t *h);

/* Fill caller-allocated arrays (each at least dhash_size() elements).
 * Either output pointer may be NULL.  Returns number of entries written. */
int      dhash_to_array(const dhash_t *h,
                        dhash_key_t *keys_out, void **vals_out);

#endif /* HEXPR_DHASH_H_ */

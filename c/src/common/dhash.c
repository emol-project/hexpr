#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "dhash.h"

#define DHASH_EMPTY    0
#define DHASH_OCCUPIED 1
#define DHASH_DELETED  2

typedef struct {
    dhash_key_t key;
    void       *val;
    int         state;
} dhash_bucket_t;

struct dhash {
    dhash_bucket_t *buckets;
    int             n_buckets;
    int             n_occupied; /* live entries */
    int             n_total;    /* occupied + deleted (load-factor denominator) */
};

/* FNV-1a over the two int64 fields */
/**
 * @brief Compute the FNV-1a hash of a key over its 16 raw bytes (k0 then k1).
 * @param k key to hash (passed by value).
 * @return the 64-bit FNV-1a digest.
 * @note Hashes the in-memory byte representation, so the digest is
 *       byte-order (endianness) dependent; this is harmless because the table
 *       is never serialized -- the hash is only ever used within one process.
 */
static uint64_t hash_key(dhash_key_t k)
{
    uint64_t h = 14695981039346656037ULL;
    const unsigned char *p = (const unsigned char *)&k.k0;
    int i;
    for (i = 0; i < 8; i++) { h ^= *p++; h *= 1099511628211ULL; }
    p = (const unsigned char *)&k.k1;
    for (i = 0; i < 8; i++) { h ^= *p++; h *= 1099511628211ULL; }
    return h;
}

/**
 * @brief Test two keys for component-wise equality.
 * @param a,b keys to compare (by value).
 * @return non-zero if @c k0 and @c k1 both match, 0 otherwise.
 */
static int key_eq(dhash_key_t a, dhash_key_t b)
{
    return a.k0 == b.k0 && a.k1 == b.k1;
}

/* Probe for key.  insert=0: return NULL if absent.
 *                 insert=1: return slot to write into. */
/**
 * @brief Linear-probe the open-addressing table for @p key.
 * @param h      table to probe.
 * @param key    key being located.
 * @param insert 0 = lookup (return NULL if absent); 1 = locate a slot to write.
 * @return for lookup: the occupied bucket holding @p key, else NULL. For insert:
 *         the bucket to write (the key's own bucket if present, otherwise a
 *         reusable tombstone or the first empty slot), or NULL only if the table
 *         is full with no empty/tombstone slot.
 * @note Internal helper underlying @ref dhash_get / @ref dhash_set /
 *       @ref dhash_delete. Worst case O(n_buckets).
 */
static dhash_bucket_t *probe(dhash_t *h, dhash_key_t key, int insert)
{
    int nb = h->n_buckets;
    int idx = (int)(hash_key(key) % (uint64_t)nb);
    dhash_bucket_t *first_del = NULL;
    int i;
    for (i = 0; i < nb; i++) {
        dhash_bucket_t *b = &h->buckets[(idx + i) % nb];
        /* A DELETED slot is a tombstone: lookups must probe PAST it (the key may
         * sit further along), but an insert prefers to reuse the first one seen
         * rather than consume a fresh empty slot -- hence remember first_del and
         * only fall back to the empty bucket `b` if none was found. */
        if (b->state == DHASH_EMPTY)
            return insert ? (first_del ? first_del : b) : NULL;
        if (b->state == DHASH_DELETED) {
            if (!first_del) first_del = b;
            continue;
        }
        if (key_eq(b->key, key)) return b;
    }
    /* Full sweep with no empty slot: insert can still land on a tombstone. */
    return insert ? first_del : NULL;
}

/**
 * @brief Allocate an empty hash table (initial capacity 16 buckets).
 * @return a new table; release with @ref dhash_free. Pairs with @ref dhash_free.
 */
dhash_t *dhash_create(void)
{
    dhash_t *h = calloc(1, sizeof(dhash_t));
    h->n_buckets = 16;
    h->buckets   = calloc((size_t)h->n_buckets, sizeof(dhash_bucket_t));
    return h;
}

/**
 * @brief Free the table and its bucket array.
 * @param h        table to free; NULL is a no-op.
 * @param free_val optional value-destructor; if non-NULL it is invoked once on
 *                 every live (occupied) value before the table is released.
 * @note Stored values are owned by the caller unless @p free_val is supplied to
 *       reclaim them here. Pairs with @ref dhash_create.
 */
void dhash_free(dhash_t *h, void (*free_val)(void *))
{
    if (!h) return;
    if (free_val) {
        int i;
        for (i = 0; i < h->n_buckets; i++)
            if (h->buckets[i].state == DHASH_OCCUPIED)
                free_val(h->buckets[i].val);
    }
    free(h->buckets);
    free(h);
}

/**
 * @brief Double the bucket count and reinsert every live entry.
 * @param h table to grow.
 * @note Internal helper for @ref dhash_set. Reinserting only OCCUPIED entries
 *       drops all tombstones, so @c n_total is rebuilt to equal @c n_occupied
 *       (this is the only place the deleted-slot count is reclaimed). Values are
 *       relocated by pointer; the pointees are not touched.
 */
static void dhash_rehash(dhash_t *h)
{
    int old_nb = h->n_buckets;
    dhash_bucket_t *old_b = h->buckets;
    int i;
    h->n_buckets = old_nb * 2;
    h->buckets   = calloc((size_t)h->n_buckets, sizeof(dhash_bucket_t));
    /* Reset both counters: the loop below recounts from live entries only, so
     * tombstones vanish and n_total drops back to n_occupied. */
    h->n_occupied = 0;
    h->n_total    = 0;
    for (i = 0; i < old_nb; i++) {
        if (old_b[i].state == DHASH_OCCUPIED) {
            dhash_bucket_t *b = probe(h, old_b[i].key, 1);
            *b = old_b[i];
            h->n_occupied++;
            h->n_total++;
        }
    }
    free(old_b);
}

/**
 * @brief Look up the value stored under @p key.
 * @param h   table to search.
 * @param key key to find.
 * @return the stored value pointer, or NULL if @p key is absent.
 * @note Returns a borrowed pointer owned by the caller/table, not a copy.
 *       Because NULL doubles as "absent", values stored via @ref dhash_set are
 *       required to be non-NULL. O(1) average.
 */
void *dhash_get(const dhash_t *h, dhash_key_t key)
{
    dhash_bucket_t *b = probe((dhash_t *)h, key, 0);
    return (b && b->state == DHASH_OCCUPIED) ? b->val : NULL;
}

/**
 * @brief Test whether @p key is present.
 * @param h   table to query.
 * @param key key to test.
 * @return 1 if present, 0 otherwise.
 * @note Thin predicate over @ref dhash_get.
 */
int dhash_has(const dhash_t *h, dhash_key_t key)
{
    return dhash_get(h, key) != NULL;
}

/**
 * @brief Insert a new entry or overwrite the value of an existing key.
 * @param h   table to modify.
 * @param key key to insert/update.
 * @param val value to store; MUST be non-NULL (NULL is reserved to mean absent
 *            in @ref dhash_get). Stored by reference -- ownership stays with the
 *            caller; the pointee is not copied.
 * @note May trigger @ref dhash_rehash, which moves buckets and so invalidates
 *       any bucket pointer or iteration in progress (values/keys themselves
 *       survive). O(1) amortized.
 */
void dhash_set(dhash_t *h, dhash_key_t key, void *val)
{
    /* Load factor uses n_total (occupied + tombstones), NOT n_occupied, so a
     * table churned by many deletes still rehashes and reclaims its tombstones. */
    /* Resize when total (occupied+deleted) >= 70% capacity */
    if ((h->n_total + 1) * 10 > h->n_buckets * 7)
        dhash_rehash(h);
    dhash_bucket_t *b = probe(h, key, 1);
    if (b->state == DHASH_OCCUPIED) {
        b->val = val;
    } else {
        /* Bump n_total only for a fresh EMPTY slot; reusing a tombstone keeps
         * n_total unchanged because that slot was already counted in it. */
        if (b->state == DHASH_EMPTY) h->n_total++;
        b->key   = key;
        b->val   = val;
        b->state = DHASH_OCCUPIED;
        h->n_occupied++;
    }
}

/**
 * @brief Remove @p key from the table.
 * @param h   table to modify.
 * @param key key to remove; a no-op if absent.
 * @note Marks the slot as a tombstone (DHASH_DELETED) rather than clearing it,
 *       preserving probe chains; the slot is only reclaimed by a later
 *       @ref dhash_rehash. The stored value is NOT freed -- the caller still
 *       owns it.
 */
void dhash_delete(dhash_t *h, dhash_key_t key)
{
    dhash_bucket_t *b = probe(h, key, 0);
    if (b && b->state == DHASH_OCCUPIED) {
        b->state = DHASH_DELETED;
        h->n_occupied--;
    }
}

/**
 * @brief Number of live entries.
 * @param h table to query. @return live (non-tombstone) entry count.
 * @note Excludes tombstones; this is the count to size the arrays for
 *       @ref dhash_to_array.
 */
int dhash_size(const dhash_t *h)
{
    return h->n_occupied;
}

/**
 * @brief Copy all live keys and/or values into caller-provided arrays.
 * @param h        table to enumerate.
 * @param keys_out destination for keys, or NULL to skip keys.
 * @param vals_out destination for values, or NULL to skip values.
 * @return number of entries written (equals @ref dhash_size).
 * @note Each non-NULL output array must hold at least @ref dhash_size elements.
 *       Keys/values are written in bucket order, which is unspecified (callers
 *       that need a defined order must sort the result).
 */
int dhash_to_array(const dhash_t *h, dhash_key_t *keys_out, void **vals_out)
{
    int count = 0, i;
    for (i = 0; i < h->n_buckets; i++) {
        if (h->buckets[i].state == DHASH_OCCUPIED) {
            if (keys_out) keys_out[count] = h->buckets[i].key;
            if (vals_out) vals_out[count] = h->buckets[i].val;
            count++;
        }
    }
    return count;
}

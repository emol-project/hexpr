/* loopgen/drt.c -- DRT construction, shared SOCI/FOCI logic */
/**
 * @file drt.c
 * @brief Construction and accessors for the Distinct Row Table (DRT) -- the
 *        upstream root of the GUGA pipeline.
 *
 * A DRT is a layered directed graph whose nodes carry the triple
 * @c (orb, nve, is) and whose four arc slots per node correspond to the four
 * GUGA step types (indices 0..3, rendered as the characters @c 'e','u','d','f').
 * Everything downstream -- @c expr.c's coupling-coefficient generator and the
 * @c csf/ lookup code -- walks this graph and reads its arc weights and
 * path-weight lists.
 *
 * Two public builders produce a ::hexpr_drt_t:
 *   - @ref hexpr_drt_create   : parametric SOCI/FOCI construction from MO counts
 *     and a target state, in phases (core nodes, external nodes, valence nodes,
 *     prune, sort, wire arcs, compute weights, pack CSFs).
 *   - @ref hexpr_drt_from_csfset : build a sub-DRT from an explicit set of input
 *     CSFs, validating each against the Paldus non-negativity constraints.
 * Release with @ref hexpr_drt_destroy: that is the name declared in
 * @c hexpr.h. @ref hexpr_drt_free is the internal implementation it forwards
 * to; in-tree callers may use either, but public documentation and the
 * language bindings name only @c hexpr_drt_destroy.
 *
 * The accessors near the bottom (@ref hexpr_drt_norb, @ref hexpr_drt_arc_weight,
 * @ref hexpr_drt_lower_path_weights, ...) are the cross-link hubs referenced
 * from @c expr.c; they borrow DRT-internal storage and never transfer ownership.
 *
 * @note In-body comments describe only the data flow visible in the code; the
 *       GUGA / Paldus-Shavitt intent behind node admissibility, the arc-weight
 *       convention, and the step ordering is deferred to <tt>\@todo domain:</tt>
 *       markers rather than asserted here.
 */
#include "drt.h"
#include "hexpr_error.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */
/* Internal node representation                                        */
/* ------------------------------------------------------------------ */

#define ARC_NIL (-1)

typedef struct {
    int orb, nve, is;          /* code = [orb, nve, is] */
    int status;                /* 1=valid, 0=invalid */
    int upper_arc[4];          /* ARC_NIL or node index */
    int arc_weight[4];         /* ARC_NIL or weight */
    int number_upper_arcs;     /* ARC_NIL until set */
    int lower_arc[4];          /* ARC_NIL or node index */
    int number_lower_arcs;     /* ARC_NIL until set */
} node_t;

/* ------------------------------------------------------------------ */
/* Hash table for node lookup during construction                      */
/* Maps (orb,nve,is) -> node_t*                                       */
/* ------------------------------------------------------------------ */

#define HTAB_INIT 4096

typedef struct htab_entry {
    int orb, nve, is;
    node_t *node;
    struct htab_entry *next;
} htab_entry_t;

typedef struct {
    htab_entry_t **buckets;
    int nbuckets;
    int size;
    /* keep insertion order for iteration */
    node_t **order;
    int order_cap;
    int order_size;
} htab_t;

/**
 * @brief Mix a node code @c (orb,nve,is) into an unsigned bucket hash.
 * @param orb,nve,is the three node-code components.
 * @return the (unreduced) hash; callers take it modulo @c nbuckets.
 */
static unsigned int htab_hash(int orb, int nve, int is)
{
    unsigned int h = (unsigned int)(orb * 100003 + nve * 1009 + is);
    return h;
}

/**
 * @brief Allocate an empty construction hash table (code -> node, with an
 *        insertion-order list for stable iteration).
 * @return a new table; release with @ref htab_free.
 * @note Internal to construction only (not the public ::hexpr_drt_t).
 */
static htab_t *htab_new(void)
{
    htab_t *h = calloc(1, sizeof(htab_t));
    h->nbuckets = HTAB_INIT;
    h->buckets = calloc(h->nbuckets, sizeof(htab_entry_t *));
    h->order_cap = 256;
    h->order = malloc(h->order_cap * sizeof(node_t *));
    return h;
}

/**
 * @brief Free a construction table, its chains, AND every node it owns.
 * @param h table to free.
 * @note Calls @c free on each stored @c node_t (the table owns inserted nodes;
 *       see @ref htab_insert). Used to release the transient tables and, after
 *       the sorted node array is copied out, the nodes themselves.
 */
static void htab_free(htab_t *h)
{
    for (int i = 0; i < h->nbuckets; i++) {
        htab_entry_t *e = h->buckets[i];
        while (e) {
            htab_entry_t *next = e->next;
            free(e->node);
            free(e);
            e = next;
        }
    }
    free(h->buckets);
    free(h->order);
    free(h);
}

/**
 * @brief Find the node with code @c (orb,nve,is).
 * @param h        table to search.
 * @param orb,nve,is node code to look up.
 * @return the stored node, or NULL if absent (borrowed; owned by @p h).
 */
static node_t *htab_lookup(htab_t *h, int orb, int nve, int is)
{
    unsigned int idx = htab_hash(orb, nve, is) % (unsigned int)h->nbuckets;
    htab_entry_t *e = h->buckets[idx];
    while (e) {
        if (e->orb == orb && e->nve == nve && e->is == is)
            return e->node;
        e = e->next;
    }
    return NULL;
}

/* Insert a new node; does NOT check for duplicates.
 * Takes ownership of node. */
/**
 * @brief Append @p node to the table and its insertion-order list.
 * @param h    table to insert into.
 * @param node node to store; the table TAKES OWNERSHIP (freed by @ref htab_free).
 * @note Does not check for an existing code -- the caller must dedupe (see
 *       @ref h_store). Grows the order[] array by doubling when full.
 */
static void htab_insert(htab_t *h, node_t *node)
{
    if (h->order_size == h->order_cap) {
        h->order_cap *= 2;
        h->order = realloc(h->order, h->order_cap * sizeof(node_t *));
    }
    h->order[h->order_size++] = node;
    h->size++;

    unsigned int idx = htab_hash(node->orb, node->nve, node->is) % (unsigned int)h->nbuckets;
    htab_entry_t *e = malloc(sizeof(htab_entry_t));
    e->orb = node->orb;
    e->nve = node->nve;
    e->is  = node->is;
    e->node = node;
    e->next = h->buckets[idx];
    h->buckets[idx] = e;
}

/* ------------------------------------------------------------------ */
/* MO classification                                                   */
/* ------------------------------------------------------------------ */

typedef struct {
    int norb, ncore, nval, next;
    int core_first, core_last;   /* 1-based ranges (inclusive) */
    int val_first,  val_last;
    int ext_first,  ext_last;
} mo_t;

/**
 * @brief Build the MO classification: counts plus 1-based inclusive orbital
 *        ranges for the core, valence, and external blocks.
 * @param ncore     number of core orbitals.
 * @param nval      number of valence orbitals.
 * @param nexternal number of external orbitals.
 * @return a populated ::mo_t with @c norb = ncore+nval+nexternal and the three
 *         contiguous ranges laid out core, then valence, then external.
 */
static mo_t mo_init(int ncore, int nval, int nexternal)
{
    mo_t m;
    m.ncore = ncore;
    m.nval  = nval;
    m.next  = nexternal;
    m.norb  = ncore + nval + nexternal;
    m.core_first = 1;
    m.core_last  = ncore;
    m.val_first  = ncore + 1;
    m.val_last   = ncore + nval;
    m.ext_first  = ncore + nval + 1;
    m.ext_last   = ncore + nval + nexternal;
    return m;
}

/* ------------------------------------------------------------------ */
/* StateInfo                                                           */
/* ------------------------------------------------------------------ */

typedef struct { int nve, is; } info_t;

/* ------------------------------------------------------------------ */
/* Node helpers                                                        */
/* ------------------------------------------------------------------ */

/**
 * @brief Allocate a node for code @c (orb,nve,is) and set its admissibility
 *        @c status, with all arcs/weights/counts initialized to @ref ARC_NIL.
 * @param orb,nve,is the node code.
 * @param mo   MO classification (supplies @c norb and the valence boundary).
 * @param info target state (supplies the total @c nve and spin).
 * @return a freshly @c calloc'd node (caller/owner is whoever inserts it).
 *
 * Role: the single factory for construction-time nodes. Every generator and the
 * valence recursion routes allocation through here so the @c check_node
 * admissibility test is applied uniformly. A node is marked invalid
 * (@c status=0) when its code is out of range (negative spin, electron count
 * past the target, orbital past @c norb) or when it lies past the valence block
 * yet cannot still reach the target electron count (the @c orb >= val_last and
 * @c nve <= info->nve-3 clause); otherwise it is provisionally valid.
 *
 * DOMAIN (Shavitt/Paldus DRT; used by the paper, not re-derived there):
 * A node is inadmissible when it can no longer reach the target electron count.
 * Past the valence block (orb >= val_last) only the external region remains, and
 * each external orbital adds at most two electrons; if nve <= info->nve-3 the node
 * is already more than a double excitation short of the target, so no completion
 * exists and it is pruned. This is the Paldus/Shavitt distinct-row reachability
 * bound (refs [2,4,5]): only rows on a walk that completes to the target head
 * survive. The other clauses (is<0, nve>target, orb>norb) reject out-of-range
 * rows. The paper builds the DRT following Shavitt (Sec. 4.1) without restating
 * the bound.
 */
static node_t *node_new(int orb, int nve, int is,
                        const mo_t *mo, const info_t *info)
{
    node_t *n = calloc(1, sizeof(node_t));
    n->orb = orb;
    n->nve = nve;
    n->is  = is;
    /* check_node */
    /* Provisional validity: out-of-range codes and nodes too far past the
     * valence block to still reach info->nve are rejected here; reachability
     * pruning (@ref prune_unreachable) tightens this later. */
    if (is < 0 || nve > info->nve || orb > mo->norb ||
        (orb >= mo->val_last && nve <= info->nve - 3))
        n->status = 0;
    else
        n->status = 1;
    n->number_upper_arcs = ARC_NIL;
    for (int i = 0; i < 4; i++) {
        n->upper_arc[i]  = ARC_NIL;
        n->arc_weight[i] = ARC_NIL;
        n->lower_arc[i]  = ARC_NIL;
    }
    n->number_lower_arcs = ARC_NIL;
    return n;
}

/**
 * @brief Test whether node @p n can step (orbital @c orb+1) onto any node that
 *        already exists in table @p dest.
 * @param n    source node.
 * @param dest table of candidate successor nodes (e.g. the external-first set).
 * @return 1 if at least one of the four step successors is present in @p dest,
 *         else 0.
 *
 * Role: the join test used by @ref generate_valence_node at the valence/external
 * boundary -- a valence node is only kept if it connects into the
 * already-built downstream graph. The four candidate codes are the node's four
 * step successors at @c orb+1: the @c (ƒ¢nve,ƒ¢is) offsets @c (0,0), @c (1,+1),
 * @c (1,-1), @c (2,0).
 *
 * DOMAIN (Shavitt/Paldus DRT; used by the paper, not re-derived there):
 * The four (Delta nve, Delta is) offsets (0,0), (1,+1), (1,-1), (2,0) are the
 * canonical GUGA step vectors for the four one-orbital occupations e/u/d/f (step
 * codes 0..3): empty adds nothing, a single up/down electron adds one electron
 * and +/-1 to the spin, and double occupation adds two electrons at unchanged
 * spin. These are the Paldus row-difference steps between adjacent DRT levels
 * (refs [2,4,5]); step encoding cf. Sec. 4.2. The paper references the DRT
 * construction (Sec. 4.1) without restating the step convention.
 */
static int node_if_connect_area(const node_t *n, htab_t *dest)
{
    int targets[4][3] = {
        {n->orb+1, n->nve,   n->is  },
        {n->orb+1, n->nve+1, n->is+1},
        {n->orb+1, n->nve+1, n->is-1},
        {n->orb+1, n->nve+2, n->is  }
    };
    for (int i = 0; i < 4; i++)
        if (htab_lookup(dest, targets[i][0], targets[i][1], targets[i][2]))
            return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* DRT public struct                                                   */
/* ------------------------------------------------------------------ */

struct hexpr_drt {
    /* scalar parameters */
    info_t info;
    mo_t   mo;
    hexpr_drt_kind_t kind;

    /* sorted node array */
    node_t *nodes;    /* flat array, allocated */
    int     nnodes;
    int     top_index;

    /* path weights: lower_pw[i] is a flat array of size lower_pw_n[i] */
    int **lower_pw;
    int  *lower_pw_n;
    int **upper_pw;
    int  *upper_pw_n;

    /* packed CSF buffer: ncsfs * norb bytes, values 0..3 (E/U/D/F).
     * This is the CANONICAL GUGA-walk list (a superset for CSFSET sub-DRTs);
     * exposed publicly via hexpr_drt_walk_csfs / hexpr_drt_walk_ncsfs. */
    unsigned char *csf_packed;
    int            ncsfs_cached;

    /* Retained user input for CSFSET sub-DRTs (from_csfset): the n_csfs_input
     * rows, in the order supplied, n_csfs_input * norb bytes. NULL / 0 for
     * SOCI/FOCI (no input set). The user-facing hexpr_drt_csfs / _ncsfs /
     * csf_index / csf_steps speak this input system for CSFSET. */
    unsigned char *csf_input;
    int            n_csfs_input;
};

/* ------------------------------------------------------------------ */
/* Node generation                                                     */
/* ------------------------------------------------------------------ */

/* Store a node into h; if code already exists, do nothing. */
/**
 * @brief Deduplicating insert: store @p n unless its code already exists.
 * @param h table to insert into.
 * @param n node to store; if a node with the same code is already present, @p n
 *          is freed here and dropped, otherwise the table takes ownership.
 * @note The dedupe wrapper around @ref htab_insert used by every generator, so
 *       repeated @c node_new(...) of the same code is safe and leak-free.
 */
static void h_store(htab_t *h, node_t *n)
{
    if (htab_lookup(h, n->orb, n->nve, n->is)) {
        free(n);
        return;
    }
    htab_insert(h, n);
}

/**
 * @brief Emit the SOCI core-orbital nodes into @p h_nodes.
 * @param h_nodes destination construction table.
 * @param mo      MO classification (supplies the core orbital range).
 * @param info    target state (passed through to @ref node_new).
 *
 * Role: phase 1 of @ref hexpr_drt_create for the SOCI kind. Seeds the bottom
 * node @c (0,0,0) then walks each core orbital @c k, emitting a fixed family of
 * @c (nve,is) codes per orbital -- a special triple for the first orbital and a
 * four-node fan @c (nve-2,0)/(nve-2,2)/(nve-1,1)/(nve,0) for the rest, where
 * @c nve=2k. All inserts go through @ref h_store, so duplicates are dropped.
 *
 * DOMAIN (per paper, Honda & Noro): SOCI admits at most double holes in the core
 * region (Sec. 4.1, second-order CI). This emits the full (nve, is) occupancy fan
 * reachable by removing up to two electrons from the reference core -- the
 * nve/nve-1/nve-2 rows (with the distinct k==1 boundary case) -- i.e. the
 * second-order (double-hole) core space. See Sec. 4.1 (SOCI = double hole).
 */
static void soci_generate_core_node(htab_t *h_nodes, const mo_t *mo, const info_t *info)
{
    h_store(h_nodes, node_new(0, 0, 0, mo, info));   /* bottom */
    for (int k = mo->core_first; k <= mo->core_last; k++) {
        if (k == 1) {
            h_store(h_nodes, node_new(1, 0, 0, mo, info));
            h_store(h_nodes, node_new(1, 1, 1, mo, info));
            h_store(h_nodes, node_new(1, 2, 0, mo, info));
        } else {
            int nve = 2 * k;
            h_store(h_nodes, node_new(k, nve-2, 0, mo, info));
            h_store(h_nodes, node_new(k, nve-2, 2, mo, info));
            h_store(h_nodes, node_new(k, nve-1, 1, mo, info));
            h_store(h_nodes, node_new(k, nve,   0, mo, info));
        }
    }
}

/**
 * @brief Emit the FOCI core-orbital nodes into @p h_nodes.
 * @param h_nodes destination construction table.
 * @param mo      MO classification (supplies the core orbital range).
 * @param info    target state (passed through to @ref node_new).
 *
 * Role: phase 1 of @ref hexpr_drt_create for the FOCI kind. Same bottom-node
 * seed and per-orbital loop as @ref soci_generate_core_node, but emits a
 * STRICTLY SMALLER family per core orbital -- only @c (nve-1,1) and @c (nve,0)
 * (plus the @c k==1 special pair) -- i.e. FOCI omits the doubly-vacated
 * @c nve-2 rows SOCI keeps.
 *
 * DOMAIN (per paper, Honda & Noro): FOCI admits at most a single hole in the core
 * region (Sec. 4.1, first-order CI). It therefore drops the (nve-2,*) rows -- the
 * double-hole core configurations SOCI keeps -- retaining only up to one core
 * hole. This single- vs double-hole restriction is what distinguishes first-order
 * from second-order CI in core generation. See Sec. 4.1 (FOCI = single hole).
 */
static void foci_generate_core_node(htab_t *h_nodes, const mo_t *mo, const info_t *info)
{
    h_store(h_nodes, node_new(0, 0, 0, mo, info));   /* bottom */
    for (int k = mo->core_first; k <= mo->core_last; k++) {
        if (k == 1) {
            h_store(h_nodes, node_new(1, 1, 1, mo, info));
            h_store(h_nodes, node_new(1, 2, 0, mo, info));
        } else {
            int nve = 2 * k;
            h_store(h_nodes, node_new(k, nve-1, 1, mo, info));
            h_store(h_nodes, node_new(k, nve,   0, mo, info));
        }
    }
}

/**
 * @brief Emit the SOCI external-orbital nodes into @p h_nodes.
 * @param h_nodes destination construction table.
 * @param mo      MO classification (supplies the external orbital range).
 * @param info    target state: @c nve and spin seed the top of the external block.
 *
 * Role: phase 2 of @ref hexpr_drt_create for SOCI. Working down from the target
 * @c (nve,sp2), each external orbital @c k below the last gets up to five codes
 * (the @c nve-2 spin fan @c sp2/sp2-2/sp2+2, plus @c nve-1 at @c sp2}1), while
 * every @c k gets the on-target @c (nve,sp2) node. A trailing guard handles the
 * degenerate no-external case by seeding the top node directly.
 *
 * DOMAIN (per paper, Honda & Noro): SOCI admits at most two particles in the
 * external region (Sec. 4.1, second-order = double excitation). This generates
 * both the single-excitation and the (nve-2) double-excitation spin fan (sp2,
 * sp2+-2), covering up to two external excitations. See Sec. 4.1 (SOCI = double).
 */
static void soci_generate_external_node(htab_t *h_nodes, const mo_t *mo, const info_t *info)
{
    int nve = info->nve;
    int sp2 = info->is;
    for (int k = mo->ext_first; k <= mo->ext_last; k++) {
        if (k < mo->ext_last) {
            h_store(h_nodes, node_new(k, nve-2, sp2, mo, info));
            if (k < mo->ext_last - 1) {
                h_store(h_nodes, node_new(k, nve-2, sp2-2, mo, info));
                h_store(h_nodes, node_new(k, nve-2, sp2+2, mo, info));
            }
            h_store(h_nodes, node_new(k, nve-1, sp2-1, mo, info));
            h_store(h_nodes, node_new(k, nve-1, sp2+1, mo, info));
        }
        h_store(h_nodes, node_new(k, nve, sp2, mo, info));
    }
    /* no-external case: top node */
    if (mo->ext_last < mo->ext_first)
        h_store(h_nodes, node_new(mo->ext_last, nve, sp2, mo, info));
}

/**
 * @brief Emit the FOCI external-orbital nodes into @p h_nodes.
 * @param h_nodes destination construction table.
 * @param mo      MO classification (supplies the external orbital range).
 * @param info    target state: @c nve and spin seed the top of the external block.
 *
 * Role: phase 2 of @ref hexpr_drt_create for FOCI. Same descending loop and
 * no-external guard as @ref soci_generate_external_node, but each non-last
 * external orbital gets only the single-excitation @c nve-1 spin pair
 * @c (sp2-1, sp2+1) -- FOCI omits the @c nve-2 double-excitation fan.
 *
 * DOMAIN (per paper, Honda & Noro): FOCI admits at most one particle in the
 * external region (Sec. 4.1, first-order = single excitation). External
 * generation is therefore limited to the single-excitation (nve-1) pair and omits
 * the (nve-2) double-excitation fan that SOCI keeps -- the FOCI counterpart of the
 * SOCI external rule. See Sec. 4.1 (FOCI = single).
 */
static void foci_generate_external_node(htab_t *h_nodes, const mo_t *mo, const info_t *info)
{
    int nve = info->nve;
    int sp2 = info->is;
    for (int k = mo->ext_first; k <= mo->ext_last; k++) {
        if (k < mo->ext_last) {
            h_store(h_nodes, node_new(k, nve-1, sp2-1, mo, info));
            h_store(h_nodes, node_new(k, nve-1, sp2+1, mo, info));
        }
        h_store(h_nodes, node_new(k, nve, sp2, mo, info));
    }
    if (mo->ext_last < mo->ext_first)
        h_store(h_nodes, node_new(mo->ext_last, nve, sp2, mo, info));
}

/* Recursive valence node generator.
 * Returns the status of the node at (orb, nve, is). */
/**
 * @brief Recursively materialize the valence-block nodes reachable from
 *        @c (orb,nve,is) and return that node's validity.
 * @param orb,nve,is node code to generate.
 * @param dest    downstream (external-first) table used as the join target.
 * @param h_nodes construction table being populated (memoizes results).
 * @param mo,info MO classification and target state (passed to @ref node_new).
 * @return 1 if the node is valid (itself admissible and connected through to a
 *         valid descendant or, at the valence boundary, into @p dest), else 0.
 *
 * Role: phase 3 of @ref hexpr_drt_create. Builds the valence layer top-down by
 * a memoized DFS: @p h_nodes doubles as the visited set, so an already-generated
 * code returns its cached @c status immediately and the recursion stays linear
 * in the number of distinct codes.
 *
 * Control flow has three exits: (1) a code already present returns its status;
 * (2) an inadmissible node is recorded with @c status=0 and not expanded;
 * (3) at the last valence orbital @c val_last, validity is decided by
 * @ref node_if_connect_area against @p dest. Otherwise it recurses on the four
 * GUGA steps and is valid iff ANY child is valid (the @c e||u||d||f join).
 *
 * DOMAIN (per paper, Honda & Noro): The valence (active) region admits all
 * configurations (Sec. 4.1), so every valence orbital allows all four GUGA step
 * types -- (Delta nve, Delta is) = (0,0),(1,+1),(1,-1),(2,0), i.e. empty,
 * singly-up, singly-down, doubly occupied. A node is valid iff at least one step
 * reaches a valid descendant (the e||u||d||f join) because only walks that
 * complete to the top contribute CSFs; dead branches carry none. Unlike
 * core/external, no hole/particle count bounds the valence steps. See Sec. 4.1
 * (valence = all configurations); step encoding Sec. 4.2.
 */
static int generate_valence_node(int orb, int nve, int is,
                                 htab_t *dest, htab_t *h_nodes,
                                 const mo_t *mo, const info_t *info)
{
    /* Memoization: h_nodes is also the visited set, so a previously generated
     * code short-circuits with its cached status (no re-expansion). */
    node_t *existing = htab_lookup(h_nodes, orb, nve, is);
    if (existing)
        return existing->status;

    node_t *nn = node_new(orb, nve, is, mo, info);

    if (!nn->status) {
        htab_insert(h_nodes, nn);
        return 0;
    }

    int orb_last = mo->val_last;

    /* Valence/external boundary: keep this node only if it joins into the
     * already-built downstream graph held in `dest`. */
    if (orb == orb_last) {
        if (node_if_connect_area(nn, dest))
            nn->status = 1;
        else
            nn->status = 0;
        htab_insert(h_nodes, nn);
        return nn->status;
    }

    /* Interior valence orbital: recurse on the four steps; the node is valid
     * iff at least one child is reachable/valid. */
    int e = generate_valence_node(orb+1, nve,   is,   dest, h_nodes, mo, info);
    int u = generate_valence_node(orb+1, nve+1, is+1, dest, h_nodes, mo, info);
    int d = generate_valence_node(orb+1, nve+1, is-1, dest, h_nodes, mo, info);
    int f = generate_valence_node(orb+1, nve+2, is,   dest, h_nodes, mo, info);

    nn->status = (e || u || d || f) ? 1 : 0;
    htab_insert(h_nodes, nn);
    return nn->status;
}

/* ------------------------------------------------------------------ */
/* Prune unreachable nodes (BFS from bottom and from top)              */
/* ------------------------------------------------------------------ */

/**
 * @brief Mark every node not on a bottom-to-top path as invalid (@c status=0).
 * @param h   construction table holding all generated nodes (in @c order[]).
 * @param norb,nve,is code of the top node (the target state).
 *
 * Role: phase 5 of @ref hexpr_drt_create. A node belongs in the final DRT only
 * if it lies on at least one complete walk, i.e. it is reachable forward from
 * the bottom @c (0,0,0) AND backward from the top @c (norb,nve,is). This runs
 * two breadth-first searches over the node graph and clears @c status on any
 * node missing from EITHER reachable set.
 *
 * Two passes share one @c queue buffer: the forward pass follows the four step
 * successors @c (ƒ¢orb=+1) seeded at the bottom; the backward pass follows the
 * four predecessors @c (ƒ¢orb=-1) seeded at the top. Reachability is tracked in
 * @c rfb (from-bottom) and @c rft (from-top) flag arrays indexed by @c order[]
 * position. Early-return guards bail out cleanly (freeing scratch) if the
 * bottom or top code is absent.
 *
 * @note Pure graph reachability -- no GUGA-specific intent; the step/predecessor
 *       offsets mirror those documented at @ref node_if_connect_area. The inner
 *       @c order[] linear scan to map a node pointer back to its index makes
 *       this O(n^2) in node count, which is acceptable for DRT sizes here.
 *
 * @note Both BFS passes only cross into a neighbor @c nb whose own @c status
 *       is already 1 (set by the generation phase before this runs). Without
 *       that check, a node reachable only through an already-invalid
 *       intermediate (one @ref generate_valence_node marked @c status=0) would
 *       still be marked reachable here, then survive into the final sorted
 *       @c nodes[] array with no real arc back to that intermediate in the
 *       arc-wiring phase that follows (arcs are only wired among survivors).
 *       Such a node has an incoming arc but no real path to the top, and
 *       @ref store_number_upper_arcs's @c imax(child_n,1) floor (there to
 *       count the top node's own single walk) would then also credit that
 *       dead node's parent with one phantom walk -- inflating
 *       @c nodes[0].number_upper_arcs (hexpr_drt_ncsfs()/walk_ncsfs() for
 *       SOCI/FOCI) past what @ref hexpr_drt_mk_csf_list actually enumerates,
 *       and past what @c csf_packed is allocated to hold.
 *
 * @note Symmetrically, each loop also skips expanding a dequeued @c cur
 *       whose own @c status is already 0 -- an invalid node radiates no
 *       reachability of its own, not just failing to be crossed into as a
 *       neighbor. This matters for the seed itself: @ref node_new can mark
 *       the bottom @c (0,0,0) invalid (its own "past the valence block and
 *       too far short of the target" bound applies to the origin whenever
 *       @c ncore=nval=0, since then @c val_last<=0 makes @c orb>=val_last
 *       trivially true even at @c orb=0). Before this check, the forward
 *       BFS still unconditionally seeded @c rfb at the bottom and expanded
 *       its successors regardless of that invalidity, so a state @p info
 *       that genuinely has no valid walk from the vacuum (e.g. requesting
 *       more particles in the external region than SOCI's <=2 / FOCI's <=1
 *       admits, when there is no core or valence reference to supply any of
 *       them) could still find a spurious top-reaching path and return a
 *       DRT built one orbital short of the true bottom, instead of
 *       correctly finding no top-reaching walk at all.
 */
static void prune_unreachable(htab_t *h, int norb, int nve, int is)
{
    /* reachable from bottom */
    int *rfb = calloc(h->order_size, sizeof(int));
    /* queue as dynamic array */
    int *queue = malloc(h->order_size * sizeof(int));
    int qhead = 0, qtail = 0;

    /* find bottom index */
    node_t *bottom = htab_lookup(h, 0, 0, 0);
    if (!bottom) { free(rfb); free(queue); return; }

    /* map node pointer -> order index */
    /* We'll use a secondary hash keyed on pointer... or just use h->order[] scan.
     * For correctness, BFS on actual codes stored in order[] */

    /* Use a separate visited-flag array indexed by h->order position */
    int n = h->order_size;
    for (int i = 0; i < n; i++) rfb[i] = 0;

    /* find index of bottom in order[] */
    int bottom_idx = -1;
    for (int i = 0; i < n; i++) {
        if (h->order[i] == bottom) { bottom_idx = i; break; }
    }
    if (bottom_idx < 0) { free(rfb); free(queue); return; }

    rfb[bottom_idx] = 1;
    queue[qtail++] = bottom_idx;
    while (qhead < qtail) {
        int ci = queue[qhead++];
        node_t *cur = h->order[ci];
        if (!cur->status) continue;   /* an invalid node radiates no reachability of
                                        * its own either -- see the doc comment above
                                        * (this is what makes the bottom-seed case
                                        * correct) */
        int successors[4][3] = {
            {cur->orb+1, cur->nve,   cur->is  },
            {cur->orb+1, cur->nve+1, cur->is+1},
            {cur->orb+1, cur->nve+1, cur->is-1},
            {cur->orb+1, cur->nve+2, cur->is  }
        };
        for (int s = 0; s < 4; s++) {
            node_t *nb = htab_lookup(h, successors[s][0], successors[s][1], successors[s][2]);
            if (!nb || !nb->status) continue;   /* don't route reachability through an
                                                  * already-invalid node -- see note below */
            /* find index */
            for (int j = 0; j < n; j++) {
                if (h->order[j] == nb && !rfb[j]) {
                    rfb[j] = 1;
                    queue[qtail++] = j;
                    break;
                }
            }
        }
    }

    /* reachable from top */
    int *rft = calloc(n, sizeof(int));
    node_t *top = htab_lookup(h, norb, nve, is);
    if (!top) { free(rfb); free(rft); free(queue); return; }
    int top_idx = -1;
    for (int i = 0; i < n; i++) {
        if (h->order[i] == top) { top_idx = i; break; }
    }
    if (top_idx < 0) { free(rfb); free(rft); free(queue); return; }

    qhead = qtail = 0;
    rft[top_idx] = 1;
    queue[qtail++] = top_idx;
    while (qhead < qtail) {
        int ci = queue[qhead++];
        node_t *cur = h->order[ci];
        if (!cur->status) continue;   /* symmetric with the forward pass above */
        int predecessors[4][3] = {
            {cur->orb-1, cur->nve,   cur->is  },
            {cur->orb-1, cur->nve-1, cur->is-1},
            {cur->orb-1, cur->nve-1, cur->is+1},
            {cur->orb-1, cur->nve-2, cur->is  }
        };
        for (int s = 0; s < 4; s++) {
            node_t *nb = htab_lookup(h, predecessors[s][0], predecessors[s][1], predecessors[s][2]);
            if (!nb || !nb->status) continue;   /* see note at the forward pass above */
            for (int j = 0; j < n; j++) {
                if (h->order[j] == nb && !rft[j]) {
                    rft[j] = 1;
                    queue[qtail++] = j;
                    break;
                }
            }
        }
    }

    /* mark unreachable nodes as status=0 so they'll be filtered out */
    for (int i = 0; i < n; i++) {
        if (!rfb[i] || !rft[i])
            h->order[i]->status = 0;
    }

    free(rfb); free(rft); free(queue);
}

/* ------------------------------------------------------------------ */
/* Sorting                                                             */
/* ------------------------------------------------------------------ */

/**
 * @brief qsort comparator ordering nodes lexicographically by @c (orb,nve,is).
 * @param a,b pointers to @c node_t* (note: pointer-to-pointer).
 * @return negative/zero/positive for the (orb, then nve, then is) ordering.
 * @note Establishes the canonical node index order of the assembled DRT.
 */
static int node_cmp(const void *a, const void *b)
{
    const node_t *na = *(const node_t **)a;
    const node_t *nb = *(const node_t **)b;
    if (na->orb != nb->orb) return na->orb - nb->orb;
    if (na->nve != nb->nve) return na->nve - nb->nve;
    return na->is - nb->is;
}

/* ------------------------------------------------------------------ */
/* store_number_upper_arcs (recursive)                                 */
/* ------------------------------------------------------------------ */

/** @brief Larger of two integers. @param a,b operands @return max(a,b). */
static int imax(int a, int b) { return a > b ? a : b; }

/**
 * @brief Recursively compute each node's upper-arc count and per-arc weights,
 *        memoized in the node fields.
 * @param nodes     flat sorted node array.
 * @param node      index of the node to process.
 * @param top_index index of the top node (recursion base case).
 * @return @c number_upper_arcs for @p node (number of distinct upward walks to
 *         the top).
 *
 * Role: phase 8 of @ref hexpr_drt_create. Fills @c number_upper_arcs for every
 * node and assigns @c arc_weight[i] along the way; those weights are what
 * @ref tree_search_lower later sums into lower path weights, and what the
 * downstream coupling code reads via @ref hexpr_drt_arc_weight.
 *
 * Structure: a memoized post-order recursion (an already-set
 * @c number_upper_arcs short-circuits; the top node anchors at 0). For each
 * present upper arc it recurses, adds @c max(child_count,1) to this node's
 * count, records the running @c weight into @c arc_weight[i], then advances
 * @c weight by the child's upper-arc count -- so @c arc_weight[i] is the offset
 * of arc @c i within this node's enumerated upward walks.
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * This sets up the reverse-lexical arc weights that csf_lookup.c later sums to
 * address CSFs (Paldus-Shavitt, refs [2,4,5]). arc_weight[i] is the prefix sum of
 * the sibling subtree counts, i.e. the number of canonical walks passing through
 * lexically-earlier arcs at this node, so summing arc weights along a walk yields
 * its canonical index. The max(child,1) guard counts the top node (one CSF though
 * its subtree count is 0). This is the defining construction of the arc-weight
 * convention; the paper builds the DRT "following Shavitt" (Sec. 4.1) without
 * restating it.
 */
static int store_number_upper_arcs(node_t *nodes, int node, int top_index)
{
    if (nodes[node].number_upper_arcs != ARC_NIL)
        return nodes[node].number_upper_arcs;

    if (node == top_index) {
        nodes[node].number_upper_arcs = 0;
        return 0;
    }

    int num = 0, weight = 0;
    for (int i = 0; i < 4; i++) {
        int v = nodes[node].upper_arc[i];
        if (v != ARC_NIL) {
            int child_n = store_number_upper_arcs(nodes, v, top_index);
            num += imax(child_n, 1);
            /* arc_weight[i] = offset of this arc within the node's upward walks;
             * advance the running weight by the child's own walk count. */
            nodes[node].arc_weight[i] = weight;
            weight += nodes[v].number_upper_arcs;
        }
    }
    nodes[node].number_upper_arcs = num;
    return num;
}

/* ------------------------------------------------------------------ */
/* store_number_lower_arcs (recursive)                                 */
/* ------------------------------------------------------------------ */

/**
 * @brief Recursively compute each node's lower-arc count, memoized in the node.
 * @param nodes flat sorted node array.
 * @param node  index of the node to process.
 * @return @c number_lower_arcs for @p node (number of distinct downward walks to
 *         the bottom).
 *
 * Role: phase 8 of @ref hexpr_drt_create, the downward counterpart of
 * @ref store_number_upper_arcs. Same memoized post-order recursion, anchored at
 * the bottom node (index 0) instead of the top, summing @c max(child,1) over the
 * present lower arcs. Unlike the upper pass it assigns no per-arc weight -- only
 * the count is needed downstream (exposed via @ref hexpr_drt_number_lower_arcs).
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * The lower-arc count is the number of downward walks reachable from a node, the
 * companion to the upper-arc count, used to size the lower path weights. The
 * bottom node (index 0) is the anchor because lower walks are enumerated upward
 * from it, mirroring how upper walks anchor at the top -- the two ends of the
 * Shavitt walk space (refs [4,5]). See csf_lookup.c for the addressing use; the
 * paper references the DRT (Sec. 4.1) without restating it.
 */
static int store_number_lower_arcs(node_t *nodes, int node)
{
    if (nodes[node].number_lower_arcs != ARC_NIL)
        return nodes[node].number_lower_arcs;

    if (node == 0) {
        nodes[node].number_lower_arcs = 0;
        return 0;
    }

    int num = 0;
    for (int i = 0; i < 4; i++) {
        int v = nodes[node].lower_arc[i];
        if (v != ARC_NIL)
            num += imax(store_number_lower_arcs(nodes, v), 1);
    }
    nodes[node].number_lower_arcs = num;
    return num;
}

/* ------------------------------------------------------------------ */
/* Path weight computation                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    int *data;
    int  size;
    int  cap;
} ivec_t;

/**
 * @brief Append an int to a growable int vector, doubling capacity as needed.
 * @param v   target vector.
 * @param val value to append.
 * @note May @c realloc and move @c v->data. Backing store for the per-node
 *       path-weight lists assembled during construction.
 */
static void ivec_push(ivec_t *v, int val)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 4;
        v->data = realloc(v->data, v->cap * sizeof(int));
    }
    v->data[v->size++] = val;
}

/**
 * @brief Accumulate lower path weights into every node by walking upward and
 *        summing arc weights.
 * @param nodes     flat sorted node array.
 * @param inode     current node index (start at bottom = 0).
 * @param top_index index of the top node (recursion stop).
 * @param lw        array of per-node int vectors; @c lw[j] receives one entry
 *                  per upward walk that reaches node @c j.
 * @param weight    running arc-weight sum on the path from the bottom to @p inode.
 *
 * Role: phase 9 of @ref hexpr_drt_create. A DFS over @c upper_arc that, on each
 * descent into successor @c next, pushes the cumulative weight @c weight+arc_weight
 * onto @c lw[next] before recursing -- so each node collects the weights of all
 * lower (bottom-side) partial walks ending at it. These become the lists served
 * by @ref hexpr_drt_lower_path_weights. Arcs with a NIL node or NIL weight are
 * skipped.
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * A lower path weight indexes a partial walk in the lower (below-node) portion of
 * the DRT: it is the running sum of arc weights along the upward path, i.e. the
 * reverse-lexical offset of that partial walk among the lower walks (Shavitt,
 * refs [4,5]). Paired with an upper path weight it identifies one full CSF (cf.
 * expand_expr, csf_lookup.c). The paper references the DRT (Sec. 4.1) as is.
 */
static void tree_search_lower(node_t *nodes, int inode, int top_index,
                              ivec_t *lw, int weight)
{
    if (inode == top_index) return;
    for (int step = 0; step < 4; step++) {
        int next = nodes[inode].upper_arc[step];
        int aw   = nodes[inode].arc_weight[step];
        if (next != ARC_NIL && aw != ARC_NIL) {
            int nw = weight + aw;
            ivec_push(&lw[next], nw);
            tree_search_lower(nodes, next, top_index, lw, nw);
        }
    }
}

/**
 * @brief Enumerate the upper path weights of one node: the arc-weight sums of
 *        all walks from it up to the top.
 * @param nodes     flat sorted node array.
 * @param inode     node whose upper path weights are being collected.
 * @param top_index index of the top node (recursion base case).
 * @param weight    running arc-weight sum accumulated so far.
 * @param out       int vector receiving one weight per completed walk to the top.
 *
 * Role: phase 10 of @ref hexpr_drt_create, called once per node. Differs from
 * @ref tree_search_lower in WHERE it records: it emits only at the top (one
 * entry per full upper walk) into a single @p out vector, rather than depositing
 * partial sums into every visited node. Results back @ref hexpr_drt_upper_path_weights.
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * An upper path weight is the reverse-lexical offset of a partial walk among the
 * upper (above-node) walks through a node (Shavitt, refs [4,5]). A node's lower
 * and upper path weights combine additively into the full walk's canonical index
 * (lower offset + upper offset), which is the CSF index used downstream. The
 * paper builds the DRT following Shavitt (Sec. 4.1) without restating the scheme.
 */
static void tree_search_upper(node_t *nodes, int inode, int top_index,
                              int weight, ivec_t *out)
{
    if (inode == top_index) {
        ivec_push(out, weight);
        return;
    }
    for (int step = 0; step < 4; step++) {
        int next = nodes[inode].upper_arc[step];
        int aw   = nodes[inode].arc_weight[step];
        if (next != ARC_NIL && aw != ARC_NIL)
            tree_search_upper(nodes, next, top_index, weight + aw, out);
    }
}

/* ------------------------------------------------------------------ */
/* CSF list generation                                                 */
/* ------------------------------------------------------------------ */

static const char WALK_CHAR[4] = {'e','u','d','f'};

typedef struct {
    char **list;
    int    size;
    int    cap;
    char  *csf_buf;   /* reused buffer of length norb */
    int    norb;
    int    depth;
} csf_ctx_t;

/**
 * @brief Snapshot the current walk buffer as a new NUL-terminated CSF string and
 *        append it to the context's list.
 * @param ctx enumeration context; @c ctx->csf_buf[0..norb-1] holds the current
 *            walk's step characters.
 * @note Copies @c norb bytes into a fresh @c malloc'd @c norb+1 string (caller of
 *       @ref hexpr_drt_mk_csf_list owns and frees each). Grows @c ctx->list by
 *       doubling. Called at each completed walk by @ref tree_search_csf.
 */
static void csf_push(csf_ctx_t *ctx)
{
    if (ctx->size == ctx->cap) {
        ctx->cap = ctx->cap ? ctx->cap * 2 : 64;
        ctx->list = realloc(ctx->list, ctx->cap * sizeof(char *));
    }
    char *s = malloc(ctx->norb + 1);
    memcpy(s, ctx->csf_buf, ctx->norb);
    s[ctx->norb] = '\0';
    ctx->list[ctx->size++] = s;
}

/**
 * @brief Enumerate every bottom-to-top walk, emitting one CSF step-string per
 *        walk into the context list.
 * @param nodes     flat sorted node array.
 * @param inode     current node index (start at bottom = 0).
 * @param top_index index of the top node (a completed walk).
 * @param ctx       enumeration context (walk buffer + output list + depth).
 *
 * Role: the recursive core of @ref hexpr_drt_mk_csf_list. A DFS over
 * @c upper_arc that records the step character @c WALK_CHAR[step] ('e','u','d','f')
 * at the current @c depth, recurses, then backtracks (@c depth--). Reaching the
 * top finalizes the buffered walk via @ref csf_push. The shared @c csf_buf is
 * mutated in place and restored on unwind, so the enumeration is O(walks*norb)
 * without per-walk buffer allocation.
 *
 * DOMAIN (per paper, Honda & Noro): The step characters e,u,d,f encode the four
 * one-orbital occupations of the genealogical CSF representation and map to step
 * codes 0,1,2,3 as fixed by the paper (Sec. 4.2, cf. the segment-walk table):
 * e=empty, u=singly-occupied spin-up, d=singly-occupied spin-down, f=doubly
 * (full) occupied. This ordering is the API-level convention the explicit-CSF
 * input mode uses. See Sec. 4.2.
 */
static void tree_search_csf(node_t *nodes, int inode, int top_index, csf_ctx_t *ctx)
{
    if (inode == top_index) {
        csf_push(ctx);
        return;
    }
    for (int step = 0; step < 4; step++) {
        int next = nodes[inode].upper_arc[step];
        if (next != ARC_NIL) {
            ctx->csf_buf[ctx->depth] = WALK_CHAR[step];
            ctx->depth++;
            tree_search_csf(nodes, next, top_index, ctx);
            ctx->depth--;
        }
    }
}

/* ------------------------------------------------------------------ */
/* hexpr_drt_create                                                    */
/* ------------------------------------------------------------------ */

/**
 * @brief Build a parametric SOCI/FOCI Distinct Row Table from MO counts and a
 *        target state.
 * @param kind    @c HEXPR_DRT_SOCI or @c HEXPR_DRT_FOCI (selects the generators).
 * @param nve     target number of valence electrons (the top node's @c nve).
 * @param spin_x2 target spin times two (the top node's @c is).
 * @param ncore,nval,next core / valence / external orbital counts.
 * @return a fully assembled ::hexpr_drt_t (release with
 *         @ref hexpr_drt_destroy), or NULL if no valid top node exists for the
 *         requested state.
 *
 * The construction pipeline, in the numbered phases below: generate core nodes
 * (@ref soci_generate_core_node / @ref foci_generate_core_node), generate
 * external nodes, grow the valence layer between them
 * (@ref generate_valence_node), drop invalid and unreachable nodes
 * (@ref prune_unreachable), sort survivors into canonical order
 * (@ref node_cmp) into a flat array, wire each node's upper/lower arcs, compute
 * arc counts and weights (@ref store_number_upper_arcs /
 * @ref store_number_lower_arcs), enumerate the lower/upper path weights
 * (@ref tree_search_lower / @ref tree_search_upper), and finally pack the
 * canonical CSF list. Several transient @c htab_t tables are freed before return.
 *
 * @note The valence layer is bracketed by two collected node sets -- the
 *       core-last nodes (seeds) and the external-first nodes (join targets) --
 *       so the recursion only has to bridge the valence block.
 *
 * DOMAIN (per paper, Honda & Noro): The overall SOCI/FOCI construction expands
 * the orbital-partition specification (inactive / core / valence / external) into
 * the complete CSF list of the corresponding first- or second-order CI space and
 * builds the DRT over it (Sec. 4.1). The core, valence, and external node
 * generators (documented individually above) implement the per-region occupancy
 * bounds -- single holes/particles for FOCI, double for SOCI, all configurations
 * in the valence region -- and this routine assembles them into one DRT. See
 * Sec. 4.1 (FOCI/SOCI partition).
 */
hexpr_drt_t *hexpr_drt_create(hexpr_drt_kind_t kind,
                              int nve, int spin_x2,
                              int ncore, int nval, int next)
{
    mo_t   mo   = mo_init(ncore, nval, next);
    info_t info = {nve, spin_x2};

    htab_t *h_nodes = htab_new();

    /* 1. generate core nodes */
    if (kind == HEXPR_DRT_SOCI)
        soci_generate_core_node(h_nodes, &mo, &info);
    else
        foci_generate_core_node(h_nodes, &mo, &info);

    /* collect nodes where orb == core_range.last (or bottom if ncore==0) */
    htab_t *nodes_core_last = htab_new();
    for (int i = 0; i < h_nodes->order_size; i++) {
        node_t *nd = h_nodes->order[i];
        if (nd->orb == mo.core_last)
            h_store(nodes_core_last, node_new(nd->orb, nd->nve, nd->is, &mo, &info));
    }

    /* 2. generate external nodes */
    if (kind == HEXPR_DRT_SOCI)
        soci_generate_external_node(h_nodes, &mo, &info);
    else
        foci_generate_external_node(h_nodes, &mo, &info);

    /* collect nodes where orb == ext_range.first (or top if next==0) */
    htab_t *nodes_ext_first = htab_new();
    for (int i = 0; i < h_nodes->order_size; i++) {
        node_t *nd = h_nodes->order[i];
        int want_orb = (mo.next > 0) ? mo.ext_first : mo.ext_last;
        if (nd->orb == want_orb)
            h_store(nodes_ext_first, node_new(nd->orb, nd->nve, nd->is, &mo, &info));
    }

    /* 3. generate valence nodes from each core_last node */
    for (int i = 0; i < nodes_core_last->order_size; i++) {
        node_t *nd = nodes_core_last->order[i];
        generate_valence_node(nd->orb+1, nd->nve,   nd->is,   nodes_ext_first, h_nodes, &mo, &info);
        generate_valence_node(nd->orb+1, nd->nve+1, nd->is+1, nodes_ext_first, h_nodes, &mo, &info);
        generate_valence_node(nd->orb+1, nd->nve+1, nd->is-1, nodes_ext_first, h_nodes, &mo, &info);
        generate_valence_node(nd->orb+1, nd->nve+2, nd->is,   nodes_ext_first, h_nodes, &mo, &info);
    }

    /* 4. remove invalid nodes */
    for (int i = 0; i < h_nodes->order_size; i++)
        if (!h_nodes->order[i]->status)
            h_nodes->order[i]->status = 0;  /* already marked */

    /* 5. prune unreachable */
    prune_unreachable(h_nodes, mo.norb, nve, spin_x2);

    /* 6. collect valid nodes and sort */
    int valid_count = 0;
    for (int i = 0; i < h_nodes->order_size; i++)
        if (h_nodes->order[i]->status) valid_count++;

    node_t **tmp_ptrs = malloc(valid_count * sizeof(node_t *));
    int vi = 0;
    for (int i = 0; i < h_nodes->order_size; i++)
        if (h_nodes->order[i]->status)
            tmp_ptrs[vi++] = h_nodes->order[i];
    qsort(tmp_ptrs, valid_count, sizeof(node_t *), node_cmp);

    node_t *nodes = malloc(valid_count * sizeof(node_t));
    for (int i = 0; i < valid_count; i++)
        nodes[i] = *tmp_ptrs[i];
    free(tmp_ptrs);

    /* find top_index */
    int top_index = -1;
    for (int i = 0; i < valid_count; i++) {
        if (nodes[i].orb == mo.norb && nodes[i].nve == nve && nodes[i].is == spin_x2) {
            top_index = i;
            break;
        }
    }
    if (top_index < 0) {
        hexpr_set_error(
            "hexpr_drt_create: no node exists for the requested target "
            "state (kind=%d, nve=%d, spin_x2=%d, ncore=%d, nval=%d, next=%d)",
            (int)kind, nve, spin_x2, ncore, nval, next);
        free(nodes);
        htab_free(nodes_core_last);
        htab_free(nodes_ext_first);
        htab_free(h_nodes);
        return NULL;
    }

    /* 7. build code2num mapping and set upper/lower arcs */
    /* (reuse a simple flat array of sorted nodes to binary-search) */
    /* Actually, just rebuild a lookup htab on the sorted nodes */
    htab_t *code2num_h = htab_new();
    for (int i = 0; i < valid_count; i++) {
        /* store a node_t* whose address encodes the index; we abuse the pointer */
        node_t *proxy = calloc(1, sizeof(node_t));
        proxy->orb = nodes[i].orb;
        proxy->nve = nodes[i].nve;
        proxy->is  = nodes[i].is;
        /* store index in number_upper_arcs field (arbitrary reuse) */
        proxy->number_upper_arcs = i;
        proxy->status = 1;
        for (int w = 0; w < 4; w++) { proxy->upper_arc[w] = ARC_NIL; proxy->arc_weight[w] = ARC_NIL; proxy->lower_arc[w] = ARC_NIL; }
        proxy->number_lower_arcs = ARC_NIL;
        htab_insert(code2num_h, proxy);
    }

    /* Wire arcs reciprocally: for each node, look up its four step successors
     * (codes succ[w]) via code2num_h to get the target index jnode; record it as
     * this node's upper_arc[w] and, symmetrically, set the target's lower_arc[w]
     * back to this node. The proxy's stashed index (number_upper_arcs) is the
     * code -> sorted-index map. */
    for (int inode = 0; inode < valid_count; inode++) {
        node_t *nd = &nodes[inode];
        int succ[4][3] = {
            {nd->orb+1, nd->nve,   nd->is  },
            {nd->orb+1, nd->nve+1, nd->is+1},
            {nd->orb+1, nd->nve+1, nd->is-1},
            {nd->orb+1, nd->nve+2, nd->is  }
        };
        for (int w = 0; w < 4; w++) {
            node_t *p = htab_lookup(code2num_h, succ[w][0], succ[w][1], succ[w][2]);
            int jnode = p ? p->number_upper_arcs : ARC_NIL;
            nd->upper_arc[w] = jnode;
            if (jnode != ARC_NIL)
                nodes[jnode].lower_arc[w] = inode;
        }
    }

    /* 8. store_number_upper_arcs and store_number_lower_arcs */
    store_number_upper_arcs(nodes, 0, top_index);
    /* NOTE: this assumes the top node sorts last in nodes[] (node_cmp's
     * ascending (orb,nve,is) order), unchecked -- the mirror of the
     * nodes[0]-is-bottom assumption asserted below, but for the top end.
     * hexpr_drt_from_csfset's analogous call passes top_index explicitly
     * instead of assuming a position; this one does not. Left as found,
     * not fixed here. */
    store_number_lower_arcs(nodes, valid_count - 1);

    /* 9. set_lower_path_weights */
    ivec_t *lower_pw = calloc(valid_count, sizeof(ivec_t));
    ivec_push(&lower_pw[0], 0);
    tree_search_lower(nodes, 0, top_index, lower_pw, 0);

    /* 10. set_upper_path_weights */
    ivec_t *upper_pw = calloc(valid_count, sizeof(ivec_t));
    for (int i = 0; i < valid_count; i++) {
        if (nodes[i].number_upper_arcs != ARC_NIL)
            tree_search_upper(nodes, i, top_index, 0, &upper_pw[i]);
    }

    /* Assemble result */
    hexpr_drt_t *drt = calloc(1, sizeof(hexpr_drt_t));
    drt->info      = info;
    drt->mo        = mo;
    drt->kind      = kind;
    drt->nodes     = nodes;
    drt->nnodes    = valid_count;
    drt->top_index = top_index;

    drt->lower_pw   = malloc(valid_count * sizeof(int *));
    drt->lower_pw_n = malloc(valid_count * sizeof(int));
    drt->upper_pw   = malloc(valid_count * sizeof(int *));
    drt->upper_pw_n = malloc(valid_count * sizeof(int));
    for (int i = 0; i < valid_count; i++) {
        drt->lower_pw[i]   = lower_pw[i].data;
        drt->lower_pw_n[i] = lower_pw[i].size;
        drt->upper_pw[i]   = upper_pw[i].data;
        drt->upper_pw_n[i] = upper_pw[i].size;
    }
    free(lower_pw);
    free(upper_pw);

    /* htab_free frees all node_t* via the bucket chain */
    htab_free(code2num_h);
    htab_free(nodes_core_last);
    htab_free(nodes_ext_first);
    htab_free(h_nodes);

    /* build packed CSF buffer */
    {
        static const unsigned char CHAR_TO_STEP[256] = {
            ['e'] = 0, ['u'] = 1, ['d'] = 2, ['f'] = 3
        };
        char **csfs = NULL;
        int ncsf = hexpr_drt_mk_csf_list(drt, &csfs);
        int norb  = drt->mo.norb;

        /* Construction-time invariant: nodes[0] must be the walk origin
         * (0,0,0). Every accessor that treats "index 0 in the sorted node
         * array" as "the bottom" -- hexpr_drt_ncsfs()/walk_ncsfs() via
         * nodes[0].number_upper_arcs, hexpr_drt_walk_csf_index/_steps'
         * descent starting at node 0, and hexpr_drt_mk_csf_list's DFS above
         * -- depends on this and none of them checks it. The ncsf-vs-
         * number_upper_arcs guard just below cannot catch a violation of
         * this on its own: mk_csf_list and store_number_upper_arcs both
         * start from the SAME nodes[0], so if it is the wrong node both
         * sides shift together and still agree with each other. node_new
         * can mark the true (0,0,0) inadmissible (its "past the valence
         * block, too far short of target" bound applies to the origin
         * whenever ncore=nval=0), and prior to the prune_unreachable fix
         * above, a spurious path from some OTHER surviving node could still
         * reach the top, leaving that other node -- not (0,0,0) -- as the
         * smallest sorted survivor. */
        if (drt->nodes[0].orb != 0 || drt->nodes[0].nve != 0 || drt->nodes[0].is != 0) {
            hexpr_set_error(
                "hexpr_drt_create: internal invariant violated (kind=%d, "
                "nve=%d, spin_x2=%d, ncore=%d, nval=%d, next=%d): nodes[0] "
                "is (orb=%d,nve=%d,is=%d), not the walk origin (0,0,0); "
                "refusing to build a DRT every accessor would misread",
                (int)kind, nve, spin_x2, ncore, nval, next,
                drt->nodes[0].orb, drt->nodes[0].nve, drt->nodes[0].is);
            for (int i = 0; i < ncsf; i++) free(csfs[i]);
            free(csfs);
            hexpr_drt_free(drt);
            return NULL;
        }

        /* Construction-time consistency guard: nodes[0].number_upper_arcs
         * (what hexpr_drt_ncsfs()/hexpr_drt_walk_ncsfs() report for
         * SOCI/FOCI) is computed independently, via store_number_upper_arcs's
         * imax(child_n, 1) recursion, from the walk count mk_csf_list just
         * enumerated by actually following upper_arc to the top. The two
         * must agree; a mismatch means a node survived prune_unreachable
         * without a real path to the top (dead node), which inflates the
         * reported count past what csf_packed below is sized to hold. Fail
         * loudly here rather than returning a DRT whose own count doesn't
         * match its own CSF buffer. */
        if (ncsf != drt->nodes[0].number_upper_arcs) {
            hexpr_set_error(
                "hexpr_drt_create: canonical CSF count mismatch (kind=%d, "
                "nve=%d, spin_x2=%d, ncore=%d, nval=%d, next=%d): "
                "mk_csf_list enumerated %d walk(s) but "
                "nodes[0].number_upper_arcs reports %d; refusing to build "
                "an internally inconsistent DRT",
                (int)kind, nve, spin_x2, ncore, nval, next, ncsf,
                drt->nodes[0].number_upper_arcs);
            for (int i = 0; i < ncsf; i++) free(csfs[i]);
            free(csfs);
            hexpr_drt_free(drt);
            return NULL;
        }

        drt->csf_packed    = malloc((size_t)ncsf * (size_t)norb);
        drt->ncsfs_cached  = ncsf;
        for (int i = 0; i < ncsf; i++) {
            for (int j = 0; j < norb; j++)
                drt->csf_packed[i * norb + j] = CHAR_TO_STEP[(unsigned char)csfs[i][j]];
            free(csfs[i]);
        }
        free(csfs);
    }

    return drt;
}

/* ------------------------------------------------------------------ */
/* hexpr_drt_from_csfset                                               */
/* ------------------------------------------------------------------ */

/**
 * @brief Build a sub-DRT from an explicit set of input CSFs (CSFSET kind).
 * @param norb   number of orbitals (length of each CSF row).
 * @param csfs   @p n_csfs by @p norb byte matrix of step codes 0..3, row-major;
 *               must be non-NULL. Retained verbatim for the user-facing accessors.
 * @param n_csfs number of input CSFs (>= 1).
 * @return a ::hexpr_drt_t of kind @c HEXPR_DRT_CSFSET (release with
 *         @ref hexpr_drt_destroy), or NULL on a validation error (with the
 *         cause recorded via @c hexpr_set_error).
 *
 * Role: the alternative builder to @ref hexpr_drt_create. Instead of generating
 * a parametric node space, it walks each supplied CSF orbital by orbital,
 * accumulating the running code @c (orb,nve,is) from the per-step deltas
 * @c DELTA_NVE / @c DELTA_IS and inserting each distinct node visited. The union
 * of all walks forms the graph, which is then sorted, arc-wired, and weighted by
 * the SAME machinery @ref hexpr_drt_create uses (@ref store_number_upper_arcs,
 * @ref tree_search_lower, ...).
 *
 * Two consistency checks gate construction: every step must be in range and stay
 * on an admissible node (the @c a/b/c constraints below), and all CSFs must end
 * at the SAME top code (else they do not share one DRT). The input set is copied
 * into @c csf_input so the user-facing list keeps the caller's order even though
 * the canonical packed list may be a superset.
 *
 * DOMAIN (Shavitt/Paldus DRT; used by the paper, not re-derived there):
 * The three guards are the Paldus ABC non-negativity conditions (Paldus,
 * J.Chem.Phys. 61, ref [2]). In the Paldus row parametrization a = doubly-
 * occupied count, b = singly-occupied (spin) count, c = empty count, expressed
 * via (orb,nve,is) as b=is, a=(nve-is)/2, c=orb-(nve+is)/2. Non-negativity
 * a,b,c >= 0 becomes is>=0 (b>=0), nve>=is (a>=0), and nve+is<=2*orb (c>=0). A
 * running code violating any of these is not an admissible DRT row and the CSF
 * is rejected. The paper builds the DRT following Shavitt (Sec. 4.1) without
 * restating the ABC conditions.
 */
hexpr_drt_t *hexpr_drt_from_csfset(int norb,
                                    const unsigned char *csfs,
                                    int n_csfs)
{
    static const int DELTA_NVE[4] = {0, 1, 1, 2};
    static const int DELTA_IS [4] = {0, 1,-1, 0};

    if (!csfs) {
        hexpr_set_error("hexpr_drt_from_csfset: csfs is NULL");
        return NULL;
    }
    if (norb < 1) {
        hexpr_set_error("hexpr_drt_from_csfset: norb=%d < 1", norb);
        return NULL;
    }
    if (n_csfs < 1) {
        hexpr_set_error("hexpr_drt_from_csfset: n_csfs=%d < 1", n_csfs);
        return NULL;
    }

    htab_t *h_nodes = htab_new();

    /* Insert bottom node (0,0,0) */
    {
        node_t *n = calloc(1, sizeof(node_t));
        n->status = 1;
        n->number_upper_arcs = ARC_NIL;
        n->number_lower_arcs = ARC_NIL;
        for (int w = 0; w < 4; w++)
            n->upper_arc[w] = n->arc_weight[w] = n->lower_arc[w] = ARC_NIL;
        htab_insert(h_nodes, n);
    }

    int top_nve = -1, top_is = -1, top_set = 0;

    for (int ci = 0; ci < n_csfs; ci++) {
        const unsigned char *csf = csfs + (size_t)ci * (size_t)norb;
        int orb = 0, nve = 0, is = 0;

        for (int j = 0; j < norb; j++) {
            int step = (int)csf[j];
            if (step < 0 || step > 3) {
                hexpr_set_error(
                    "hexpr_drt_from_csfset: CSF[%d] step[%d]=%d out of range 0..3",
                    ci, j, step);
                htab_free(h_nodes);
                return NULL;
            }
            nve += DELTA_NVE[step];
            is  += DELTA_IS [step];
            orb++;

            if (is < 0) {
                hexpr_set_error(
                    "hexpr_drt_from_csfset: CSF[%d] b<0: is=%d at orb=%d",
                    ci, is, orb);
                htab_free(h_nodes);
                return NULL;
            }
            if (nve < is) {
                hexpr_set_error(
                    "hexpr_drt_from_csfset: CSF[%d] a<0: nve=%d < is=%d at orb=%d",
                    ci, nve, is, orb);
                htab_free(h_nodes);
                return NULL;
            }
            if (nve + is > 2 * orb) {
                hexpr_set_error(
                    "hexpr_drt_from_csfset: CSF[%d] c<0: nve+is=%d > 2*orb=%d at orb=%d",
                    ci, nve + is, 2 * orb, orb);
                htab_free(h_nodes);
                return NULL;
            }

            if (!htab_lookup(h_nodes, orb, nve, is)) {
                node_t *n = calloc(1, sizeof(node_t));
                n->orb = orb; n->nve = nve; n->is = is;
                n->status = 1;
                n->number_upper_arcs = ARC_NIL;
                n->number_lower_arcs = ARC_NIL;
                for (int w = 0; w < 4; w++)
                    n->upper_arc[w] = n->arc_weight[w] = n->lower_arc[w] = ARC_NIL;
                htab_insert(h_nodes, n);
            }
        }

        if (!top_set) {
            top_nve = nve; top_is = is; top_set = 1;
        } else if (nve != top_nve || is != top_is) {
            hexpr_set_error(
                "hexpr_drt_from_csfset: CSF[%d] top=(%d,%d,%d) differs from CSF[0] top=(%d,%d,%d)",
                ci, norb, nve, is, norb, top_nve, top_is);
            htab_free(h_nodes);
            return NULL;
        }
    }

    /* Sort nodes lex by (orb, nve, is) */
    int ncount = h_nodes->order_size;
    node_t **tmp_ptrs = malloc(ncount * sizeof(node_t *));
    for (int i = 0; i < ncount; i++)
        tmp_ptrs[i] = h_nodes->order[i];
    qsort(tmp_ptrs, ncount, sizeof(node_t *), node_cmp);

    node_t *nodes = malloc(ncount * sizeof(node_t));
    for (int i = 0; i < ncount; i++)
        nodes[i] = *tmp_ptrs[i];
    free(tmp_ptrs);
    htab_free(h_nodes);   /* frees all node_t* owned by htab */

    /* Find top_index (always ncount-1 after lex sort, but scan to verify) */
    int top_index = -1;
    for (int i = ncount - 1; i >= 0; i--) {
        if (nodes[i].orb == norb) { top_index = i; break; }
    }
    if (top_index < 0) { free(nodes); return NULL; }   /* shouldn't happen */

    /* Build code2num and wire upper/lower arcs */
    htab_t *code2num_h = htab_new();
    for (int i = 0; i < ncount; i++) {
        node_t *proxy = calloc(1, sizeof(node_t));
        proxy->orb = nodes[i].orb;
        proxy->nve = nodes[i].nve;
        proxy->is  = nodes[i].is;
        proxy->number_upper_arcs = i;   /* store index in this field */
        proxy->status = 1;
        for (int w = 0; w < 4; w++)
            proxy->upper_arc[w] = proxy->arc_weight[w] = proxy->lower_arc[w] = ARC_NIL;
        proxy->number_lower_arcs = ARC_NIL;
        htab_insert(code2num_h, proxy);
    }

    for (int inode = 0; inode < ncount; inode++) {
        node_t *nd = &nodes[inode];
        int succ[4][3] = {
            {nd->orb+1, nd->nve,   nd->is  },
            {nd->orb+1, nd->nve+1, nd->is+1},
            {nd->orb+1, nd->nve+1, nd->is-1},
            {nd->orb+1, nd->nve+2, nd->is  }
        };
        for (int w = 0; w < 4; w++) {
            node_t *p = htab_lookup(code2num_h, succ[w][0], succ[w][1], succ[w][2]);
            int jnode = p ? p->number_upper_arcs : ARC_NIL;
            nd->upper_arc[w] = jnode;
            if (jnode != ARC_NIL)
                nodes[jnode].lower_arc[w] = inode;
        }
    }
    htab_free(code2num_h);

    /* Compute arc counts, arc weights, and path weights */
    store_number_upper_arcs(nodes, 0, top_index);
    store_number_lower_arcs(nodes, top_index);

    ivec_t *lower_pw = calloc(ncount, sizeof(ivec_t));
    ivec_push(&lower_pw[0], 0);
    tree_search_lower(nodes, 0, top_index, lower_pw, 0);

    ivec_t *upper_pw = calloc(ncount, sizeof(ivec_t));
    for (int i = 0; i < ncount; i++) {
        if (nodes[i].number_upper_arcs != ARC_NIL)
            tree_search_upper(nodes, i, top_index, 0, &upper_pw[i]);
    }

    /* Assemble result */
    hexpr_drt_t *drt = calloc(1, sizeof(hexpr_drt_t));
    drt->mo        = mo_init(0, norb, 0);  /* ncore=0, nval=norb, next=0 */
    drt->info.nve  = top_nve;
    drt->info.is   = top_is;
    drt->kind      = HEXPR_DRT_CSFSET;
    drt->nodes     = nodes;
    drt->nnodes    = ncount;
    drt->top_index = top_index;

    drt->lower_pw   = malloc(ncount * sizeof(int *));
    drt->lower_pw_n = malloc(ncount * sizeof(int));
    drt->upper_pw   = malloc(ncount * sizeof(int *));
    drt->upper_pw_n = malloc(ncount * sizeof(int));
    for (int i = 0; i < ncount; i++) {
        drt->lower_pw[i]   = lower_pw[i].data;
        drt->lower_pw_n[i] = lower_pw[i].size;
        drt->upper_pw[i]   = upper_pw[i].data;
        drt->upper_pw_n[i] = upper_pw[i].size;
    }
    free(lower_pw);
    free(upper_pw);

    /* Build packed CSF buffer (may contain more paths than n_csfs if nodes shared) */
    {
        static const unsigned char CHAR_TO_STEP[256] = {
            ['e'] = 0, ['u'] = 1, ['d'] = 2, ['f'] = 3
        };
        char **csf_list = NULL;
        int ncsf = hexpr_drt_mk_csf_list(drt, &csf_list);

        /* Construction-time invariant -- mirror of the one in
         * hexpr_drt_create: nodes[0] must be the walk origin (0,0,0). Every
         * supplied CSF is a step sequence accumulated from the vacuum, so
         * the union graph's bottom is always (0,0,0); nothing else can
         * legitimately occupy nodes[0]. */
        if (drt->nodes[0].orb != 0 || drt->nodes[0].nve != 0 || drt->nodes[0].is != 0) {
            hexpr_set_error(
                "hexpr_drt_from_csfset: internal invariant violated: "
                "nodes[0] is (orb=%d,nve=%d,is=%d), not the walk origin "
                "(0,0,0); refusing to build a DRT every accessor would misread",
                drt->nodes[0].orb, drt->nodes[0].nve, drt->nodes[0].is);
            for (int i = 0; i < ncsf; i++) free(csf_list[i]);
            free(csf_list);
            hexpr_drt_free(drt);
            return NULL;
        }

        /* Construction-time consistency guard -- mirror of the one in
         * hexpr_drt_create: nodes[0].number_upper_arcs (hexpr_drt_walk_ncsfs())
         * must equal the walk count mk_csf_list actually enumerated. This is
         * independent of the documented CSFSET superset/dedup relationship
         * between n_csfs (input) and ncsf (canonical) -- it catches a dead
         * node surviving prune_unreachable in the union graph, not a
         * legitimate superset. */
        if (ncsf != drt->nodes[0].number_upper_arcs) {
            hexpr_set_error(
                "hexpr_drt_from_csfset: canonical CSF count mismatch: "
                "mk_csf_list enumerated %d walk(s) but "
                "nodes[0].number_upper_arcs reports %d; refusing to build "
                "an internally inconsistent DRT",
                ncsf, drt->nodes[0].number_upper_arcs);
            for (int i = 0; i < ncsf; i++) free(csf_list[i]);
            free(csf_list);
            hexpr_drt_free(drt);
            return NULL;
        }

        drt->csf_packed   = malloc((size_t)ncsf * (size_t)norb);
        drt->ncsfs_cached = ncsf;
        for (int i = 0; i < ncsf; i++) {
            for (int j = 0; j < norb; j++)
                drt->csf_packed[i * norb + j] = CHAR_TO_STEP[(unsigned char)csf_list[i][j]];
            free(csf_list[i]);
        }
        free(csf_list);
    }

    /* Retain the user input set verbatim (input order) for the user-facing
     * accessors; the canonical csf_packed above may be a superset. */
    drt->csf_input    = malloc((size_t)n_csfs * (size_t)norb);
    drt->n_csfs_input = n_csfs;
    memcpy(drt->csf_input, csfs, (size_t)n_csfs * (size_t)norb);

    return drt;
}

/* ------------------------------------------------------------------ */
/* hexpr_drt_free                                                      */
/* ------------------------------------------------------------------ */

/**
 * @brief Free a DRT and all storage it owns (per-node path-weight lists, node
 *        array, packed/input CSF buffers).
 * @param drt table to free; NULL is a no-op.
 * @note Invalidates every pointer previously returned by the accessors
 *       (@ref hexpr_drt_lower_path_weights, @ref hexpr_drt_walk_csfs, ...).
 *       Pairs with @ref hexpr_drt_create / @ref hexpr_drt_from_csfset.
 */
void hexpr_drt_free(hexpr_drt_t *drt)
{
    if (!drt) return;
    for (int i = 0; i < drt->nnodes; i++) {
        free(drt->lower_pw[i]);
        free(drt->upper_pw[i]);
    }
    free(drt->lower_pw);
    free(drt->lower_pw_n);
    free(drt->upper_pw);
    free(drt->upper_pw_n);
    free(drt->nodes);
    free(drt->csf_packed);
    free(drt->csf_input);   /* NULL for SOCI/FOCI; no-op */
    free(drt);
}

/** @brief Public alias for @ref hexpr_drt_free -- the name declared in
 *  @c hexpr.h and the only one the bindings call.
 *  @param drt table to free. */
void hexpr_drt_destroy(hexpr_drt_t *drt) { hexpr_drt_free(drt); }

/* ------------------------------------------------------------------ */
/* Public accessors                                                    */
/* ------------------------------------------------------------------ */

/** @brief Total orbital count (core+valence+external). @param drt DRT. @return norb. */
int hexpr_drt_norb(const hexpr_drt_t *drt)   { return drt->mo.norb; }
/** @brief Number of nodes in the assembled DRT. @param drt DRT. @return node count. */
int hexpr_drt_nnodes(const hexpr_drt_t *drt)  { return drt->nnodes; }
/* User-facing count: input count for CSFSET, canonical otherwise.
 * Keeps the invariant len(hexpr_drt_csfs()) == hexpr_drt_ncsfs(). */
/**
 * @brief User-facing CSF count.
 * @param drt DRT.
 * @return for a CSFSET sub-DRT, the retained input count; otherwise the
 *         canonical walk count (the bottom node's @c number_upper_arcs).
 * @note Matches the row count of @ref hexpr_drt_csfs. For the canonical count of
 *       a CSFSET sub-DRT use @ref hexpr_drt_walk_ncsfs instead.
 */
int hexpr_drt_ncsfs(const hexpr_drt_t *drt)
{
    if (drt->kind == HEXPR_DRT_CSFSET) return drt->n_csfs_input;
    return drt->nodes[0].number_upper_arcs;
}
/** @brief DRT kind (SOCI / FOCI / CSFSET). @param drt DRT. @return the kind enum. */
hexpr_drt_kind_t hexpr_drt_kind(const hexpr_drt_t *drt) { return drt->kind; }

/** @brief Core orbital count.     @param drt DRT. @return ncore. */
int hexpr_drt_ncore(const hexpr_drt_t *drt)  { return drt->mo.ncore; }
/** @brief Valence orbital count.  @param drt DRT. @return nval. */
int hexpr_drt_nval(const hexpr_drt_t *drt)   { return drt->mo.nval; }
/** @brief External orbital count. @param drt DRT. @return next. */
int hexpr_drt_next(const hexpr_drt_t *drt)   { return drt->mo.next; }
/** @brief Target-state valence-electron count. @param drt DRT. @return nve. */
int hexpr_drt_nve_state(const hexpr_drt_t *drt) { return drt->info.nve; }
/** @brief Target-state spin times two. @param drt DRT. @return is (=2S). */
int hexpr_drt_spin_x2(const hexpr_drt_t *drt)   { return drt->info.is; }

/**
 * @brief Read a node's code @c (orb,nve,is) by index.
 * @param drt DRT.
 * @param idx node index in @c [0, nnodes).
 * @param[out] out_orb,out_nve,out_is receive the node code (untouched on error).
 * @return ::HEXPR_OK, or ::HEXPR_ERR_RANGE if @p idx is out of range.
 */
hexpr_status_t hexpr_drt_node_info(const hexpr_drt_t *drt, int idx,
                                   int *out_orb, int *out_nve, int *out_is)
{
    if (idx < 0 || idx >= drt->nnodes) return HEXPR_ERR_RANGE;
    *out_orb = drt->nodes[idx].orb;
    *out_nve = drt->nodes[idx].nve;
    *out_is  = drt->nodes[idx].is;
    return HEXPR_OK;
}

/* User-facing CSF list: the retained input set (input order) for CSFSET
 * sub-DRTs; the canonical buffer for SOCI/FOCI (where canonical == the
 * user's parametric result, no input set exists). */
/**
 * @brief User-facing packed CSF buffer (step codes 0..3, row-major @c norb cols).
 * @param drt DRT.
 * @return for CSFSET, the retained input rows in caller order; otherwise the
 *         canonical packed buffer. Borrowed -- owned by @p drt, valid until
 *         @ref hexpr_drt_destroy.
 * @note Row count is @ref hexpr_drt_ncsfs. See @ref hexpr_drt_walk_csfs for the
 *       always-canonical view.
 */
const unsigned char *hexpr_drt_csfs(const hexpr_drt_t *drt)
{
    if (drt->kind == HEXPR_DRT_CSFSET) return drt->csf_input;
    return drt->csf_packed;
}

/* Canonical GUGA-walk CSF list and count for ALL kinds. For a CSFSET
 * sub-DRT this is the internal (possibly superset) view; for SOCI/FOCI it
 * equals hexpr_drt_csfs / hexpr_drt_ncsfs. */
/**
 * @brief Canonical GUGA-walk packed CSF buffer (all kinds).
 * @param drt DRT.
 * @return the canonical packed buffer (step codes 0..3, @c norb cols). Borrowed;
 *         valid until @ref hexpr_drt_destroy. For a CSFSET sub-DRT this may
 *         be a SUPERSET of the input; for SOCI/FOCI it equals
 *         @ref hexpr_drt_csfs.
 * @note Row count is @ref hexpr_drt_walk_ncsfs.
 */
const unsigned char *hexpr_drt_walk_csfs(const hexpr_drt_t *drt)
{
    return drt->csf_packed;
}

/**
 * @brief Canonical GUGA-walk CSF count (all kinds).
 * @param drt DRT.
 * @return the number of canonical walks (the bottom node's @c number_upper_arcs);
 *         the row count of @ref hexpr_drt_walk_csfs.
 */
int hexpr_drt_walk_ncsfs(const hexpr_drt_t *drt)
{
    return drt->nodes[0].number_upper_arcs;
}

/** @brief Orbital of node @p i.  @param drt DRT. @param i node index. @return orb. */
int hexpr_drt_orb(const hexpr_drt_t *drt, int i)  { return drt->nodes[i].orb; }
/** @brief Electron count of node @p i. @param drt DRT. @param i node index. @return nve. */
int hexpr_drt_nve(const hexpr_drt_t *drt, int i)  { return drt->nodes[i].nve; }
/** @brief Spin label of node @p i as @c is+1 (1-based). @param drt DRT. @param i node index. @return is+1. */
int hexpr_drt_spin(const hexpr_drt_t *drt, int i) { return drt->nodes[i].is + 1; }

/** @brief Upper-arc target index for node @p i, step @p w (or @ref ARC_NIL).
 *  @param drt DRT. @param i node index. @param w step 0..3. @return node index or ARC_NIL. */
int hexpr_drt_upper_arc(const hexpr_drt_t *drt, int i, int w)
{ return drt->nodes[i].upper_arc[w]; }

/** @brief Lower-arc target index for node @p i, step @p w (or @ref ARC_NIL).
 *  @param drt DRT. @param i node index. @param w step 0..3. @return node index or ARC_NIL. */
int hexpr_drt_lower_arc(const hexpr_drt_t *drt, int i, int w)
{ return drt->nodes[i].lower_arc[w]; }

/** @brief Weight of node @p i's upper arc @p w (or @ref ARC_NIL).
 *  @param drt DRT. @param i node index. @param w step 0..3. @return arc weight or ARC_NIL.
 *  @note Set by @ref store_number_upper_arcs; see that routine for the weight
 *  convention (deferred to its domain marker). */
int hexpr_drt_arc_weight(const hexpr_drt_t *drt, int i, int w)
{ return drt->nodes[i].arc_weight[w]; }

/** @brief Number of upper arcs (upward walks) from node @p i.
 *  @param drt DRT. @param i node index. @return number_upper_arcs. */
int hexpr_drt_number_upper_arcs(const hexpr_drt_t *drt, int i)
{ return drt->nodes[i].number_upper_arcs; }

/** @brief Number of lower arcs (downward walks) from node @p i.
 *  @param drt DRT. @param i node index. @return number_lower_arcs. */
int hexpr_drt_number_lower_arcs(const hexpr_drt_t *drt, int i)
{ return drt->nodes[i].number_lower_arcs; }

/**
 * @brief Borrow node @p i's lower path-weight list.
 * @param drt DRT.
 * @param i   node index.
 * @param[out] out_n receives the list length (0 if the node has none).
 * @return pointer to the internal int array (length @c *out_n), or possibly NULL
 *         when empty. Borrowed -- owned by @p drt, valid until @ref hexpr_drt_free;
 *         do not free. Cross-link hub read by @c expr.c.
 */
const int *hexpr_drt_lower_path_weights(const hexpr_drt_t *drt, int i, int *out_n)
{
    *out_n = drt->lower_pw_n[i];
    return drt->lower_pw[i];
}

/**
 * @brief Borrow node @p i's upper path-weight list.
 * @param drt DRT.
 * @param i   node index.
 * @param[out] out_n receives the list length (0 if the node has none).
 * @return pointer to the internal int array (length @c *out_n), or possibly NULL
 *         when empty. Borrowed -- owned by @p drt, valid until @ref hexpr_drt_free;
 *         do not free. Cross-link hub read by @c expr.c.
 */
const int *hexpr_drt_upper_path_weights(const hexpr_drt_t *drt, int i, int *out_n)
{
    *out_n = drt->upper_pw_n[i];
    return drt->upper_pw[i];
}

/* ------------------------------------------------------------------ */
/* hexpr_drt_mk_csf_list                                               */
/* ------------------------------------------------------------------ */

/**
 * @brief Enumerate the DRT's canonical CSFs as an array of step-strings.
 * @param drt DRT to walk.
 * @param[out] out_steps receives a freshly @c malloc'd array of @c norb+1 byte
 *             strings (chars from {'e','u','d','f'}); the CALLER owns the array
 *             and every string and must free them.
 * @return the number of CSFs written (array length).
 *
 * Role: the public driver over @ref tree_search_csf. Sets up a one-shot
 * enumeration context with a reusable @c norb+1 work buffer, runs the DFS from
 * the bottom node to @c top_index collecting one string per complete walk, frees
 * the work buffer, and hands the accumulated list out. Used during construction
 * to build the packed canonical buffer and by @c expr.c to label its output.
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * The DFS enumeration order is the canonical CSF order because it visits upper
 * arcs in the same lexical order the arc weights encode (Shavitt, refs [4,5]):
 * the i-th walk it produces is the walk with canonical index i, so walk-DFS order
 * matches the path-weight address from tree_search_lower/_upper. This is the
 * ordering csf_lookup.c and expr.c index against. The paper references the
 * canonical DRT ordering (Sec. 4.1) without restating it.
 */
int hexpr_drt_mk_csf_list(const hexpr_drt_t *drt, char ***out_steps)
{
    csf_ctx_t ctx;
    ctx.list    = NULL;
    ctx.size    = 0;
    ctx.cap     = 0;
    ctx.norb    = drt->mo.norb;
    ctx.depth   = 0;
    ctx.csf_buf = malloc(ctx.norb + 1);

    tree_search_csf(drt->nodes, 0, drt->top_index, &ctx);
    free(ctx.csf_buf);
    *out_steps = ctx.list;
    return ctx.size;
}

/* ------------------------------------------------------------------ */
/* hexpr_drt_show                                                      */
/* ------------------------------------------------------------------ */

/**
 * @brief Print the DRT as the reference's formatted tables (nodes + upper arcs +
 *        weights, then lower arcs) to a stream.
 * @param drt DRT to print (must be SOCI or FOCI).
 * @param f   destination stream.
 * @return ::HEXPR_OK, or ::HEXPR_ERR_INVALID_ARG for a CSFSET sub-DRT, which
 *         has no meaningful tabular form.
 *
 * Role: a human-readable dump mirroring the original Fortran reference layout
 * (title, STATE/MO summary, counts, then two tables). Both tables iterate nodes
 * in DESCENDING index order (top first), and every @ref ARC_NIL field is printed
 * as a dash. Pure formatting -- it only reads node fields filled in during
 * construction.
 *
 * @note No GUGA intent here; the column meanings (orb/nve/is, arcs, weights) are
 *       defined by the construction routines that populate them.
 */
hexpr_status_t hexpr_drt_show(const hexpr_drt_t *drt, FILE *f)
{
    /* hexpr_drt_show is not meaningful for sub-DRTs built from CSF sets */
    if (drt->kind == HEXPR_DRT_CSFSET)
        return HEXPR_ERR_INVALID_ARG;

    const mo_t   *mo   = &drt->mo;
    const info_t *info = &drt->info;
    int n = drt->nnodes;

    /* title */
    if (drt->kind == HEXPR_DRT_SOCI)
        fprintf(f, " \n  D I S T I N C T  R O W  T A B L E --- S O C I                                      coded by T. Noro \n \n");
    else
        fprintf(f, " \n  D I S T I N C T  R O W  T A B L E --- F O C I                                      coded by T. Noro \n \n");

    fprintf(f, "   STATE information   NVE = %3d,  S*2 = %3d \n", info->nve, info->is);
    fprintf(f, "   MO classification   NORB = %3d   ( NCORE=%3d, NVAL=%3d, NEXT=%3d )\n\n",
            mo->norb, mo->ncore, mo->nval, mo->next);
    fprintf(f, "   Number of generated nodes = %7d \n", n);
    fprintf(f, "   Number of generated csfs  = %7d \n", hexpr_drt_ncsfs(drt));

    fprintf(f, "\n *** Nodes, Upper Arc Table, and Arc Weights\n");
    fprintf(f, "                             upper arcs                #arcs         upper weight\n");
    fprintf(f, "     n    mo  el  is         e      u      d      f                  e        u        d        f\n");

    for (int node = n - 1; node >= 0; node--) {
        const node_t *nd = &drt->nodes[node];
        fprintf(f, "%6d ( %3d %3d %3d ) ", node, nd->orb, nd->nve, nd->is);

        /* upper arcs */
        for (int w = 0; w < 4; w++) {
            if (nd->upper_arc[w] == ARC_NIL)
                fprintf(f, "      -");
            else
                fprintf(f, "  %5d", nd->upper_arc[w]);
        }
        /* number_upper_arcs */
        if (nd->number_upper_arcs == ARC_NIL)
            fprintf(f, "        -");
        else
            fprintf(f, "  %7d ", nd->number_upper_arcs);
        /* arc weights */
        for (int w = 0; w < 4; w++) {
            if (nd->arc_weight[w] == ARC_NIL)
                fprintf(f, "        -");
            else
                fprintf(f, "  %7d", nd->arc_weight[w]);
        }
        fprintf(f, "\n");
    }

    fprintf(f, "\n *** Lower Arc Table\n");
    fprintf(f, "            lower arcs\n");
    fprintf(f, "    n       e      u      d      f\n");
    for (int node = n - 1; node >= 0; node--) {
        const node_t *nd = &drt->nodes[node];
        fprintf(f, "%5d ", node);
        for (int w = 0; w < 4; w++) {
            if (nd->lower_arc[w] == ARC_NIL)
                fprintf(f, "      -");
            else
                fprintf(f, "  %5d", nd->lower_arc[w]);
        }
        if (nd->number_lower_arcs == ARC_NIL)
            fprintf(f, "        -");
        else
            fprintf(f, "  %7d ", nd->number_lower_arcs);
        fprintf(f, "\n");
    }
    return HEXPR_OK;
}

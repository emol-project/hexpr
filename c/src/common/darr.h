/* common/darr.h -- generic growable array of void* pointers */
#ifndef HEXPR_DARR_H_
#define HEXPR_DARR_H_

typedef struct {
    void **data;
    int    size;
    int    cap;
} darr_t;

void  darr_init (darr_t *a);
void  darr_push (darr_t *a, void *item);
void *darr_get  (const darr_t *a, int i);
int   darr_size (const darr_t *a);
void  darr_clear(darr_t *a);           /* reset size to 0, keep allocation */
void  darr_free (darr_t *a);           /* free backing store */

#endif /* HEXPR_DARR_H_ */

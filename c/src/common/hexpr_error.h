#ifndef HEXPR_ERROR_H_
#define HEXPR_ERROR_H_
/* Internal: set the thread-local error string returned by hexpr_last_error(). */
void hexpr_set_error(const char *fmt, ...);
#endif /* HEXPR_ERROR_H_ */

#ifndef AMR_UTIL_H
#define AMR_UTIL_H

#include <p4est_to_p8est.h>
#include <p8est.h>

/**
 * Get the coordinates of the midpoint of a quadrant
 */
void get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3]);

#endif //AMR_UTIL_H

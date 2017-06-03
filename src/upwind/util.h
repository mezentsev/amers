#ifndef AMR_UTIL_H
#define AMR_UTIL_H

#include <p4est_to_p8est.h>
#include <p8est.h>

int ipow(int base, int exp);

/**
 * Get the coordinates of the midpoint of a quadrant
 */
void get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3]);

/**
 * Detect boundary
 * @param p4est
 * @param tt
 * @param node
 * @return
 */
int is_quadrant_on_face_boundary (p8est_t * p4est,
                                  p4est_topidx_t treeid,
                                  int face,
                                  const p8est_quadrant_t * q);

/**
 * Вывод информации о ячейке
 * @param q
 * @param is_ghost
 */
void quadrant_pprint (p8est_quadrant_t * q, int is_ghost);

#endif //AMR_UTIL_H

#ifndef AMR_UTIL_H
#define AMR_UTIL_H

#include <p4est_to_p8est.h>
#include <p8est.h>

#define TO_LEFT 0
#define TO_RIGHT 1
#define TO_BOTTOM 2
#define TO_TOP 3
#define TO_FRONT 4
#define TO_BEHIND 5

#define FROM_RIGHT 0
#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_BOTTOM 3
#define FROM_BEHIND 4
#define FROM_FRONT 5


int ipow(int base, int exp);

/**
 * Получение центра ячейки
 *
 * @param p8est
 * @param tree
 * @param q
 * @param xyz
 */
void get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3]);

/**
 * Проверка на границу
 *
 * @param p4est
 * @param tt
 * @param node
 * @return
 */
int is_quadrant_on_face_boundary (p8est_t *p4est,
                                  p4est_topidx_t treeid,
                                  int face,
                                  const p8est_quadrant_t * q);

/**
 * Вывод информации о ячейке
 *
 * @param q
 * @param is_ghost
 */
void quadrant_pprint (p8est_quadrant_t * q, int is_ghost);

/**
 * Декодирует face относительно следующего соседа
 *
 * @param next_face следующий face
 * @param next_subface следующий subface
 * @return
 */
int get_neighbour_face_by_next_one(const int next_face, const int next_subface);

#endif //AMR_UTIL_H

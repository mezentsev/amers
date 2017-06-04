#ifndef AMR_SOLVER_H
#define AMR_SOLVER_H

#include "data.h"
#include "upwind.h"

/**
 * Решить шаг N
 *
 * @param p8est
 * @param step шаг
 */
void solver_step(p8est_t *p8est,
                 int step);

/**
 * Возвращает данные с границы
 *
 *        e-------f
 *       /|      /|
 *      / |     / |
 *     a--|----b  |
 *     |  g----|--h
 *     | /     | /
 *     c-------d
 *
 *           3   5
 *           +  +
 *           | /
 *           |/
 *  1 - -----.---- + 0
 *          /|
 *         / |
 *        -  -
 *       4   2
 *
 * @param quad граничная ячейка
 * @param face номер относительно границы
 * @return заполненные данные
 */
element_data_t get_boundary_data_by_face(p8est_t *p8est, p8est_quadrant_t *quad, int face);

/**
 * Обход по сетке для вычисления потока
 *
 * @param mesh сетка
 */
void calc_flux_mesh_iter(p8est_t *p8est,
                         p8est_mesh_t *mesh,
                         p8est_ghost_t *ghost,
                         void *ghost_data);


/** Approximate the flux across a boundary between quadrants.
 *
 * We use a very simple upwind numerical flux.
 *
 * This function matches the p4est_iter_face_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info the information about the quadrants on either side of the
 *                  interface, populated by p4est_iterate()
 * \param [in] user_data the user_data given to p4est_iterate(): in this case,
 *                       it points to the ghost_data array, which contains the
 *                       element_data_t data for all of the ghost cells, which
 *                       was populated by p4est_ghost_exchange_data()
 */
void calc_flux_face_iter(p8est_iter_face_info_t *face_info, void *user_data);

/**
* Сброс производных
*
* @param p8est
* @param mesh сетка
*/
void reset_derivatives(p8est_t *p8est,
                       p8est_mesh_t *mesh);


/**
 * Approximate the divergence of (vu) on each quadrant
 *
 * We use piecewise constant approximations on each quadrant, so the value is
 * always 0.
 */
void quad_divergence (p8est_t *p8est,
                      p8est_mesh_t *mesh);


/** For two quadrants on either side of a face, estimate the derivative normal
 * to the face.
 *
 * This function matches the p4est_iter_face_t prototype used by
 * p4est_iterate().
 *
 * \param [in] info          the information about this quadrant that has been
 *                           populated by p4est_iterate()
 * \param [in] user_data the user_data given to p4est_iterate(): in this case,
 *                       it points to the ghost_data array, which contains the
 *                       elements_data_t data for all of the ghost cells, which
 *                       was populated by p4est_ghost_exchange_data()
 */
void minmod_estimate (p8est_iter_face_info_t * info, void *user_data);


void timestep_update_volume_iter(p8est_iter_volume_info_t *info,
                                 void *user_data);


/**
 * Получить новое значение dt
 *
 * @param p8est
 * @param mesh
 * @return dt
 */
double timestep_update_mesh_iter(p8est_t *p8est,
                                 p8est_mesh_t *mesh);

/**
 * Инициализация солвера
 *
 * @param p8est
 * @param which_tree
 * @param quad
 */
void init_solver(p8est_t *p8est,
                 p4est_topidx_t which_tree,
                 p8est_quadrant_t *quad);

/**
 * Созранить решение в VTK
 * @param p8est
 * @param step
 */
void write_solution(p8est_t *p8est,
                    int step);


#endif //AMR_SOLVER_H

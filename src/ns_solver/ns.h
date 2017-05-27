//
// Created by Vadim Mezentsev on 01.01.17.
//

#ifndef AMR_NS_H
#define AMR_NS_H

#include <p4est_to_p8est.h>
#include <p8est_connectivity.h>
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <p8est_iterate.h>
#include <time.h>
#include <p8est_iterate.h>
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_mesh.h>
#include <p8est_search.h>
#include <p8est_ghost.h>

#include "data.h"
#include "solver.h"

/**
 * Созранить решение в VTK
 * @param p8est
 * @param step
 */
void write_solution(p8est_t *p8est,
                    int step);

/**
 * Решить шаг N
 *
 * @param p8est
 * @param step
 */
void solver_step(p8est_t *p8est,
                 int step);

/**
 * Обход всех ячеек для вычисления временного шага
 *
 * @param p8est_iter_volume_info_t
 * @param user_data
 */
void calc_cfl_timestep(p8est_iter_volume_info_t *info,
                       void *user_data);

/**
 * Обход по всем соседям
 *
 * @param p8est
 * @param ghost
 * @param mesh
 * @param ghost_data
 */
void mesh_neighbors_iter(p8est_t *p8est,
                         p8est_ghost_t *ghost,
                         p8est_mesh_t *mesh,
                         void *ghost_data);

/**
 * Функция безусловного деления ячейки
 *
 * @param p8est
 * @param which_tree
 * @param q
 * @return
 */
int refine_always (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);

/**
 * Деление ячейки с приведением к геометрии
 *
 * @param p8est
 * @param which_tree
 * @param q
 * @return
 */
int refine_fn (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);

/**
 * Запись решения в VTK
 *
 * @param p8est
 * @param step
 */
void write_vtk(p8est_t *p8est, int step);

#endif //AMR_NS_H

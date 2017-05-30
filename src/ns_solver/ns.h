//
// Created by Vadim Mezentsev on 01.01.17.
//

#ifndef AMR_NS_H
#define AMR_NS_H

#include <p4est_to_p8est.h>
#include <p8est.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <p8est_iterate.h>
#include <p8est_bits.h>
#include <p8est_mesh.h>
#include <p8est_search.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_bits.h>
#include <p8est_algorithms.h>
#include <p8est_balance.h>
#include <p8est_communication.h>
#include <p8est_geometry.h>
#include <p8est_io.h>
#include <p8est_lnodes.h>
#include <p8est_wrap.h>
#include <p8est_points.h>
#include <p8est_nodes.h>

#include <time.h>


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

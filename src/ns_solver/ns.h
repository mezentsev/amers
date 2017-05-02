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
 * Write solution in step to VTK
 * @param p8est
 * @param step
 */
void write_solution(p8est_t *p8est,
                    int step);

/**
 * Run Solver
 * @param p8est
 * @param step
 */
void solve(p8est_t *p8est,
           int step);

/**
 * Iterate throw each cell in mesh
 * @param p8est
 * @param mesh
 */
void mesh_iter(p8est_t *p8est,
               p8est_mesh_t *mesh);

/**
 * Iterate throw each neighbor's cell
 * with their ghost data
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
 * Refine all cells
 * @param p8est
 * @param which_tree
 * @param q
 * @return
 */
int refine_always (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);

/**
 * Refine conditional function
 * @param p8est
 * @param which_tree
 * @param q
 * @return
 */
int refine_fn (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);

void write_vtk(p8est_t *p8est, int step);

#endif //AMR_NS_H

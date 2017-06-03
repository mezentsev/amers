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

#define FROM_RIGHT 0
#define FROM_LEFT 1
#define FROM_TOP 2
#define FROM_BOTTOM 3
#define FROM_BEHIND 4
#define FROM_FRONT 5

/**
 * Созранить решение в VTK
 * @param p8est
 * @param step
 */
void write_solution(p8est_t *p8est,
                    int step);

#endif //AMR_NS_H

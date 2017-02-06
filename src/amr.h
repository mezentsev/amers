//
// Created by Vadim Mezentsev on 01.01.17.
//

#ifndef AMR_AMR_H
#define AMR_AMR_H

#include <p4est_to_p8est.h>
#include <p8est_connectivity.h>
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <p8est_iterate.h>


/**
 * Data set for each cell
 */
typedef struct data
{
    double              f1; /* first function value */
    double              f2; /* second function value */
    double              v;
    double              e;  /* error estinate - fabs(f1-f2) */
} data_t;

/**
 * Solver function for 3 params
 * @param x
 * @param y
 * @param z
 * @return
 */
typedef double (*t_func_3)(double, double, double);
typedef double (*t_func_5)(double, double, double, double, double);

/**
 * Describes current problem
 */
typedef struct ctx
{
    double center[P4EST_DIM];  // coordinates of the center
    double width;

    t_func_3 f;
} ctx_t;

/**
 * Solver function
 * @param x
 * @param y
 * @param z
 * @return
 */
double s_func(double x, double y, double z)
{
    return pow(x, 2.) + pow(y, 2.) + pow(z, 2.);
}

double grad(t_func_3, double, double, double, double, double, double);
double laplacian(t_func_3, double, double, double);
void get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3]);

void init(p8est_t *p4est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
int refine_always (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
int refine_fn (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *q);
int coarsen_fn (p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t **children);
void volume_iter(p8est_iter_volume_info_t *info, void *user_data);
void face_iter_f1(p8est_iter_face_info_t *info, void *user_data);
void calc_error_iter(p8est_iter_volume_info_t *info, void *user_data);

void get_solution(p8est_iter_volume_info_t *info, void *user_data, int f);
void get_solution_f1(p8est_iter_volume_info_t *info, void *user_data);
void get_solution_f2(p8est_iter_volume_info_t *info, void *user_data);
void get_error_estimate(p8est_iter_volume_info_t *info, void *user_data);

void write_solution(p8est_t *p8est, int step);
void solve(p8est_t *p8est, int step);

/**
 * Debug print quadrant
 * @param q
 * @param rank
 */
void quadrant_print(p8est_quadrant_t *, int);

#endif //AMR_AMR_H

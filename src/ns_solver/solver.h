#ifndef AMR_SOLVER_H
#define AMR_SOLVER_H

#include "data.h"

/**
 * Initialize solver
 * @param q
 * @param ctx
 */
void init_solver(data_t *q, context_t *ctx);

/**
 * Calculate CFL for quadrant with length h
 * @param data
 * @param ctx
 * @param h
 */
void cflq(data_t *data, context_t *ctx, double h);

/**
 * Calculate flow
 * @param data
 * @param nx
 * @param ny
 * @param nz
 * @return
 */
double flow(data_t *data, double nx, double ny, double nz);

/**
 * Transform Z to Q
 * @param data
 */
void setq(data_t *data);

#endif //AMR_SOLVER_H

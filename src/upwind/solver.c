#include <math.h>
#include <sc.h>
#include "solver.h"
#include "../util.h"

void init_solver(p8est_t *p8est, element_data_t *data) {
    double          midpoint[3];
    context_t       *ctx = (context_t *) p8est->user_pointer;
    double          r2, dx, dy, dz;
    double          arg, retval;

    get_midpoint(p8est, data->quad.p.which_tree, &data->quad, midpoint);

    data->dux = -1;
    data->duy = -1;
    data->duz = -1;

    // Compute the value and derivatives of the initial condition.
    r2 = 0.;

    // расстояние до центра
    dx = midpoint[0] - ctx->center[0];
    dy = midpoint[1] - ctx->center[1];
    dz = midpoint[2] - ctx->center[2];

    // радиус
    r2 += dx * dx + dy * dy + dz * dz;

    data->u = exp(-(1./2.) * r2 / (ctx->width * ctx->width));

    data->dux = -(1./ (ctx->width * ctx->width)) * dx * data->u;
    data->duy = -(1./ (ctx->width * ctx->width)) * dy * data->u;
    data->duz = -(1./ (ctx->width * ctx->width)) * dz * data->u;

    data->dummy = 0;
}

void init_empty_solver(p8est_t *p8est, element_data_t *data) {
    data->dux = 0;
    data->duy = 0;
    data->duz = 0;

    data->dummy = 0;
}

element_data_t get_boundary_data_by_face(p8est_t *p8est,
                                         p8est_quadrant_t *q,
                                         int face) {
    double e;
    context_t *ctx = (context_t *) p8est->user_pointer;
    element_data_t *data = (element_data_t *) q->p.user_data;
    element_data_t boundary_data;

    init_empty_solver(p8est, &boundary_data);
    boundary_data.dummy = 1;

    /* координаты ячейки сейчас не учитываются */

    /* направление face от соседа к текущей ячейке */
    switch (face) {
        case 0:                      /* -x side */
            SC_PRODUCTION("-x side\n");
            // втекает
            boundary_data.dux        = data->dux;
            boundary_data.duy        = data->duy;
            boundary_data.duz        = data->duz;
            break;
        case 1:                      /* +x side */
            // вытекает
            SC_PRODUCTION("+x side\n");
            boundary_data.dux        = data->dux;
            boundary_data.duy        = data->duy;
            boundary_data.duz        = data->duz;
            break;
        case 2:                      /* -y side */
            SC_PRODUCTION("-y side\n");
            // стенка
            boundary_data.dux        = data->dux;
            boundary_data.duy        = -data->duy;
            boundary_data.duz        = data->duz;
            break;
        case 3:                      /* +y side */
            SC_PRODUCTION("+y side\n");
            // стенка
            boundary_data.dux        = data->dux;
            boundary_data.duy        = -data->duy;
            boundary_data.duz        = data->duz;
            break;
        case 4:                      /* -z side */
            SC_PRODUCTION("-z side\n");
            // стенка
            boundary_data.dux        = data->dux;
            boundary_data.duy        = data->duy;
            boundary_data.duz        = -data->duz;
            break;
        case 5:                      /* +z side */
            SC_PRODUCTION("+z side\n");
            // стенка
            boundary_data.dux        = data->dux;
            boundary_data.duy        = data->duy;
            boundary_data.duz        = -data->duz;
            break;
        default:
            SC_ABORT("Wrong face");
    }

    return boundary_data;
}

element_data_t calc_flux(p8est_t *p8est, p8est_quadrant_t *cur_quad, p8est_quadrant_t *n_quad, int nface) {
    element_data_t q_new;
    init_solver(p8est, &q_new);

    int nx = 0;
    int ny = 0;
    int nz = 0;

    /* направление face от соседа к текущей ячейке */
    switch (nface) {
        case 0:                      /* -x side */
            SC_PRODUCTION("-x side\n");
            nx = 1;
            break;
        case 1:                      /* +x side */
            SC_PRODUCTION("+x side\n");
            nx = -1;
            break;
        case 2:                      /* -y side */
            SC_PRODUCTION("-y side\n");
            ny = 1;
            break;
        case 3:                      /* +y side */
            SC_PRODUCTION("+y side\n");
            ny = -1;
            break;
        case 4:                      /* -z side */
            SC_PRODUCTION("-z side\n");
            nz = 1;
            break;
        case 5:                      /* +z side */
            SC_PRODUCTION("+z side\n");
            nz = -1;
            break;
        default:
            SC_ABORT("Wrong face");
    }

    return q_new;
}

void solver_step(p8est_t *p8est,
                 int step) {
    SC_PRODUCTIONF("Start solve step %d\n", step);

    int                 mpiret;
    element_data_t      *ghost_data;
    p8est_ghost_t       *ghost;
    p8est_mesh_t        *mesh;
    context_t           *ctx = (context_t *) p8est->user_pointer;

    /* выделение гостового слоя */
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC (element_data_t, ghost->ghosts.elem_count);
    mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FACE);

    /* exchange ghost data */
    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    /* calc timestep */
    //SC_PRODUCTION("Cell iter started\n");
    //p8est_iterate(p8est, NULL, NULL,        /* слой гостовых ячеек не нужен, доп. параметры тоже */
    //              calc_cfl_timestep_iter,        /* вычисление временного шага, проходя по всем ячейкам */
    //              NULL, NULL, NULL);
    //SC_PRODUCTION("Cell iter ended\n");

    /* calc min dt */
    //SC_PRODUCTIONF("dt old: %.20lf\n", ctx->dt);
    //mpiret = sc_MPI_Allreduce(MPI_IN_PLACE, &ctx->dt, 1, MPI_DOUBLE, MPI_MIN, p8est->mpicomm);
    //SC_CHECK_MPI(mpiret);
    //SC_PRODUCTIONF("dt new: %.20lf\n", ctx->dt);

    /* exchange ghost data */
    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    SC_PRODUCTION("Neighbors iter started\n");
    //p8est_iterate(p8est, ghost, (void *) ghost_data,        /* вкладываем гостовый слой */
    //              NULL,//calc_flux_volume_iter,
    //              calc_flux_face_iter,                      /* обход по фейсам каждой ячейки для высчитывания потока */
    //              NULL,
    //              NULL);
    calc_flux_mesh_iter(p8est, mesh, ghost, ghost_data);
    SC_PRODUCTION("Neighbors iter ended\n");

    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    /* generate vtk and print solution */
    SC_PRODUCTION("Write solution started\n");
    write_solution(p8est, step);
    SC_PRODUCTION("Write solution ended\n");

    /* очистка всех выделенных данных */
    P4EST_FREE (ghost_data);
    p8est_ghost_destroy(ghost);
    p8est_mesh_destroy(mesh);

    SC_PRODUCTIONF("End solve step %d\n", step);
}

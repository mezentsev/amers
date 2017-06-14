#include <math.h>
#include <sc.h>
#include "solver.h"
#include "../util.h"

void init_solver(p8est_t *p8est, element_data_t *data) {
    context_t *ctx = (context_t *) p8est->user_pointer;
    double e;

    data->Density = 1;
    data->Pressure = 1;
    data->u1 = 1.5;
    data->u2 = 0.1;
    data->u3 = 0;
    data->dummy = 0;

    e = data->Pressure / ((ctx->Adiabatic - 1) * data->Density);
    data->E = e + (pow(data->u1, 2) + pow(data->u2, 2) + pow(data->u3, 2))/2;

    updateQ(p8est, data);
}

void init_empty_solver(p8est_t *p8est, element_data_t *data) {
    data->Density = 0;
    data->Pressure = 0;
    data->u1 = 0;
    data->u2 = 0;
    data->u3 = 0;
    data->dummy = 0;

    data->E = 0;

    updateQ(p8est, data);
}

void updateQ(p8est_t *p8est, element_data_t *data) {
    data->Q.D   = data->Density;
    data->Q.Du1 = data->Density * data->u1;
    data->Q.Du2 = data->Density * data->u2;
    data->Q.Du3 = data->Density * data->u3;

    data->Q.PE  = data->Pressure * data->E;
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
            SC_LDEBUG("-x side\n");
            // втекает
            boundary_data.Density   = data->Density;
            boundary_data.Pressure  = data->Pressure;
            boundary_data.u1        = data->u1;
            boundary_data.u2        = data->u2;
            boundary_data.u3        = data->u3;
            break;
        case 1:                      /* +x side */
            // вытекает
            SC_LDEBUG("+x side\n");
            boundary_data.Density   = data->Density;
            boundary_data.Pressure  = data->Pressure;
            boundary_data.u1        = data->u1;
            boundary_data.u2        = data->u2;
            boundary_data.u3        = data->u3;
            break;
        case 2:                      /* -y side */
            SC_LDEBUG("-y side\n");
            // стенка
            boundary_data.Density   = data->Density;
            boundary_data.Pressure  = data->Pressure;
            boundary_data.u1        = data->u1;
            boundary_data.u2        = -data->u2;
            boundary_data.u3        = data->u3;
            break;
        case 3:                      /* +y side */
            SC_LDEBUG("+y side\n");
            // стенка
            boundary_data.Density   = data->Density;
            boundary_data.Pressure  = data->Pressure;
            boundary_data.u1        = data->u1;
            boundary_data.u2        = -data->u2;
            boundary_data.u3        = data->u3;
            break;
        case 4:                      /* -z side */
            SC_LDEBUG("-z side\n");
            // стенка
            boundary_data.Density   = data->Density;
            boundary_data.Pressure  = data->Pressure;
            boundary_data.u1        = data->u1;
            boundary_data.u2        = data->u2;
            boundary_data.u3        = -data->u3;
            break;
        case 5:                      /* +z side */
            SC_LDEBUG("+z side\n");
            // стенка
            boundary_data.Density   = data->Density;
            boundary_data.Pressure  = data->Pressure;
            boundary_data.u1        = data->u1;
            boundary_data.u2        = data->u2;
            boundary_data.u3        = -data->u3;
            break;
        default:
            SC_ABORT("Wrong face");
    }

    // нужно пересчитать энергию по новым данным
    e = boundary_data.Pressure / ((ctx->Adiabatic - 1) * boundary_data.Density);
    boundary_data.E = e + (pow(boundary_data.u1, 2) + pow(boundary_data.u2, 2) + pow(boundary_data.u3, 2))/2;

    // Z -> Q
    updateQ(p8est, &boundary_data);

    return boundary_data;
}

/**
 * Условие устойчивости (CFL - Куранта-Фридриха-Леви)
 *
 * @param data данные ячейки
 * @param ctx контекст. в dt записывается курантовый шаг
 * @param length сторона ячейки
 */
void cflq(element_data_t *data, context_t *ctx, double length) {
    P4EST_ASSERT(length > 0);

    double t1;
    double t2;
    double t3;
    double dt;
    double speed_sound = calc_speed_sound(data->Density, data->Pressure, ctx->Adiabatic);

    t1 = length/(fabs(data->u1) + speed_sound);
    t2 = length/(fabs(data->u2) + speed_sound);
    t3 = length/(fabs(data->u3) + speed_sound);

    dt = 1/(1/t1 + 1/t2 + 1/t3);

    ctx->dt = (ctx->dt == 0) ? dt : SC_MIN(ctx->dt, dt);
}

double calc_speed_sound(double density, double pressure, double adiabatic) {
    return sqrt(adiabatic * pressure / density);
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

    /* calc cfl */
    SC_PRODUCTION("Cell iter started\n");
    p8est_iterate(p8est, NULL, NULL,        /* слой гостовых ячеек не нужен, доп. параметры тоже */
                  calc_cfl_timestep_iter,   /* вычисление временного шага, проходя по всем ячейкам */
                  NULL, NULL, NULL);
    SC_PRODUCTION("Cell iter ended\n");

    /* calc min dt */
    SC_PRODUCTIONF("dt old: %.20lf\n", ctx->dt);
    mpiret = sc_MPI_Allreduce(MPI_IN_PLACE, &ctx->dt, 1, MPI_DOUBLE, MPI_MIN, p8est->mpicomm);
    SC_CHECK_MPI(mpiret);
    SC_PRODUCTIONF("dt new: %.20lf\n", ctx->dt);

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

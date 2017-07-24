#include "solver.h"
#include "../util.h"

void cflq(element_data_t *data, context_t *ctx, double length) {
    P4EST_ASSERT(length > 0);

    double t1;
    double t2;
    double t3;
    double dt;
    double speed = calc_speed(data->Z.Density, data->Z.Pressure, ctx->Adiabatic);

    t1 = length/(fabs(data->Z.u1) + speed);
    t2 = length/(fabs(data->Z.u2) + speed);
    t3 = length/(fabs(data->Z.u3) + speed);

    if (isnan(t1)) {
        t1 = 1;
    }

    if (isnan(t2)) {
        t2 = 1;
    }

    if (isnan(t3)) {
        t3 = 1;
    }

    dt = 1/(1/t1 + 1/t2 + 1/t3);

    ctx->dt = (ctx->dt == 0) ? dt : SC_MIN(ctx->dt, dt);
}

void calc_cfl_timestep_iter(p8est_iter_volume_info_t *info,
                            void *user_data) { // user_data (в т.ч. госты не используется при вычислении временного шага)
    p8est_t             *p8est     = info->p4est;
    context_t           *ctx       = (context_t *) p8est->user_pointer;  /* весь контекст на этом процессоре */
    p8est_quadrant_t    *q         = info->quad;                         /* текущая ячейка */

    element_data_t      *data;

    double              h;
    double              p[3];

    /* Длина стороны */
    h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;
    get_midpoint(p8est, info->treeid, q, p);

    data = (element_data_t *) q->p.user_data;

    /* Вычисление курантового шага */
    cflq(data, ctx, h);
}
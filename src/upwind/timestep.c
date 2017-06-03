#include "solver.h"
#include "../util.h"

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

    /* Вычисление шага */
    //cflq(data, ctx, h);
}

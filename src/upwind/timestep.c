#include "solver.h"
#include "../util.h"

void timestep_update_volume_iter(p8est_iter_volume_info_t *info,
                                 void *user_data) {
    p8est_quadrant_t   *q = info->quad;
    element_data_t     *data = (element_data_t *) q->p.user_data;
    double             dt = *((double *) user_data);
    double             h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;

    data->u += dt * data->dudt / (h * h * h);
}

double timestep_update_mesh_iter(p8est_t *p8est,
                                 p8est_mesh_t *mesh) {
    int                 qumid;
    int                 which_tree;

    p4est_locidx_t      quadrant_id;
    p8est_quadrant_t    *q;
    element_data_t      *data;

    context_t           *ctx = (context_t *) p8est->user_pointer;

    double              h;

    p4est_topidx_t      t, flt, llt;
    p8est_tree_t        *tree;
    int                 max_level, global_max_level;
    int                 mpiret, i;
    double              min_h, vnorm;
    double              dt;

    /* compute the timestep by finding the smallest quadrant */
    max_level = 0;

    // прохождение по всем локальным ячейкам
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        which_tree = -1;
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);

        max_level = SC_MAX(max_level, q->level);
    }

    mpiret = sc_MPI_Allreduce (&max_level, &global_max_level, 1, sc_MPI_INT,
                              sc_MPI_MAX, p8est->mpicomm);
    SC_CHECK_MPI (mpiret);

    min_h = (double) P8EST_QUADRANT_LEN (global_max_level) / (double) P8EST_ROOT_LEN;

    vnorm = 0;
    for (i = 0; i < P8EST_DIM; i++) {
        vnorm += ctx->v[i] * ctx->v[i];
    }
    vnorm = sqrt (vnorm);

    dt = min_h / 2. / vnorm;

    // прохождение по всем локальным ячейкам
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        which_tree = -1;
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);

        h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;

        data = (element_data_t *) q->p.user_data;
        data->u += dt * data->dudt / (h * h * h);
    }

    return dt;
}

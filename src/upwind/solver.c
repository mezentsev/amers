#include <math.h>
#include <sc.h>
#include "solver.h"

void init_solver(p8est_t *p8est,
                 p4est_topidx_t which_tree,
                 p8est_quadrant_t *quad) {
    double          midpoint[3];
    context_t       *ctx = (context_t *) p8est->user_pointer;
    element_data_t  *data = (element_data_t *) quad->p.user_data;
    double          r2, dx, dy, dz;

    get_midpoint(p8est, which_tree, quad, midpoint);

    data->du[0] = -1;
    data->du[1] = -1;
    data->du[2] = -1;

    // Compute the value and derivatives of the initial condition.
    r2 = 0.;

    // расстояние до центра
    dx = midpoint[0] - ctx->center[0];
    dy = midpoint[1] - ctx->center[1];
    dz = midpoint[2] - ctx->center[2];

    // радиус
    r2 += dx * dx + dy * dy + dz * dz;

    data->u = exp(-(1./2.) * r2 / (ctx->width * ctx->width));

    data->du[0] = -(1./ (ctx->width * ctx->width)) * dx * data->u;
    data->du[1] = -(1./ (ctx->width * ctx->width)) * dy * data->u;
    data->du[2] = -(1./ (ctx->width * ctx->width)) * dz * data->u;

    data->dummy = 0;
}

element_data_t get_boundary_data_by_face(p8est_t *p8est,
                                         p8est_quadrant_t *q,
                                         int face) {
    context_t *ctx = (context_t *) p8est->user_pointer;
    element_data_t *data = (element_data_t *) q->p.user_data;
    element_data_t boundary_data;

    boundary_data.dummy = 1;

    /* координаты ячейки сейчас не учитываются */

    /* направление face от соседа к текущей ячейке */
    switch (face) {
        case 0:                      /* -x side */
            SC_PRODUCTION("-x side\n");
            break;
        case 1:                      /* +x side */
            SC_PRODUCTION("+x side\n");
            break;
        case 2:                      /* -y side */
            SC_PRODUCTION("-y side\n");
            break;
        case 3:                      /* +y side */
            SC_PRODUCTION("+y side\n");
            break;
        case 4:                      /* -z side */
            SC_PRODUCTION("-z side\n");
            break;
        case 5:                      /* +z side */
            SC_PRODUCTION("+z side\n");
            break;
        default:
            SC_ABORT("Wrong face");
    }

    return boundary_data;
}

void reset_derivatives(p8est_t *p8est,
                       p8est_mesh_t *mesh) {
    int qumid;
    int which_tree;

    p4est_locidx_t quadrant_id;
    p8est_quadrant_t *q;
    element_data_t *data;

    // прохождение по всем локальным ячейкам
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        which_tree = -1;
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);

        data = (element_data_t *) q->p.user_data;
        data->du[0] = -1;
        data->du[1] = -1;
        data->du[2] = -1;
    }
}

void quad_divergence (p8est_t *p8est,
                      p8est_mesh_t *mesh)
{
    int qumid;
    int which_tree;

    p4est_locidx_t quadrant_id;
    p8est_quadrant_t *q;
    element_data_t *data;

    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        which_tree = -1;
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);

        data = (element_data_t *) q->p.user_data;
        data->dudt = 0.;
    }
}

void minmod_estimate(p8est_iter_face_info_t *info, void *user_data)
{
    int                 i, j;
    p8est_iter_face_side_t *side[2];
    sc_array_t          *sides = &(info->sides);
    element_data_t      *ghost_data = (element_data_t *) user_data;
    element_data_t      *udata;
    p8est_quadrant_t    *quad;
    double              uavg[2];
    double              h[2];
    double              du_est, du_old;
    int                 which_dir;

    side[0] = p8est_iter_fside_array_index_int (sides, 0);
    if (sides->elem_count == 1) {
        /* наткнулись на границу, копируем текущую сторону */
        side[1] = p8est_iter_fside_array_index_int (sides, 0);
    } else {
        side[1] = p8est_iter_fside_array_index_int (sides, 1);
    }

    which_dir = side[0]->face / 2;        /* 0 == x, 1 == y, 2 == z */

    for (i = 0; i < 2; i++) {
        uavg[i] = 0;
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P8EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];
                h[i] =
                        (double) P8EST_QUADRANT_LEN (quad->level) / (double) P8EST_ROOT_LEN;
                if (side[i]->is.hanging.is_ghost[j]) {
                    udata = &ghost_data[side[i]->is.hanging.quadid[j]];
                }
                else {
                    udata = (element_data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                }
                uavg[i] += udata->u;
            }
            uavg[i] /= P4EST_HALF;
        }
        else {
            quad = side[i]->is.full.quad;
            h[i] = (double) P8EST_QUADRANT_LEN (quad->level) / (double) P8EST_ROOT_LEN;
            if (side[i]->is.full.is_ghost) {
                udata = &ghost_data[side[i]->is.full.quadid];
            }
            else {
                udata = (element_data_t *) side[i]->is.full.quad->p.user_data;
            }
            uavg[i] = udata->u;
        }
    }
    du_est = (uavg[1] - uavg[0]) / ((h[0] + h[1]) / 2.);
    for (i = 0; i < 2; i++) {
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P8EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];
                if (!side[i]->is.hanging.is_ghost[j]) {
                    udata = (element_data_t *) quad->p.user_data;
                    du_old = udata->du[which_dir];

                    // TODO wtf?
                    if (du_old == du_old) {
                        /* there has already been an update */
                        if (du_est * du_old >= 0.) {
                            if (fabs (du_est) < fabs (du_old)) {
                                udata->du[which_dir] = du_est;
                            }
                        }
                        else {
                            udata->du[which_dir] = 0.;
                        }
                    }
                    else {
                        udata->du[which_dir] = du_est;
                    }
                }
            }
        }
        else {
            quad = side[i]->is.full.quad;
            if (!side[i]->is.full.is_ghost) {
                udata = (element_data_t *) quad->p.user_data;
                du_old = udata->du[which_dir];

                // TODO wtf?
                if (du_old == du_old) {
                    /* there has already been an update */
                    if (du_est * du_old >= 0.) {
                        if (fabs (du_est) < fabs (du_old)) {
                            udata->du[which_dir] = du_est;
                        }
                    }
                    else {
                        udata->du[which_dir] = 0.;
                    }
                }
                else {
                    udata->du[which_dir] = du_est;
                }
            }
        }
    }
}

void interpolate_solution(p8est_iter_volume_info_t *info,
                          void *user_data) {
    sc_array_t          *u_interp = (sc_array_t *) user_data;      /* we passed the array of values to fill as the user_data in the call to p4est_iterate */
    p8est_t             *p4est = info->p4est;
    p8est_quadrant_t    *q = info->quad;
    p4est_topidx_t      which_tree = info->treeid;
    p4est_locidx_t      local_id = info->quadid;  /* this is the index of q *within its tree's numbering*.  We want to convert it its index for all the quadrants on this process, which we do below */
    p8est_tree_t        *tree;
    element_data_t      *data = (element_data_t *) q->p.user_data;
    double              h;
    p4est_locidx_t      arrayoffset;
    double              this_u;
    double              *this_u_ptr;
    int                 i, j;

    tree = p8est_tree_array_index (p4est->trees, which_tree);
    local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
    arrayoffset = P8EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */
    //h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;

    for (i = 0; i < P8EST_CHILDREN; i++) {
        this_u = data->u;
        /* loop over the derivative components and linearly interpolate from the
         * midpoint to the corners */
        //for (j = 0; j < P8EST_DIM; j++) {
            /* In order to know whether the direction from the midpoint to the corner is
             * negative or positive, we take advantage of the fact that the corners
             * are in z-order.  If i is an odd number, it is on the +x side; if it
             * is even, it is on the -x side.  If (i / 2) is an odd number, it is on
             * the +y side, etc. */
        //    this_u += (h / 2) * data->du[j] * ((i & (1 << j)) ? 1. : -1.);
        //}
        this_u_ptr = (double *) sc_array_index (u_interp, arrayoffset + i);
        this_u_ptr[0] = this_u;
    }

}

void write_solution(p8est_t *p8est,
                    int step) {
    char                filename[BUFSIZ] = { '\0' };
    p4est_locidx_t      numquads;
    sc_array_t          *data;

    snprintf (filename, 12, "solution_%02d", step);
    numquads = p8est->local_num_quadrants;

    /* create a vector with one value for the corner of every local quadrant
     * (the number of children is always the same as the number of corners) */
    data = sc_array_new_size (sizeof (double), numquads * P8EST_CHILDREN);

    /* Use the iterator to visit every cell and fill in the solution values.
     * Using the iterator is not absolutely necessary in this case: we could
     * also loop over every tree (there is only one tree in this case) and loop
     * over every quadrant within every tree */
    p8est_iterate(p8est, NULL,   /* we don't need any ghost quadrants for this loop */
                  (void *) data,     /* pass in boundary so that we can fill it */
                  interpolate_solution,
                  NULL,          /* there is no callback for the faces between quadrants */
                  NULL,          /* there is no callback for the edges between quadrants */
                  NULL);         /* there is no callback for the corners between quadrants */

    /* create VTK output context and set its parameters */
    p8est_vtk_context_t *context = p8est_vtk_context_new (p8est, filename);
    p8est_vtk_context_set_scale (context, 0.99);  /* quadrant at almost full scale */

    /* begin writing the output files */
    context = p8est_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL,
                    P8EST_STRING "_vtk: Error writing vtk header");

    /* do not write the tree id's of each quadrant
     * (there is only one tree in this example) */
    context = p8est_vtk_write_cell_dataf(context, 0, 1,  /* do write the refinement level of each quadrant */
                                         1,      /* do write the mpi process id of each quadrant */
                                         0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                         0,      /* there is no custom cell scalar data. */
                                         0,      /* there is no custom cell vector data. */
                                         context);       /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P8EST_STRING "_vtk: Error writing cell data");

    /* write one scalar field: the solution value */
    context = p8est_vtk_write_point_dataf(context, 1, 0,
                                          "u",
                                          data,
                                          context);
    SC_CHECK_ABORT (context != NULL,
                    P8EST_STRING "_vtk: Error writing cell data");

    const int retval = p8est_vtk_write_footer(context);
    SC_CHECK_ABORT (!retval,
                    P8EST_STRING "_vtk: Error writing footer");

    sc_array_destroy(data);
}

void solver_step(p8est_t *p8est,
                 int step) {
    SC_PRODUCTIONF("Start solve step %d\n", step);

    int                 mpiret;
    element_data_t      *ghost_data;
    p8est_ghost_t       *ghost;
    p8est_mesh_t        *mesh;
    double              dt = 0.;
    context_t           *ctx = (context_t *) p8est->user_pointer;

    /* выделение гостового слоя */
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC (element_data_t, ghost->ghosts.elem_count);
    mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FACE);

    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    SC_PRODUCTION("Reset derivatives started\n");
    reset_derivatives(p8est, mesh);
    SC_PRODUCTION("Reset derivatives ended\n");

    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    SC_PRODUCTION("Minmod estimate started\n");
    p8est_iterate(p8est, ghost, (void *) ghost_data,
                  NULL,
                  minmod_estimate,         /*  */
                  NULL, NULL);
    SC_PRODUCTION("Minmod estimate ended\n");

    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    /* calc timestep */
    //SC_PRODUCTION("Cell iter started\n");
    //p8est_iterate(p8est, NULL, NULL,        /* слой гостовых ячеек не нужен, доп. параметры тоже */
    //              timestep_update_volume_iter,        /* вычисление временного шага, проходя по всем ячейкам */
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

    SC_PRODUCTION("Divergence started\n");
    quad_divergence(p8est, mesh);
    SC_PRODUCTION("Divergence ended\n");

    /* update du/dx estimate */
    SC_PRODUCTION("Calc flux started\n");
    p8est_iterate(p8est, ghost, (void *) ghost_data,        /* вкладываем гостовый слой */
                  NULL,
                  calc_flux_face_iter,                      /* обход по фейсам каждой ячейки для высчитывания потока */
                  NULL,
                  NULL);
    //calc_flux_mesh_iter(p8est, mesh, ghost, ghost_data);
    SC_PRODUCTION("Calc flux ended\n");

    /* exchange ghost data */
    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    /* update u */
    SC_PRODUCTION("Calc flux started\n");
    dt = timestep_update_mesh_iter(p8est, mesh);
    SC_PRODUCTIONF("dt is %lf\n", dt);
    SC_PRODUCTION("Calc flux started\n");

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

    ghost_data = NULL;
    ghost = NULL;

    /* перераспределение между процессорами */
    SC_PRODUCTION("Partition started\n");
    p8est_partition (p8est, 1, NULL);
    SC_PRODUCTION("Partition ended\n");

    SC_PRODUCTIONF("End solve step %d\n", step);
}
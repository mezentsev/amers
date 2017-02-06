#include <p8est_iterate.h>
#include "amr.h"

/**
 * Debug print quadrant
 * @param q
 * @param rank
 */
void
quadrant_print(p8est_quadrant_t *q, int rank)
{
    p4est_qcoord_t x = (q->x) >> (P4EST_MAXLEVEL - q->level);
    p4est_qcoord_t y = (q->y) >> (P4EST_MAXLEVEL - q->level);
    p4est_qcoord_t z = (q->z) >> (P4EST_MAXLEVEL - q->level);
    printf("[p4est %d] x=%d y=%d z=%d level=%d\n", rank, x, y, z, q->level);
}

/**
 * Calculate laplacian for input function f
 * http://mathworld.wolfram.com/Laplacian.html
 *
 * @param f
 * @param x
 * @param y
 * @param z
 * @param dt
 * @return
 */
double
laplacian(t_func_3 f,
          double x, double y, double z)
{
    const int       n = 10;
    const double    dt = 0.01;
    double          sum = 0;
    int             i;

    for (i = 0; i < n; ++i)
    {
        sum += 1/pow(dt,2.) * (f(x + dt, y, z) + f(x - dt, y, z) +
                               f(x, y + dt, z) + f(x, y - dt, z) +
                               f(x, y, z + dt) + f(x, y, z - dt) -
                               6 * f(x, y, z));
    }

    sum /= n;

    return sum;
}

/**
 * Calculate gradient
 * http://mathworld.wolfram.com/Gradient.html
 *
 * @return
 */
double
grad(t_func_3 f,
     double x, double y, double z,
     double nx, double ny, double nz)
{
    const int       n = 10;
    const double    dt = 0.01;
    double          sum = 0;
    int             i;

    for (i = 0; i < n; ++i)
    {
        sum +=  (nx * (f(x+dt, y, z) - f(x,y,z)) / dt / 2 +
                 ny * (f(x, y+dt, z) - f(x,y,z)) / dt / 2 +
                 nz * (f(x, y, z+dt) - f(x,y,z)) / dt / 2);
    }

    sum /= n;

    return sum;
}

/**
 * Get the coordinates of the midpoint of a quadrant.
 */
void
get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3])
{
    p4est_qcoord_t half_length = P8EST_QUADRANT_LEN(q->level) / 2;

    p8est_qcoord_to_vertex (p8est->connectivity,
                            tree,
                            q->x + half_length,
                            q->y + half_length,
                            q->z + half_length,
                            xyz);
}

/**
 * Initialize each cell
 *
 * @param p4est
 * @param which_tree
 * @param q
 */
void
init(p8est_t *p4est, p4est_topidx_t which_tree, p8est_quadrant_t *q)
{
    /* the data associated with a forest is accessible by user_pointer */
    ctx_t       *ctx = (ctx_t *) p4est->user_pointer;
    /* the data associated with a quadrant is accessible by p.user_data */
    data_t      *data = (data_t *) q->p.user_data;

    /* initialize the data */
    data->f1 = 0;
    data->f2 = 0;
    data->v  = 0;
    data->e  = 0;
}

/**
 * Refine every step (non-recursive)
 *
 * @param p8est
 * @param which_tree
 * @param q
 * @return
 */
int
refine_always (p8est_t *p8est,
               p4est_topidx_t which_tree,
               p8est_quadrant_t *q)
{
    return 1;
}

/**
 * Refine rule
 *
 * @param p8est
 * @param which_tree
 * @param q
 * @return
 */
int
refine_fn (p8est_t *p8est,
           p4est_topidx_t which_tree,
           p8est_quadrant_t *q)
{
    ctx_t *ctx = (ctx_t *) p8est->user_pointer;
    double midpoint[3];
    double h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;
    double l;

    get_midpoint(p8est, which_tree, q, midpoint);

    l = pow(midpoint[0] - ctx->center[0], 2.) +
        pow(midpoint[1] - ctx->center[1], 2.) +
        pow(midpoint[2] - ctx->center[2], 2.);

    if (l < ctx->width &&
        h > pow(ctx->width, 2.))
        return 1;

    return 0;
}

int
coarsen_fn (p8est_t *p8est,
            p4est_topidx_t which_tree,
            p8est_quadrant_t **children)
{
    return 0;
}

/**
 * Calculate F1
 *
 * The function p4est_iterate() takes as an argument a p4est_iter_volume_t
 * callback function, which it executes at every local quadrant
 **/
void
volume_iter(p8est_iter_volume_info_t *info, void *user_data)
{
    p8est_quadrant_t    *q = info->quad;
    data_t              *data = (data_t *) q->p.user_data;
    double              h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;
    double              mp[3];
    ctx_t               *ctx = (ctx_t *) info->p4est->user_pointer;

    get_midpoint(info->p4est, info->treeid, q, mp);

    // Laplasian * V
    data->f2 = ctx->f(mp[0], mp[1], mp[2]) * h * h * h;
}

/**
 * The function p4est_iterate() takes as an argument a p4est_iter_face_info_t
 * callback function, which it executes for each face
 *
 * \param [in] info the information about the quadrants on either side of the
 *                  interface, populated by p4est_iterate()
 * \param [in] user_data the user_data given to p4est_iterate(): in this case,
 *                       it points to the ghost_data array, which contains the
 *                       step3_data_t data for all of the ghost cells, which
 *                       was populated by p4est_ghost_exchange_data()
 */
void
face_iter_f1(p8est_iter_face_info_t *info, void *user_data)
{
    int                 i, j;
    p8est_t             *p4est = info->p4est;
    ctx_t               *ctx = (ctx_t *) p4est->user_pointer;
    data_t              *ghost_data = (data_t *) user_data;
    data_t              *udata;
    p8est_quadrant_t    *quad;
    double              vdotn = 0.;
    double              uavg[2];
    double              q;
    double              facearea;
    p4est_qcoord_t      h;
    int                 which_side;
    int                 which_dir;
    double              x, y, z, midpoint[3];
    p8est_iter_face_side_t *side[2];
    sc_array_t          *sides = &(info->sides);

    /* because there are no boundaries, every face has two sides */
    P4EST_ASSERT (sides->elem_count == 2);

    side[0] = p8est_iter_fside_array_index_int(sides, 0);
    side[1] = p8est_iter_fside_array_index_int(sides, 1);

    which_dir = side[0]->face;

    /* Because we have non-conforming boundaries, one side of an interface can
     * either have one large ("full") quadrant or 2^(d-1) small ("hanging")
     * quadrants: we have to compute the average differently in each case.  The
     * info populated by p4est_iterate() gives us the context we need to
     * proceed. */
    for(i = 0; i < 2; ++i)
    {
        uavg[i] = 0;
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P8EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];
                if (side[i]->is.hanging.is_ghost[j]) {
                    udata = (data_t *) &ghost_data[side[i]->is.hanging.quadid[j]];
                }
                else {
                    udata = (data_t *) side[i]->is.hanging.quad[j]->p.user_data;
                }

                /* mid point from neighbor cell */
                get_midpoint(p4est,
                             quad->p.which_tree,
                             quad,
                             midpoint);
                h = P8EST_QUADRANT_LEN(side[i]->is.hanging.quad[j]->level) / 2;

                switch (which_dir) {
                    case 0:                      /* -x side */
                        x = midpoint[0] + h / 2;
                        y = midpoint[1];
                        z = midpoint[2];
                        break;
                    case 1:                      /* +x side */
                        x = midpoint[0] - h / 2;
                        y = midpoint[1];
                        z = midpoint[2];
                        break;
                    case 2:                      /* -y side */
                        x = midpoint[0];
                        y = midpoint[1] + h / 2;
                        z = midpoint[2];
                        break;
                    case 3:                      /* +y side */
                        x = midpoint[0];
                        y = midpoint[1] - h / 2;
                        z = midpoint[2];
                        break;
                    case 4:                      /* -z side */
                        x = midpoint[0];
                        y = midpoint[1];
                        z = midpoint[2] + h / 2;
                        break;
                    case 5:                      /* +z side */

                        x = midpoint[0];
                        y = midpoint[1];
                        z = midpoint[2] - h / 2;
                        break;
                }

                uavg[i] += grad(ctx->f,
                                midpoint[0], midpoint[1], midpoint[2],
                                x, y, z) * 4 * h * h; // s
            }
        }
        else {
            quad = side[i]->is.full.quad;
            if (side[i]->is.full.is_ghost) {
                udata = (data_t *) &ghost_data[side[i]->is.full.quadid];
            }
            else {
                udata = (data_t *) side[i]->is.full.quad->p.user_data;
                //quadrant_print(side[which_side]->is.full.quad, 0);
            }
            get_midpoint(p4est,
                         quad->p.which_tree,
                         side[i]->is.full.quad,
                         midpoint);
            uavg[i] = grad(ctx->f,
                           midpoint[0], midpoint[1], midpoint[2],
                           midpoint[0], midpoint[1], midpoint[2]);
        }
    }

    for (i = 0; i < 2; i++) {
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P8EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];

                if (!side[i]->is.hanging.is_ghost[j]) {
                    udata = (data_t *) quad->p.user_data;
                    udata->f1 = uavg[i];
                }
            }
        }
        else {
            quad = side[i]->is.full.quad;
            if (!side[i]->is.full.is_ghost) {
                udata = (data_t *) quad->p.user_data;
                udata->f1 = uavg[i];
            }
        }
    }
}

void
calc_error_iter(p8est_iter_volume_info_t *info, void *user_data)
{

    p8est_quadrant_t    *q = info->quad;
    data_t              *data = (data_t *) q->p.user_data;

    data->e = fabs(data->f2 - data->f1);
}

/**
 * Make 1d-array with solution from all processors
 * @param info
 * @param user_data
 * @param f
 */
void
get_solution(p8est_iter_volume_info_t *info, void *user_data, int f)
{
    sc_array_t          *u_interp = (sc_array_t *) user_data;
    /* we passed the array of values to fill as the
     * user_data in the call to p4est_iterate */

    p8est_t             *p8est = info->p4est;
    p8est_quadrant_t    *q = info->quad;
    p4est_topidx_t      which_tree = info->treeid;
    p4est_locidx_t      local_id = info->quadid;
    /* this is the index of q *within its tree's numbering*.
     * We want to convert it its index for all the quadrants on this process, which we do below */

    p8est_tree_t        *tree;
    data_t              *data = (data_t *) q->p.user_data;
    p4est_locidx_t      array_offset;
    double              *this_u_ptr;
    int                 i;

    tree = p8est_tree_array_index (p8est->trees, which_tree);
    local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
    array_offset = P8EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */

    for (i = 0; i < P8EST_CHILDREN; i++) {
        this_u_ptr = (double *) sc_array_index (u_interp, array_offset + i);

        switch (f)
        {
            case 1:
                this_u_ptr[0] = data->f1;
                break;
            case 2:
                this_u_ptr[0] = data->f2;
                break;
            case 0:
            default:
                //SC_PRODUCTIONF("e=%lf: f2=%lf, f1=%lf\n", data->e, data->f2, data->f1);
                this_u_ptr[0] = (data->e == data->f2) ? 0 : data->e;
                break;
        }
    }
}

/**
 * Fill user_data with scalar data for function F1
 *
 * The function p4est_iterate() takes as an argument a p4est_iter_volume_t
 * callback function, which it executes at every local quadrant
 **/
void
get_solution_f1(p8est_iter_volume_info_t *info, void *user_data)
{
    get_solution(info, user_data, 1);
}

/**
 * Fill user_data with scalar data for function F2
 *
 * The function p4est_iterate() takes as an argument a p4est_iter_volume_t
 * callback function, which it executes at every local quadrant
 **/
void
get_solution_f2(p8est_iter_volume_info_t *info, void *user_data)
{
    get_solution(info, user_data, 2);
}

/**
 * Fill user_data with scalar data for error
 *
 * The function p4est_iterate() takes as an argument a p4est_iter_volume_t
 * callback function, which it executes at every local quadrant
 **/
void
get_error_estimate(p8est_iter_volume_info_t *info, void *user_data)
{
    get_solution(info, user_data, 0);
}


/* Create vtk with solution */
void
write_solution(p8est_t *p8est, int step)
{
    char                filename[BUFSIZ] = { '\0' };
    size_t              i;
    sc_array_t          *f1_interp;
    sc_array_t          *f2_interp;
    sc_array_t          *error_estimate_interp;
    double              *err_ptr;
    double              error;
    p4est_locidx_t      numquads;

    snprintf (filename, 12, "solution_%02d", step);
    error = 0;

    numquads = p8est->local_num_quadrants;

    /* create a vector with one value for the corner of every local quadrant
     * (the number of children is always the same as the number of corners) */
    f1_interp = sc_array_new_size (sizeof (double), numquads * P8EST_CHILDREN);
    f2_interp = sc_array_new_size (sizeof (double), numquads * P8EST_CHILDREN);
    error_estimate_interp = sc_array_new_size (sizeof (double), numquads * P8EST_CHILDREN);

    /* Use the iterator to visit every cell and fill in the solution values.
     * Using the iterator is not absolutely necessary in this case: we could
     * also loop over every tree (there is only one tree in this case) and loop
     * over every quadrant within every tree */
    p4est_iterate(p8est, NULL,   /* we don't need any ghost quadrants for this loop */
                  (void *) f1_interp,     /* pass in u_interp so that we can fill it */
                  get_solution_f1,    /* callback function that fill u_interp */
                  NULL,          /* there is no callback for the faces between quadrants */
                  NULL,          /* there is no callback for the edges between quadrants */
                  NULL);         /* there is no callback for the corners between quadrants */

    p4est_iterate(p8est, NULL,   /* we don't need any ghost quadrants for this loop */
                  (void *) f2_interp,     /* pass in u_interp so that we can fill it */
                  get_solution_f2,    /* callback function that fill u_interp */
                  NULL,          /* there is no callback for the faces between quadrants */
                  NULL,          /* there is no callback for the edges between quadrants */
                  NULL);         /* there is no callback for the corners between quadrants */

    p4est_iterate(p8est, NULL,   /* we don't need any ghost quadrants for this loop */
                  (void *) error_estimate_interp,     /* pass in u_interp so that we can fill it */
                  get_error_estimate,    /* callback function that fill u_interp */
                  NULL,          /* there is no callback for the faces between quadrants */
                  NULL,          /* there is no callback for the edges between quadrants */
                  NULL);         /* there is no callback for the corners between quadrants */

    for (i = 0; i < error_estimate_interp->elem_count; ++i)
    {
        err_ptr = (double *) sc_array_index(error_estimate_interp, i);
        error += err_ptr[0];
    }

    printf("Error estimate: %lf\n", error);

    /* create VTK output context and set its parameters */
    p4est_vtk_context_t *context = p4est_vtk_context_new (p8est, filename);
    p4est_vtk_context_set_scale (context, 0.99);  /* quadrant at almost full scale */

    /* begin writing the output files */
    context = p4est_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing vtk header");

    /* do not write the tree id's of each quadrant
     * (there is only one tree in this example) */
    context = p8est_vtk_write_cell_dataf(context, 0, 1,  /* do write the refinement level of each quadrant */
                                          1,      /* do write the mpi process id of each quadrant */
                                          0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                          0,      /* there is no custom cell scalar data. */
                                          0,      /* there is no custom cell vector data. */
                                          context);       /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    /* write one scalar field: the solution value */
    /*context = p8est_vtk_write_point_dataf(context, 3, 0,
                                          "solution_f1",
                                          f1_interp,
                                          "solution_f2",
                                          f2_interp,
                                          "error_estimate",
                                          error_estimate_interp,
                                          context);
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");*/

    const int retval = p4est_vtk_write_footer(context);
    SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

    sc_array_destroy (f1_interp);
    sc_array_destroy (f2_interp);
    sc_array_destroy (error_estimate_interp);
}

void
solve(p8est_t *p8est, int step)
{
    data_t              *ghost_data;
    ctx_t               *ctx = (ctx_t *) p8est->user_pointer;
    p8est_ghost_t       *ghost;

    /* create the ghost quadrants */
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FULL);
    /* create space for storing the ghost data */
    ghost_data = P4EST_ALLOC (data_t, ghost->ghosts.elem_count);
    /* synchronize the ghost data */
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);

    p8est_iterate(p8est, ghost, (void *) ghost_data,
                  volume_iter,
                  face_iter_f1,
                  NULL,
                  NULL
    );

    p8est_ghost_exchange_data (p8est, ghost, ghost_data);

    p8est_iterate(p8est, ghost, (void *) ghost_data,
                  calc_error_iter,
                  NULL,
                  NULL,
                  NULL
    );

    /* Test face neighbor iterator
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        which_tree = -1;
        q = p4est_mesh_quadrant_cumulative (p4est, qumid,
                                            &which_tree, &quadrant_id);
        p4est_mesh_face_neighbor_init2 (&mfn, p4est, ghost, mesh,
                                        which_tree, quadrant_id);
        p4est_mesh_face_neighbor_init (&mfn2, p4est, ghost, mesh, which_tree, q);
        P4EST_ASSERT (mfn2.quadrant_id == quadrant_id);
        while ((q = p4est_mesh_face_neighbor_next (&mfn, &which_tree, &which_quad,
                                                   &nface, &nrank)) != NULL) {

            data_t *data;

            data = (data_t *) p4est_mesh_face_neighbor_data (&mfn, ghost_data);

            P4EST_ASSERT (p4est_quadrant_is_equal (q, &(data->quad)));
            P4EST_ASSERT (data->quad.p.which_tree == which_tree);

        }
    }*/

    write_solution(p8est, step);

    /* clear */
    P4EST_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}

int
main (int argc, char **argv)
{
    int                     mpiret;
    sc_MPI_Comm             mpicomm;
    ctx_t                   ctx;
    p8est_t                 *p8est;
    p8est_connectivity_t    *conn;

    ctx.center[0] = 0.5;
    ctx.center[1] = 0.5;
    ctx.center[2] = 0.5;

    ctx.width = 0.1;
    ctx.f = &s_func;

    /* Initialize MPI; see sc_mpi.h.
     * If configure --enable-mpi is given these are true MPI calls.
     * Else these are dummy functions that simulate a single-processor run. */
    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    /* These functions are optional.  If called they store the MPI rank as a
     * static variable so subsequent global p4est log messages are only issued
     * from processor zero.  Here we turn off most of the logging; see sc.h. */
    sc_init (mpicomm, 1, 1, NULL, SC_LP_PRODUCTION);
    p4est_init(NULL, SC_LP_PRODUCTION);

    conn = p8est_connectivity_new_periodic();


    p8est = p8est_new_ext (mpicomm, /* communicator */
                           conn,    /* connectivity */
                           0,       /* minimum quadrants per MPI process */
                           4,       /* minimum level of refinement */
                           1,       /* fill uniform */
                           sizeof (data_t),         /* data size */
                           init,  /* initializes data */
                           (void *) (&ctx));              /* context */

    p8est_refine(p8est, 1, refine_fn, init);
    p8est_coarsen(p8est, 1, coarsen_fn, init);

    p8est_balance(p8est, P4EST_CONNECT_FACE, init);
    p8est_partition (p8est, 1, NULL);

    SC_PRODUCTIONF("Start solve 0\n", NULL);
    solve(p8est, 0);
    SC_PRODUCTIONF("End solve 0\n", NULL);

    p8est_refine(p8est, 0, refine_always, init);
    p8est_partition (p8est, 0, NULL);

    SC_PRODUCTIONF("Start solve 1\n", NULL);
    solve(p8est, 1);
    SC_PRODUCTIONF("End solve 1\n", NULL);

    /*p8est_refine(p8est, 0, refine_always, init);
    p8est_partition (p8est, 0, NULL);

    SC_PRODUCTIONF("Start solve 2\n", NULL);
    solve(p8est, 2);
    SC_PRODUCTIONF("End solve 2\n", NULL);

    p8est_refine(p8est, 0, refine_always, init);
    p8est_partition (p8est, 0, NULL);

    SC_PRODUCTIONF("Start solve 3\n", NULL);
    solve(p8est, 3);
    SC_PRODUCTIONF("End solve 3\n", NULL);*/

    // clear
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);

    /* Verify that allocations internal to p4est and sc do not leak memory.
     * This should be called if sc_init () has been called earlier. */
    sc_finalize ();

    /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
}

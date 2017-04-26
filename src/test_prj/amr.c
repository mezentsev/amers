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
init(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *quad)
{
    /* the data associated with a forest is accessible by user_pointer */
    ctx_t       *ctx = (ctx_t *) p8est->user_pointer;
    /* the data associated with a quadrant is accessible by p.user_data */
    data_t      *data = (data_t *) quad->p.user_data;
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
        l > 0.08 &&
        h > 0.02)
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
 * Iterate throw cells and neighbors
 *
 * @param p8est
 * @param ghost
 * @param mesh
 * @param ghost_data
 */
void
mesh_neighbors_iter(p8est_t *p8est, p8est_ghost_t *ghost, p8est_mesh_t *mesh, void *ghost_data)
{
    p8est_quadrant_t    *neighbor  = NULL;
    p4est_topidx_t      which_tree = -1;
    ctx_t               *ctx       = (ctx_t *) p8est->user_pointer;

    p8est_mesh_face_neighbor_t mfn;
    p4est_locidx_t      qumid, quadrant_id, which_quad;
    p8est_quadrant_t    *q;                             /* current quad */
    data_t              *g_data;                        /* ghost data */
    data_t              *data;

    int                 nface;
    int                 nrank;
    double              i = 0;
    double              n_data = 0;  /* neighbors data */
    double              h, Vi;
    double              nx, ny, nz;

    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid)
    {
        q = p8est_mesh_quadrant_cumulative (p8est, qumid,
                                            &which_tree, &quadrant_id);
        data = (data_t *) q->p.user_data;
        p8est_mesh_face_neighbor_init2 (&mfn,
                                        p8est,
                                        ghost,
                                        mesh,
                                        which_tree,
                                        quadrant_id);

        while((neighbor = p8est_mesh_face_neighbor_next (&mfn, &which_tree,
                                                         &which_quad, &nface, &nrank)) != NULL)
        {
            g_data = (data_t *) p8est_mesh_face_neighbor_data (&mfn, ghost_data);

            nx = 0;
            ny = 0;
            nz = 0;

            /* neighbor orientation */
            switch (nface) {
                case 0:                      /* -x side */
                    nx = -1;
                    break;
                case 1:                      /* +x side */
                    nx = 1;
                    break;
                case 2:                      /* -y side */
                    ny = -1;
                    break;
                case 3:                      /* +y side */
                    ny = 1;
                    break;
                case 4:                      /* -z side */
                    nz = -1;
                    break;
                case 5:                      /* +z side */
                    nz = 1;
                    break;
            }

            h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;
            Vi = pow(h, 3.);

            n_data += (nx * (Vi/(Vi+g_data->V) * g_data->dfdx + g_data->V/(Vi+g_data->V) * data->dfdx) +
                       ny * (Vi/(Vi+g_data->V) * g_data->dfdy + g_data->V/(Vi+g_data->V) * data->dfdy) +
                       nz * (Vi/(Vi+g_data->V) * g_data->dfdz + g_data->V/(Vi+g_data->V) * data->dfdz)) * g_data->S;
        }

        i += 1;

        data->f1 = n_data;
        n_data = 0;
    }
}

/**
 * Iterate throw volume cells
 *
 * @param p8est
 * @param mesh
 */
void
mesh_iter(p8est_t *p8est, p8est_mesh_t *mesh)
{
    p8est_quadrant_t    *neighbor  = NULL;
    p4est_topidx_t      which_tree = -1;
    ctx_t               *ctx       = (ctx_t *) p8est->user_pointer;

    p8est_mesh_face_neighbor_t mfn;
    p4est_locidx_t      qumid, quadrant_id, which_quad;
    p8est_quadrant_t    *q;                             /* current quad */
    data_t              *g_data;                        /* ghost data */
    data_t              *data;

    double              h, Vi, Si;
    double              p[3];

    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid)
    {
        q = p8est_mesh_quadrant_cumulative (p8est, qumid,
                                            &which_tree, &quadrant_id);


        /* side length */
        h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;
        Vi = pow(h, 3.);
        Si = pow(h, 2.);
        get_midpoint(p8est, which_tree, q, p);

        data = (data_t *) q->p.user_data;

        data->V     = Vi;
        data->S     = Si;
        data->f0    = pow(p[0], 2.) + pow(p[1], 2.) + pow(p[2], 2.); /** f    = x^2+y^2+z^2          */
        data->dfdx  = 2 * p[0];                                      /** dfdx = 2x                   */
        data->dfdy  = 2 * p[1];                                      /** dfdy = 2y                   */
        data->dfdz  = 2 * p[2];                                      /** dfdz = 2z                   */

        data->f2    = 6 * Vi;                                        /** laplacian x^2 + y^2 + z^2   */
    }
}

/**
 * Iterate throw volume cells
 *
 * @param p8est
 * @param ghost
 * @param mesh
 * @param ghost_data
 */
void
error_iter(p8est_t *p8est, p8est_ghost_t *ghost, p8est_mesh_t *mesh, void *ghost_data)
{
    p4est_topidx_t      which_tree = -1;
    ctx_t               *ctx       = (ctx_t *) p8est->user_pointer;

    int                 mpiret;
    double              h, Vi;
    double              err = 0;

    p8est_quadrant_t    *q;                             /* current quad */
    p8est_mesh_face_neighbor_t mfn;
    p4est_locidx_t      qumid, quadrant_id, which_quad;
    data_t              *data;

    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid)
    {
        q = p8est_mesh_quadrant_cumulative (p8est, qumid,
                                            &which_tree, &quadrant_id);
        data = (data_t *) q->p.user_data;
        p8est_mesh_face_neighbor_init2 (&mfn,
                                        p8est,
                                        ghost,
                                        mesh,
                                        which_tree,
                                        quadrant_id);

        h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;
        Vi = pow(h, 3.);

        data->e = pow(fabs(data->f1 - data->f2), 2.) * Vi;

        err += data->e;
        ctx->err = err;
    }
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
                this_u_ptr[0] = data->e;
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

void
solve(p8est_t *p8est, int step)
{
    data_t              *ghost_data;
    ctx_t               *ctx = (ctx_t *) p8est->user_pointer;
    p8est_ghost_t       *ghost;
    p8est_mesh_t        *mesh;

    /* create the ghost quadrants */
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC (data_t, ghost->ghosts.elem_count);
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);

    mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FACE);
    SC_PRODUCTIONF("Used memory: %ld\n", p8est_mesh_memory_used(mesh));

    /* calc f2 and partials */
    mesh_iter(p8est, mesh);
    /* exchange ghost data */
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    /* calc f1 */
    mesh_neighbors_iter(p8est, ghost, mesh, ghost_data);
    /* exchange ghost data */
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);

    /* calculate error on every cell */
    error_iter(p8est, ghost, mesh, ghost_data);

    /* generate vtk and print solution */
    write_solution(p8est, step);

    /* clear */
    p8est_mesh_destroy(mesh);
    P4EST_FREE (ghost_data);
    p4est_ghost_destroy (ghost);
}

/* Create vtk with solution */
void
write_solution(p8est_t *p8est, int step)
{
    char                filename[BUFSIZ] = { '\0' };
    ctx_t               *ctx             = (ctx_t *) p8est->user_pointer;
    size_t              i;
    sc_array_t          *f1_interp;
    sc_array_t          *f2_interp;
    sc_array_t          *error_estimate_interp;
    double              *err_ptr;
    double              error = 0;
    p4est_locidx_t      numquads;
    FILE                *err_out;
    int                 mpiret;

    snprintf (filename, 12, "solution_%02d", step);

    numquads = p8est->local_num_quadrants;

    /* create a vector with one value for the corner of every local quadrant
     * (the number of children is always the same as the number of corners) */
    f1_interp = sc_array_new_size (sizeof (double), numquads * P8EST_CHILDREN);
    //f2_interp = sc_array_new_size (sizeof (double), numquads * P8EST_CHILDREN);
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

    SC_PRODUCTIONF("Quadrants %d. Error estimate: %.20lf\n", p8est->local_num_quadrants, ctx->err);

    mpiret = MPI_Allreduce(&ctx->err, &error, 1, MPI_DOUBLE, MPI_SUM, p8est->mpicomm);
    SC_CHECK_MPI(mpiret);

    SC_PRODUCTIONF("Allreduce. Quadrants %d. Error estimate: %.20lf\n", p8est->local_num_quadrants, error);

    if (p8est->mpirank == 0)
    {
        err_out = fopen(filename, "w");
        SC_CHECK_ABORT(err_out, "Can't write to file");
        fprintf(err_out, "G%d: cells=%lld; err=%.20lf\n", step, p8est->global_num_quadrants, ctx->err);
        fclose(err_out);
    }

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
    context = p8est_vtk_write_point_dataf(context, 1, 0,
                                          "solution_f1",
                                          f1_interp,
                                          //"solution_f2",
                                          //f2_interp,
                                          //"error_estimate",
                                          //error_estimate_interp,
                                          context);
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    const int retval = p4est_vtk_write_footer(context);
    SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

    sc_array_destroy (f1_interp);
    //sc_array_destroy (f2_interp);
    sc_array_destroy (error_estimate_interp);
}

int
main (int argc, char **argv)
{
    char                    filename[BUFSIZ] = { '\0' };
    int                     mpiret;
    int                     i;
    sc_MPI_Comm             mpicomm;
    ctx_t                   ctx;
    p8est_t                 *p8est;
    p8est_connectivity_t    *conn;
    clock_t                 start;
    FILE                    *file;

    start = clock();

    ctx.center[0] = 0.5;
    ctx.center[1] = 0.5;
    ctx.center[2] = 0.5;

    ctx.width = 0.2;
    ctx.count = 0;
    ctx.err = 0;
    ctx.level = 4;
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
    //p8est_coarsen(p8est, 1, coarsen_fn, init);

    p8est_balance(p8est, P4EST_CONNECT_FULL, init);
    p8est_partition (p8est, 1, NULL);

    // init time
    if(p8est->mpirank == 0) {
        file = fopen("init.time", "w");
        SC_CHECK_ABORT(file, "Can't write to file");
        fprintf(file, "Init took %f seconds\n", ((double) clock() - start) / CLOCKS_PER_SEC);
        fclose(file);
    }

    p8est_vtk_write_file (p8est, NULL, "init");

    for (i = 0; i < ctx.level; ++i)
    {
        start = clock();

        SC_PRODUCTIONF("Start solve %d step\n", i);
        solve(p8est, i);
        SC_PRODUCTIONF("End solve %d step\n", i);

        // solve time for i-step
        if(p8est->mpirank == 0)
        {
            snprintf (filename, 17, "solution_%02d.time", i);
            file = fopen(filename, "w");
            SC_CHECK_ABORT(file, "Can't write to file");
            fprintf(file, "Solution took %f seconds\n", ((double)clock() - start)/CLOCKS_PER_SEC);
            fclose(file);
        }
        //

        if (i == ctx.level - 1)
            break;

        start = clock();

        p8est_refine(p8est, 0, refine_always, init);
        p8est_partition (p8est, 1, NULL);

        // partition time for i-step
        if(p8est->mpirank == 0)
        {
            snprintf (filename, 17, "partitio_%02d.time", i);
            file = fopen(filename, "w");
            SC_CHECK_ABORT(file, "Can't write to file");
            fprintf(file, "Step took %f seconds\n", ((double)clock() - start)/CLOCKS_PER_SEC);
            fclose(file);
        }
        //
    }

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

#include "ns.h"
#include "util.h"

void init(p8est_t *p8est,
          p4est_topidx_t which_tree,
          p8est_quadrant_t *quad) {
    context_t       *ctx = (context_t *) p8est->user_pointer;
    data_t          *data = (data_t *) quad->p.user_data;

    init_solver(data, ctx);
}

int refine_always (p8est_t *p8est,
                   p4est_topidx_t which_tree,
                   p8est_quadrant_t *q) {
    return 1;
}

int refine_fn (p8est_t *p8est,
           p4est_topidx_t which_tree,
           p8est_quadrant_t *q) {
    context_t *ctx = (context_t *) p8est->user_pointer;
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

void mesh_neighbors_iter(p8est_t *p8est,
                         p8est_ghost_t *ghost,
                         p8est_mesh_t *mesh,
                         void *ghost_data)
{
    p8est_quadrant_t    *neighbor  = NULL;
    p4est_topidx_t      which_tree = -1;
    context_t           *ctx       = (context_t *) p8est->user_pointer;

    p8est_mesh_face_neighbor_t mfn;
    p4est_locidx_t      qumid, quadrant_id, which_quad;
    p8est_quadrant_t    *q;                             /* current quad */
    data_t              *g_data;                        /* ghost data */
    data_t              *data;

    int                 nface;
    int                 nrank;
    double              i = 0;
    double              n_data = 0;  /* neighbors data */
    double              h;
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

            // TODO working with neighbours
        }
    }
}

void mesh_iter(p8est_t *p8est, p8est_mesh_t *mesh)
{
    p4est_topidx_t      which_tree = -1;
    context_t           *ctx       = (context_t *) p8est->user_pointer;

    p8est_mesh_face_neighbor_t mfn;
    p4est_locidx_t      qumid, quadrant_id;
    p8est_quadrant_t    *q;                             /* current quad */
    data_t              *data;

    double              h;
    double              p[3];

    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid)
    {
        q = p8est_mesh_quadrant_cumulative (p8est, qumid,
                                            &which_tree, &quadrant_id);

        /* side length */
        h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;
        get_midpoint(p8est, which_tree, q, p);

        data = (data_t *) q->p.user_data;

        // calculate CFL
        cflq(data, ctx, h);
    }
}

void solve(p8est_t *p8est,
           int step) {
    SC_PRODUCTIONF("Start solve step %d\n", step);

    int                 mpiret;
    data_t              *ghost_data;
    context_t           *ctx = (context_t *) p8est->user_pointer;
    p8est_ghost_t       *ghost;
    p8est_mesh_t        *mesh;

    /* create the ghost quadrants */
    ghost = p8est_ghost_new (p8est, P8EST_CONNECT_FACE);
    ghost_data = P4EST_ALLOC (data_t, ghost->ghosts.elem_count);
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);

    mesh = p8est_mesh_new(p8est, ghost, P8EST_CONNECT_FACE);
    SC_PRODUCTIONF("Used memory: %ld\n", p8est_mesh_memory_used(mesh));

    /* calc f2 and partials */
    SC_PRODUCTION("Cell iter started\n");
    mesh_iter(p8est, mesh);
    SC_PRODUCTION("Cell iter ended\n");

    /* calc min dt */
    SC_PRODUCTIONF("dt old: %.20lf\n", ctx->dt);
    mpiret = MPI_Allreduce(MPI_IN_PLACE, &ctx->dt, 1, MPI_DOUBLE, MPI_MIN, p8est->mpicomm);
    SC_CHECK_MPI(mpiret);
    SC_PRODUCTIONF("dt new: %.20lf\n", ctx->dt);

    /* exchange ghost data */
    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");
    /* calc f1 */
    SC_PRODUCTION("Neighbors iter started\n");
    mesh_neighbors_iter(p8est, ghost, mesh, ghost_data);
    SC_PRODUCTION("Neighbors iter ended\n");
    /* exchange ghost data */
    SC_PRODUCTION("Exchange started\n");
    p8est_ghost_exchange_data (p8est, ghost, ghost_data);
    SC_PRODUCTION("Exchange ended\n");

    /* generate vtk and print solution */
    //write_solution(p8est, step);

    /* clear */
    p8est_mesh_destroy(mesh);
    P4EST_FREE (ghost_data);
    p4est_ghost_destroy (ghost);

    SC_PRODUCTIONF("End solve step %d\n", step);
}

int main (int argc, char **argv){
    int                     i;
    char                    filename[BUFSIZ] = { '\0' };
    int                     mpiret;
    sc_MPI_Comm             mpicomm;
    context_t               ctx;
    p8est_t                 *p8est;
    p8est_connectivity_t    *conn;
    clock_t                 start;
    FILE                    *file;

    start = clock();

    ctx.center[0] = 0.5;
    ctx.center[1] = 0.5;
    ctx.center[2] = 0.5;

    ctx.width = 0.2;
    ctx.level = 1;
    ctx.dt = 0;

    ctx.Adiabatic = 1.4; //air

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
        solve(p8est, i);

        // solve time for i-step
        if(p8est->mpirank == 0)
        {
            snprintf (filename, 17, "solution_%02d.time", i);
            file = fopen(filename, "w");
            SC_CHECK_ABORT(file, "Can't write to file");
            fprintf(file, "Solution took %f seconds\n", ((double)clock() - start)/CLOCKS_PER_SEC);
            fclose(file);
        }

        if (i == ctx.level - 1)
            break;

        // repartition
        p8est_refine(p8est, 0, refine_always, init);
        p8est_partition (p8est, 1, NULL);
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
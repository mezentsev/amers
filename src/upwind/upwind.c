#include "upwind.h"

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

void solve(p8est_t *p8est) {
    int                     i;
    clock_t                 start;
    FILE                    *file;
    char                    filename[BUFSIZ] = { '\0' };
    context_t               *ctx = (context_t *) p8est->user_pointer;

    for (i = 0; i < ctx->steps; ++i)
    {
        start = clock();
        solver_step(p8est, i);

        /* вывод в файл времени вычисление N-ого решения */
        if(p8est->mpirank == 0)
        {
            snprintf (filename, 17, "solution_%02d.time", i);
            file = fopen(filename, "w");
            SC_CHECK_ABORT(file, "Can't write to file");
            fprintf(file, "Solution took %f seconds\n", ((double)clock() - start)/CLOCKS_PER_SEC);
            fclose(file);
        }

        if (i == ctx->steps - 1)
            break;

        // TODO не нужно делить каждый шаг. Потом нужно научить делить только там где нужно
        //p8est_refine(p8est, 1, refine_fn, init);
    }
}

int main (int argc, char **argv) {
    int                     mpiret;
    sc_MPI_Comm             mpicomm;
    context_t               ctx;
    p8est_t                 *p8est;
    p8est_connectivity_t    *conn;

    ctx.center[0] = 0.7;
    ctx.center[1] = 0.7;
    ctx.center[2] = 0.7;

    //ctx.v[0] = 0.485191768970225;
    //ctx.v[1] = -0.427996381877778;
    //ctx.v[2] = 0.762501176669961;
    ctx.v[0] = -0.7;
    ctx.v[1] = -0.7;
    ctx.v[2] = -0.7;

    ctx.width = 0.1;
    ctx.steps = 100;
    ctx.dt = 0;
    ctx.get_boundary_data_by_face = get_boundary_data_by_face;

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

    conn = p8est_connectivity_new_unitcube();

    p8est = p8est_new_ext (mpicomm, /* communicator */
                           conn,    /* connectivity */
                           0,       /* minimum quadrants per MPI process */
                           4,       /* minimum level of refinement */
                           1,       /* fill uniform */
                           sizeof (element_data_t),         /* data size */
                           init_solver,  /* initializes data */
                           (void *) (&ctx));              /* context */

    p8est_refine(p8est, 1, refine_fn, init_solver);
    //p8est_coarsen(p8est, 1, coarsen_fn, init);

    p8est_balance(p8est, P8EST_CONNECT_FULL, init_solver);
    p8est_partition (p8est, 1, NULL);

    // SOLVE
    solve(p8est);

    /* очистка всех данных и проверка */
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn);

    /* проверка на закрытие всех дескрипторов */
    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}
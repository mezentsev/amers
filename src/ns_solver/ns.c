//
// Created by Vadim Mezentsev on 01.01.17.
//

#include "ns.h"
#include "data.h"

void init(p8est_t *p8est,
          p4est_topidx_t which_tree,
          p8est_quadrant_t *quad) {
    ctx_t       *ctx = (ctx_t *) p8est->user_pointer;
    data_t      *data = (data_t *) quad->p.user_data;
}

void get_midpoint(p8est_t *p8est,
                  p4est_topidx_t tree,
                  p8est_quadrant_t *q,
                  double xyz[3]) {
    p4est_qcoord_t half_length = P8EST_QUADRANT_LEN(q->level) / 2;

    p8est_qcoord_to_vertex (p8est->connectivity,
                            tree,
                            q->x + half_length,
                            q->y + half_length,
                            q->z + half_length,
                            xyz);
}

int refine_always (p8est_t *p8est,
                   p4est_topidx_t which_tree,
                   p8est_quadrant_t *q) {
    return 1;
}

int refine_fn (p8est_t *p8est,
           p4est_topidx_t which_tree,
           p8est_quadrant_t *q) {
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

void solve(p8est_t *p8est,
           int step) {
    SC_PRODUCTIONF("Start solve step %d\n", step);

    SC_PRODUCTIONF("End solve step %d\n", step);
}

int main (int argc, char **argv){
    int                     i;
    char                    filename[BUFSIZ] = { '\0' };
    int                     mpiret;
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
    ctx.level = 4;

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
#include <p4est_to_p8est.h>
#include <p8est_connectivity.h>
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

typedef struct data
{
    double              u;             // the state variable
} data_t;

// Describes current problem
typedef struct ctx
{
    double              center[P4EST_DIM];  // coordinates of the center
    double              width;
} ctx_t;

/// Compute the value of the initial condition.
static double
init_condition(double x[3], ctx_t *ctx)
{
    int                 i;
    double              *c = ctx->center;
    double              bump_width = ctx->width;
    double              r2, d[P4EST_DIM];
    double              arg, retval;

    r2 = 0.;
    for (i = 0; i < P4EST_DIM; i++) {
        d[i] = x[i] - c[i];
        r2 += d[i] * d[i];
    }

    arg = -(1. / 2.) * r2 / bump_width / bump_width;
    retval = exp (arg);

    return retval;
}

/// Get the coordinates of the midpoint of a quadrant.
static void
get_midpoint(p8est_t *p4est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3])
{
    p4est_qcoord_t half_length = P4EST_QUADRANT_LEN(q->level) / 2;

    p4est_qcoord_to_vertex (p4est->connectivity,
                            tree,
                            q->x + half_length,
                            q->y + half_length,
                            q->z + half_length,
                            xyz);
}

static void
init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q)
{
    /* the data associated with a forest is accessible by user_pointer */
    ctx_t *ctx = (ctx_t *) p4est->user_pointer;
    /* the data associated with a quadrant is accessible by p.user_data */
    data_t *data = (data_t *) q->p.user_data;
    double midpoint[3];

    get_midpoint(p4est, which_tree, q, midpoint);
    /* initialize the data */
    data->u = init_condition(midpoint, ctx);
}

/// http://p4est.github.io/api/p8est_8h.html#aeb61645ae5dbbdbcb3c4f8a8810b4ecf
static int
refine_fn (p8est_t *p8est,
           p4est_topidx_t which_tree,
           p8est_quadrant_t *q)
{
    ctx_t *ctx = (ctx_t *) p8est->user_pointer;
    double midpoint[3];
    double h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;

    get_midpoint(p8est, which_tree, q, midpoint);

    if (pow(midpoint[0] - ctx->center[0], 2.) +
        pow(midpoint[1] - ctx->center[1], 2.) +
        pow(midpoint[2] - ctx->center[2], 2.) <= ctx->width &&
            q->level < 7)
        return 1;

    return 0;
}

static int
coarsen_fn (p8est_t *p8est,
            p4est_topidx_t which_tree,
            p8est_quadrant_t **children)
{
    return 0;
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

    /* Initialize MPI; see sc_mpi.h.
     * If configure --enable-mpi is given these are true MPI calls.
     * Else these are dummy functions that simulate a single-processor run. */
    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    /* These functions are optional.  If called they store the MPI rank as a
     * static variable so subsequent global p4est log messages are only issued
     * from processor zero.  Here we turn off most of the logging; see sc.h. */
    sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
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

    p8est_vtk_write_file (p8est, NULL, "amr_test_vtk");

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

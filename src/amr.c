#include <p4est_to_p8est.h>
#include <p8est_connectivity.h>
#include <p8est.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>

typedef struct data
{
    double              u;             // the state variable
    double              du[P4EST_DIM]; // the spatial derivatives
    double              dudt;          // the time derivative
} data_t;

// Describes current problem
typedef struct ctx
{
    double              center[P4EST_DIM];  // coordinates of the center
    double              width;
} ctx_t;

/// Compute the value and derivatives of the initial condition.
static double init_condition(double x[3], double du[3], ctx_t *ctx)
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

    if (du) {
        for (i = 0; i < P4EST_DIM; i++) {
            du[i] = -(1. / bump_width / bump_width) * d[i] * retval;
        }
    }

    return retval;
}

/// Get the coordinates of the midpoint of a quadrant.
static void get_midpoint(p8est_t *p4est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3])
{
    p4est_qcoord_t half_length = P4EST_QUADRANT_LEN(q->level) / 2;

    p4est_qcoord_to_vertex (p4est->connectivity,
                            tree,
                            q->x + half_length,
                            q->y + half_length,
                            q->z + half_length,
                            xyz);
}

static void write_solution(p4est_t *p4est)
{
    char                filename[BUFSIZ] = { "amr\0" };
    sc_array_t          *u_interp;
    p4est_locidx_t      numquads;

    //snprintf (filename, 8, "ex_%04d", timestep);

    numquads = p4est->local_num_quadrants;

    /* create a vector with one value for the corner of every local quadrant
     * (the number of children is always the same as the number of corners) */
    u_interp = sc_array_new_size (sizeof (double), numquads * P4EST_CHILDREN);

    /* Use the iterator to visit every cell and fill in the solution values.
     * Using the iterator is not absolutely necessary in this case: we could
     * also loop over every tree (there is only one tree in this case) and loop
     * over every quadrant within every tree, but we are trying to demonstrate
     * the usage of p4est_iterate in this example */
    p4est_iterate(p4est, NULL,   /* we don't need any ghost quadrants for this loop */
                  (void *) u_interp,     /* pass in u_interp so that we can fill it */
                  //interpolate_solution,    /* callback function that interpolates from the cell center to the cell corners, defined above */
                  NULL,
                  NULL,          /* there is no callback for the faces between quadrants */
                  NULL,          /* there is no callback for the edges between quadrants */
                  NULL);         /* there is no callback for the corners between quadrants */

    /* create VTK output context and set its parameters */
    p4est_vtk_context_t *context = p4est_vtk_context_new (p4est, filename);
    p4est_vtk_context_set_scale (context, 0.99);  /* quadrant at almost full scale */

    /* begin writing the output files */
    context = p4est_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing vtk header");

    /* do not write the tree id's of each quadrant
     * (there is only one tree in this example) */
    context = p4est_vtk_write_cell_dataf (context, 0, 1,  /* do write the refinement level of each quadrant */
                                          1,      /* do write the mpi process id of each quadrant */
                                          0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                          0,      /* there is no custom cell scalar data. */
                                          0,      /* there is no custom cell vector data. */
                                          context);       /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    /* write one scalar field: the solution value */
    context = p4est_vtk_write_point_dataf (context, 1, 0, /* write no vector fields */
                                           "solution", u_interp, context);        /* mark the end of the variable cell data. */
    SC_CHECK_ABORT (context != NULL,
                    P4EST_STRING "_vtk: Error writing cell data");

    const int           retval = p4est_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

    sc_array_destroy (u_interp);
}

static void init(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q)
{
    /* the data associated with a forest is accessible by user_pointer */
    ctx_t *ctx = (ctx_t *) p4est->user_pointer;
    /* the data associated with a quadrant is accessible by p.user_data */
    data_t *data = (data_t *) q->p.user_data;
    double midpoint[3];

    get_midpoint(p4est, which_tree, q, midpoint);
    /* initialize the data */
    data->u = init_condition (midpoint, data->du, ctx);
}

int
main (int argc, char **argv)
{
    int                     mpiret;
    sc_MPI_Comm             mpicomm;
    ctx_t                   ctx;
    p4est_t                 *p4est;
    p4est_connectivity_t    *conn;

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


    p4est = p4est_new_ext (mpicomm, /* communicator */
                           conn,    /* connectivity */
                           0,       /* minimum quadrants per MPI process */
                           4,       /* minimum level of refinement */
                           1,       /* fill uniform */
                           sizeof (data_t),         /* data size */
                           init,  /* initializes data */
                           (void *) (&ctx));              /* context */

    //p4est_refine(p4est, 0, NULL, init);
    //p8est_coarsen(p8est, 0, NULL, init);

    //p4est_balance(p4est, P4EST_CONNECT_FACE, init);
    //p4est_partition (p4est, 0, NULL);

    write_solution(p4est);

    // clear
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);

    /* Verify that allocations internal to p4est and sc do not leak memory.
     * This should be called if sc_init () has been called earlier. */
    sc_finalize ();

    /* This is standard MPI programs.  Without --enable-mpi, this is a dummy. */
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
}

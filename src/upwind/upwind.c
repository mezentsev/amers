#include "upwind.h"
#include "../util.h"

void init(p8est_t *p8est,
          p4est_topidx_t which_tree,
          p8est_quadrant_t *quad) {
    context_t       *ctx = (context_t *) p8est->user_pointer;
    element_data_t  *data = (element_data_t *) quad->p.user_data;

    init_solver(p8est, data);
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
        p8est_partition (p8est, 1, NULL);
    }
}

void get_quad_data(p8est_iter_volume_info_t *info,
                   void *user_data) {
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
    element_data_t      *data = (element_data_t *) q->p.user_data;
    p4est_locidx_t      array_offset;
    double              *this_u_ptr;
    int                 i;

    tree = p8est_tree_array_index (p8est->trees, which_tree);
    local_id += tree->quadrants_offset;             /* now the id is relative to the MPI process */
    array_offset = P8EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */

    for (i = 0; i < P8EST_CHILDREN; i++) {
        this_u_ptr = (double *) sc_array_index (u_interp, array_offset + i);
        this_u_ptr[0] = data->dummy;
    }
}

void write_solution(p8est_t *p8est,
                    int step) {
    char                filename[BUFSIZ] = { '\0' };
    p4est_locidx_t      numquads;
    sc_array_t          *data;

    snprintf (filename, 12, "solution_%04d", step);
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
                  get_quad_data,
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
                                          "boundary",
                                          data,
                                          context);
    SC_CHECK_ABORT (context != NULL,
                    P8EST_STRING "_vtk: Error writing cell data");

    const int retval = p8est_vtk_write_footer(context);
    SC_CHECK_ABORT (!retval,
                    P8EST_STRING "_vtk: Error writing footer");

    sc_array_destroy(data);
}

int main (int argc, char **argv) {
    int                     mpiret;
    sc_MPI_Comm             mpicomm;
    context_t               ctx;
    p8est_t                 *p8est;
    p8est_connectivity_t    *conn;

    ctx.center[0] = 0.3;
    ctx.center[1] = 0.3;
    ctx.center[2] = 0.3;

    ctx.width = 0.1;
    ctx.steps = 20;
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
                           1,       /* minimum level of refinement */
                           1,       /* fill uniform */
                           sizeof (element_data_t),         /* data size */
                           init,  /* initializes data */
                           (void *) (&ctx));              /* context */

    // TODO сначала проверяем работу алгоритма на регулярной сетке, потом включаем адаптацию
    p8est_refine(p8est, 1, refine_fn, init);
    //p8est_coarsen(p8est, 1, coarsen_fn, init);

    p8est_balance(p8est, P8EST_CONNECT_FULL, init);
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
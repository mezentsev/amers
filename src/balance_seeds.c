//
// Created by Vadim Mezentsev on 18.12.16.
//

#include <p4est_to_p8est.h>
#include <sc_dmatrix.h>
#include <p8est_balance.h>
#include <p8est_bits.h>
#include <p8est_vtk.h>
#include <p8est_extended.h>

typedef struct
{
    int                 flag;
} balance_seeds_elem_t;

static const int refine_level = 8;
static p8est_quadrant_t center =
  {
          0x30000, // x
          0x30000, // y
          0x30000, // z
          4,       // refine level
          0,       // pad8
          0,       // pad16
          {NULL}   // p8est_quadrant_data
  };

static void
init_fn (p8est_t * p4est, p4est_topidx_t which_tree,
         p8est_quadrant_t * quadrant)
{
    ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -2;
}

static int
refine_fn (p8est_t * p4est, p4est_topidx_t which_tree,
           p8est_quadrant_t * quadrant)
{
    p8est_connect_type_t balance = P8EST_CONNECT_EDGE;
    p8est_quadrant_t    desc;
    int                 i;

    if (((balance_seeds_elem_t *) (quadrant->p.user_data))->flag > -2) {
        return 0;
    }

    if (p8est_quadrant_is_ancestor (quadrant, &center)) {
        return 1;
    }

    if (p8est_quadrant_is_equal (quadrant, &center)) {
        ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = center.level;
        return 0;
    }

    /*if (quadrant->x >= 0x60000 || quadrant->y >= 0x60000 || quadrant->z >= 0x60000) {
      ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -1;
      return 0;
    }*/

    for (i = 0; i < P8EST_CHILDREN; i++) {
        // check potomok
        p8est_quadrant_corner_descendant (quadrant, &desc, i, P8EST_QMAXLEVEL);
        // balance condition
        if (p8est_balance_seeds (&desc, &center, balance, NULL)) {
          break;
        }
    }
    if (i == P8EST_CHILDREN) {
        P4EST_ASSERT (!p8est_balance_seeds (quadrant, &center, balance, NULL));
        ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = -1;
        return 0;
    }

    // potomok
    p8est_quadrant_corner_descendant (quadrant, &desc, i, quadrant->level + 1);
    if (!p8est_balance_seeds (&desc, &center, balance, NULL)) {
        if (quadrant->level < refine_level) {
            return 1;
        }
    }

    ((balance_seeds_elem_t *) (quadrant->p.user_data))->flag = quadrant->level;
    return 0;
}

int
main (int argc, char **argv)
    {
    sc_MPI_Comm         mpicomm;
    int                 mpiret;
    int                 mpisize, mpirank;
    p8est_t            *p4est;
    p8est_connectivity_t *connectivity;
    sc_dmatrix_t       *vtkvec;
    p8est_tree_t       *tree;
    sc_array_t         *quadrants;
    size_t              zz, count;
    p8est_quadrant_t   *q;
    int                 i;

    char                filename[] = "p8est_balance_edge";

    /* initialize MPI */
    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;
    mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
    SC_CHECK_MPI (mpiret);

    sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
    p4est_init (NULL, SC_LP_DEFAULT);

    /* Create a forest from cube */
    connectivity = p8est_connectivity_new_unitcube ();

    p4est = p8est_new_ext (mpicomm,
                         connectivity,
                         0, // min_quadrants per processor
                         2, // Forest is refined at least to this level
                         1, // Fill the forest with a uniform mesh
                            // instead of the coarsest possible one.
                         sizeof (balance_seeds_elem_t),
                         init_fn,
                         NULL);
    /* After create */
    p8est_vtk_write_file (p4est, NULL, "cube_created");

    p8est_refine (p4est,
                1,         // Recursive refinement
                refine_fn, // Callback function that must return true if a quadrant
                           // shall be refined. If refine_recursive is true,
                           // refine_fn is called for every existing and newly created quadrant
                init_fn);  // Callback function to initialize the user_data of newly
                           // created quadrants, which is already allocated
    /* After refine */
    p8est_vtk_write_file (p4est, NULL, "cube_refined");

    /* Test vtk part writer */
    p8est_vtk_context_t *context = p8est_vtk_context_new (p4est, filename);
    p8est_vtk_context_set_scale (context, 1. - 2. * SC_EPS);
    context = p8est_vtk_write_header (context);
    SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing header");

    vtkvec = sc_dmatrix_new (p4est->local_num_quadrants, P8EST_CHILDREN);
    tree = p8est_tree_array_index (p4est->trees, 0);
    quadrants = &(tree->quadrants);
    count = quadrants->elem_count;

    for (zz = 0; zz < count; zz++) {
        q = p8est_quadrant_array_index (quadrants, zz);
        for (i = 0; i < P8EST_CHILDREN; i++) {
          vtkvec->e[zz][i] = (double) ((balance_seeds_elem_t *) (q->p.user_data))->flag;
        }
    }
    sc_array_t *level =
    sc_array_new_data ((void *) vtkvec->e[0], sizeof (double),
                       count * P8EST_CHILDREN);
    context = p8est_vtk_write_point_dataf (context, 1, 0, "level", level, context);
    SC_CHECK_ABORT (context != NULL,
                  P8EST_STRING "_vtk: Error writing point data");
    sc_array_destroy (level);

    const int retval = p8est_vtk_write_footer (context);
    SC_CHECK_ABORT (!retval, P8EST_STRING "_vtk: Error writing footer");

    /* Balance tree */
    p4est_balance (p4est, P8EST_CONNECT_EDGE, NULL);
    /* Redistribute */
    p4est_partition (p4est, 0, NULL);

    /* Write whole vtk file */
    p8est_vtk_write_file (p4est, NULL, "cube_finally");

    sc_dmatrix_destroy (vtkvec);
    p8est_destroy (p4est);
    p8est_connectivity_destroy (connectivity);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}

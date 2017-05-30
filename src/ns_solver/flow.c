#include "solver.h"
#include "util.h"

void calc_flow_mesh_iter(p8est_t *p8est,
                         p8est_mesh_t *mesh,
                         p8est_ghost_t *ghost,
                         void *ghost_data) {
    int                         qumid;
    int                         which_tree = -1;
    int                         quadrant_id;
    p8est_quadrant_t            *q;
    context_t                   *ctx       = (context_t *) p8est->user_pointer;

    element_data_t              *ndata;
    p8est_quadrant_t            *qn;
    p8est_mesh_face_neighbor_t  mfn;
    int                         nface, nrank;
    int                         nquadrant_id;
    int                         nwhich_tree = -1;

    int                         is_boundary = 0;
    element_data_t              boundary_data;


    SC_PRODUCTION("*****************\n");

    // прохождение по всем локальным ячейкам
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        // получение текущей ячейки
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);
        // инииализация итератора по соседям
        p8est_mesh_face_neighbor_init2(&mfn, p8est, ghost, mesh, which_tree, quadrant_id);

        /**
         * проход по соседям
         * если qn == q, то nface - граничный face
         */
        while ((qn = p8est_mesh_face_neighbor_next(&mfn,
                                                   &nwhich_tree,
                                                   &nquadrant_id,
                                                   &nface,
                                                   &nrank)) != NULL) {
            // если ячейки совпадают - попали на границу
            is_boundary = p8est_quadrant_is_equal(qn, q);

            // извлечение данных из теневого слоя
            ndata = (element_data_t *) p8est_mesh_face_neighbor_data(&mfn, ghost_data);

            if (is_boundary) {
                SC_PRODUCTIONF("[BOUNDARY] Data: %lf, nface: %d, nrank: %d, cur rank: %d\n",
                               ndata->Density, nface, nrank, p8est->mpirank);

                boundary_data = get_boundary_data_by_face(qn, nface);
                SC_PRODUCTIONF("New data from boundary: %lf\n", boundary_data.Pressure);
            } else {
                SC_PRODUCTIONF("Data: %lf, nface: %d, nrank: %d, cur rank: %d\n",
                               ndata->Density, nface, nrank, p8est->mpirank);
            }
        }

        SC_PRODUCTION("*****************\n");
    }
}

__deprecated
void calc_flow_face_iter(p8est_iter_face_info_t *info,
                         void *user_data) {
    p8est_quadrant_t            *quad = NULL;
    p8est_iter_face_side_t      *side = NULL;
    int                         which_face;
    int                         is_boundary = (info->tree_boundary == P8EST_CONNECT_FACE);

    SC_PRODUCTION("*****************\n");
    SC_PRODUCTIONF("Is boundary: %d\n", is_boundary);
    SC_PRODUCTIONF("Sides count: %d\n", info->sides.elem_count);
    SC_PRODUCTIONF("Orientation: %d\n", info->orientation);

    /**
     * Вычисление через граничные функции, хранящиеся в контексте.
     * Либо вычисление потока по обычной схеме
     **/
    if (is_boundary) {
        /* which of the quadrant's faces the interface touches */
        side = p8est_iter_fside_array_index(&info->sides, 0);
        which_face = side->face;
        SC_PRODUCTIONF("[BOUNDARY] Side is: %d; treeid: %d\n", side->face, side->treeid);

        switch (which_face) {
            case 0:                      /* -x side */
                SC_PRODUCTION("-x side\n");
                //vdotn = -ctx->v[0];
                break;
            case 1:                      /* +x side */
                SC_PRODUCTION("+x side\n");
                //vdotn = ctx->v[0];
                break;
            case 2:                      /* -y side */
                SC_PRODUCTION("-y side\n");
                //vdotn = -ctx->v[1];
                break;
            case 3:                      /* +y side */
                SC_PRODUCTION("+y side\n");
                //vdotn = ctx->v[1];
                break;
            case 4:                      /* -z side */
                SC_PRODUCTION("-z side\n");
                //vdotn = -ctx->v[2];
                break;
            case 5:                      /* +z side */
                SC_PRODUCTION("+z side\n");
                //vdotn = ctx->v[2];
                break;
            default:
                SC_ABORT("Wrong face");
                // TODO вычислить граничное значение
        }
    } else {
        // Это не граница, есть обе стороны
        // TODO почитать поток
        for (size_t i = 0; i < info->sides.elem_count; ++i) {
            side = p8est_iter_fside_array_index(&info->sides, i);
            SC_PRODUCTIONF("Side is: %d; treeid: %d\n", side->face, side->treeid);

            if (!side->is_hanging) {
                SC_PRODUCTION("Quad is full\n");
                quad = side->is.full.quad;
                quadrant_pprint(quad, side->is.full.is_ghost);
            } else {
                SC_PRODUCTION("Quad is hanging\n");
                for (int j = 0; j < P8EST_HALF; ++j) {
                    quad = side->is.hanging.quad[j];
                    quadrant_pprint(quad, side->is.hanging.is_ghost[j]);
                }
            }
        }
    }

    SC_PRODUCTION("*****************\n");
}

__deprecated
void calc_flow_volume_iter(p8est_iter_volume_info_t *info,
                           void *user_data) {
    int                 i;
    p8est_quadrant_t    nq[P8EST_FACES];                // neighbor quads
    p8est_quadrant_t    *q = info->quad;
    p8est_t             *p8est = info->p4est;
    p8est_ghost_t       *ghost = info->ghost_layer;

    p4est_topidx_t      which_tree = -1;
    context_t           *ctx       = (context_t *) p8est->user_pointer;

    int                 is_small_neighbors_exists;
    int                 is_same_neighbor_exists;
    int                 is_big_neighbor_exists;

    int                 is_boundary;

    SC_PRODUCTION("*****************\n");

    for (i = 0; i < P8EST_FACES; ++i) {
        p8est_quadrant_all_face_neighbors(q, i, nq);

        is_small_neighbors_exists = p8est_quadrant_exists(p8est, ghost, info->treeid, &nq[0], NULL, NULL, NULL);
        is_same_neighbor_exists = p8est_quadrant_exists(p8est, ghost, info->treeid, &nq[P8EST_HALF], NULL, NULL, NULL);
        is_big_neighbor_exists = p8est_quadrant_exists(p8est, ghost, info->treeid, &nq[P8EST_HALF + 1], NULL, NULL, NULL);

        if (is_same_neighbor_exists) {
            element_data_t *data = (element_data_t *) nq[P8EST_HALF].p.user_data;

            SC_PRODUCTIONF("Quad is (%d), neib (%d), value (%d)\n", q->p.which_tree, nq[P8EST_HALF].p.which_tree, data->Pressure);
        }

        is_boundary = !is_small_neighbors_exists && !is_same_neighbor_exists && !is_big_neighbor_exists;

        /**
         * Вычисление через граничные функции, хранящиеся в контексте.
         * Либо вычисление потока по обычной схеме
         **/
        if (is_boundary) {
            switch (i) {
                case 0:                      /* -x side */
                    SC_PRODUCTION("-x side\n");
                    //vdotn = -ctx->v[0];
                    break;
                case 1:                      /* +x side */
                    SC_PRODUCTION("+x side\n");
                    //vdotn = ctx->v[0];
                    break;
                case 2:                      /* -y side */
                    SC_PRODUCTION("-y side\n");
                    //vdotn = -ctx->v[1];
                    break;
                case 3:                      /* +y side */
                    SC_PRODUCTION("+y side\n");
                    //vdotn = ctx->v[1];
                    break;
                case 4:                      /* -z side */
                    SC_PRODUCTION("-z side\n");
                    //vdotn = -ctx->v[2];
                    break;
                case 5:                      /* +z side */
                    SC_PRODUCTION("+z side\n");
                    //vdotn = ctx->v[2];
                    break;
                default:
                    SC_ABORT("Wrong face");

                    // TODO посчитать граничное значение
            }
        } else {
            // TODO вычислить поток
            //SC_PRODUCTIONF();
        }
    }

    SC_PRODUCTION("*****************\n");
}



/*void mesh_neighbors_iter(p8est_t *p8est,
                         p8est_ghost_t *ghost,
                         p8est_mesh_t *mesh,
                         void *ghost_data) {
    p4est_topidx_t      which_tree = -1;
    context_t           *ctx       = (context_t *) p8est->user_pointer;

    p8est_mesh_face_neighbor_t mfn;
    p4est_locidx_t      qumid, quadrant_id;

    p8est_quadrant_t    *q;                             // current quad
    p8est_quadrant_t    nq[P8EST_FACES];                // neighbor quads

    element_data_t              *data;

    int                 nface;
    double              nx, ny, nz;

    int                 small_neighbors;
    int                 same_neighbor;
    int                 big_neighbor;

    // iterate for all quadrants
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        q = p8est_mesh_quadrant_cumulative (p8est,
                                            qumid,
                                            &which_tree,
                                            &quadrant_id);
        data = (element_data_t *) q->p.user_data;

        // loop for all faces
        for (nface = 0; nface < P8EST_FACES; ++nface){
            // neighbor orientation
            nx = 0;
            ny = 0;
            nz = 0;

            switch (nface) {
                case 0:                      // -x side
                    nx = -1;
                    break;
                case 1:                      // +x side
                    nx = 1;
                    break;
                case 2:                      // -y side
                    ny = -1;
                    break;
                case 3:                      // +y side
                    ny = 1;
                    break;
                case 4:                      // -z side
                    nz = -1;
                    break;
                case 5:                      // +z side
                    nz = 1;
                    break;
                default:
                    break;
            }

            p8est_quadrant_all_face_neighbors(q, nface, nq);

            small_neighbors = p8est_quadrant_exists(p8est, ghost, which_tree, &nq[0], NULL, NULL, NULL);
            same_neighbor = p8est_quadrant_exists(p8est, ghost, which_tree, &nq[P8EST_HALF], NULL, NULL, NULL);
            big_neighbor = p8est_quadrant_exists(p8est, ghost, which_tree, &nq[P8EST_HALF + 1], NULL, NULL, NULL);

            // check boundary
            if (!small_neighbors && !same_neighbor && !big_neighbor) {
                //p8est_quadrant_print(SC_LP_PRODUCTION, q);
                if (data->boundary == -1)
                    data->boundary = ipow(2, nface);
                else
                    data->boundary |= ipow(2, nface);

                data->bx = nx;
                data->by = ny;
                data->bz = nz;
            }

            //h = (double) P8EST_QUADRANT_LEN (q->level) / (double) P8EST_ROOT_LEN;

            // TODO Calculate flow
            flow(data, nx, ny, nz);
        }
    }
}*/


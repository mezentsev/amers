#include <p8est_mesh.h>
#include "solver.h"
#include "util.h"

void calc_flux_mesh_iter(p8est_t *p8est,
                         p8est_mesh_t *mesh,
                         p8est_ghost_t *ghost,
                         void *ghost_data) {
    int                         qumid;
    int                         which_tree = -1;
    int                         quadrant_id;
    p8est_quadrant_t            *q;
    element_data_t              *data;
    context_t                   *ctx       = (context_t *) p8est->user_pointer;

    element_data_t              *ndata;
    p8est_quadrant_t            *qn;

    p8est_mesh_face_neighbor_t  mfn;
    int                         nface, nrank;
    int                         nquadrant_id;
    int                         nwhich_tree = -1;
    double                      nh; // сторона соседа

    int                         is_boundary = 0;
    double                      V; // объём текущей ячейки
    double                      Sn; // площадь между соседом и текущей ячейкой
    double                      h; // сторона текущей ячейки

    p8est_quadrant_t            quad_n[P8EST_FACES];                // neighbor quads
    int                         i = 0;
    int                         is_small_neighbors_exists;
    int                         is_same_neighbor_exists;
    int                         is_big_neighbor_exists;

    element_data_t              new_data; // посчитанный поток в локальном базисе
    element_data_t              sum_flux; // суммируем поток со всех соседей

    SC_PRODUCTION("*****************\n");

    // прохождение по всем локальным ячейкам
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        // получение текущей ячейки
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);

        h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;
        V = h*h*h;

        // даные текущей ячейки
        data = (element_data_t *) q->p.user_data;
        init_empty_solver(p8est, &sum_flux); // инициализируем пустыми значениями для каждой ячейки

        // TODO проверка на существование по всем фейсам
        for (i = 0; i < P8EST_FACES; ++i) {
            p8est_quadrant_all_face_neighbors(q, i, quad_n);

            is_small_neighbors_exists = p8est_quadrant_exists(p8est, ghost, which_tree, &quad_n[0], NULL, NULL, NULL);
            is_same_neighbor_exists = p8est_quadrant_exists(p8est, ghost, which_tree, &quad_n[P8EST_HALF], NULL, NULL, NULL);
            is_big_neighbor_exists = p8est_quadrant_exists(p8est, ghost, which_tree, &quad_n[P8EST_HALF + 1], NULL, NULL, NULL);

            if (is_same_neighbor_exists) {
                element_data_t *d = (element_data_t *) quad_n[P8EST_HALF].p.user_data;

                SC_PRODUCTIONF("Quad is (%d), neib (%d), value (%d)\n", q->p.which_tree, quad_n[P8EST_HALF].p.which_tree, d->Pressure);
            }
        }

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
            // длина стороны соседа
            nh = (double) P4EST_QUADRANT_LEN (qn->level) / (double) P4EST_ROOT_LEN;

            // площадь соприкосноверия
            Sn = h < nh
                    ? h * h
                    : nh * nh;

            // если ячейки совпадают - попали на границу
            is_boundary = p8est_quadrant_is_equal(qn, q);

            // извлечение данных из соседа
            ndata = (element_data_t *) p8est_mesh_face_neighbor_data(&mfn, ghost_data);

            // Граничные значения забираем в соответствии с функцией контекста
            if (is_boundary) {
                SC_PRODUCTIONF("[BOUNDARY] Data: %lf, nface: %d, nrank: %d, cur rank: %d\n",
                               ndata->dummy, nface, nrank, p8est->mpirank);

                new_data = get_boundary_data_by_face(p8est, qn, nface);
                SC_PRODUCTIONF("New data from boundary: %lf\n", new_data.dummy);

                // TODO test
                data->dummy += new_data.dummy;
            } else {
                SC_PRODUCTIONF("Data: %lf, nface: %d, nrank: %d, cur rank: %d, face: %d, subface: %d\n",
                               ndata->dummy, nface, nrank, p8est->mpirank, mfn.face, mfn.subface);

                // TODO проблема с face - почему не тот? откуда берётся 6?
                new_data = calc_flux(p8est, q, qn, nface);

                // TODO test
                // сверху вниз
                if (nface == FROM_TOP) {
                    data->dummy += ndata->dummy;
                }
            }

            // TODO в отдельную функцию
            sum_flux.Pressure  += new_data.Pressure * Sn;
            sum_flux.Density   += new_data.Density * Sn;
            sum_flux.u1        += new_data.u1 * Sn;
            sum_flux.u2        += new_data.u2 * Sn;
            sum_flux.u3        += new_data.u3 * Sn;
        }

        // TODO в отдельную функцию
        sum_flux.Pressure       = sum_flux.Pressure * ctx->dt / V;
        sum_flux.Density        = sum_flux.Density * ctx->dt / V;
        sum_flux.u1             = sum_flux.u1 * ctx->dt / V;
        sum_flux.u2             = sum_flux.u2 * ctx->dt / V;
        sum_flux.u3             = sum_flux.u3 * ctx->dt / V;

        data->Pressure         -= sum_flux.Pressure;
        data->Density          -= sum_flux.Density;
        data->u1               -= sum_flux.u1;
        data->u2               -= sum_flux.u2;
        data->u3               -= sum_flux.u3;

        SC_PRODUCTION("*****************\n");
    }
}

__deprecated
void calc_flux_face_iter(p8est_iter_face_info_t *info,
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
void calc_flux_volume_iter(p8est_iter_volume_info_t *info,
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

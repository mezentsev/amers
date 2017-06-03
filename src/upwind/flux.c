#include <p8est_mesh.h>
#include "solver.h"
#include "../util.h"

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
    int                         nface, nrank, neib_face;
    int                         nquadrant_id;
    int                         nwhich_tree = -1;
    double                      nh; // сторона соседа

    int                         is_boundary = 0;
    double                      V; // объём текущей ячейки
    double                      Sn; // площадь между соседом и текущей ячейкой
    double                      h; // сторона текущей ячейки

    int                         qtf;
    p4est_locidx_t              qtq, quadfacecode;
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
        which_tree = -1;
        neib_face = 0;
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
            neib_face = get_neighbour_face_by_next_one(mfn.face, mfn.subface);

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

            P4EST_ASSERT (p8est_quadrant_is_equal (qn, &(ndata->quad)));
            P4EST_ASSERT (ndata->quad.p.which_tree == which_tree);

            // Граничные значения забираем в соответствии с функцией контекста
            if (is_boundary) {
                SC_PRODUCTIONF("[BOUNDARY] Data: %lf, nface: %d, nrank: %d, cur rank: %d, neib_face: %d\n",
                               ndata->dummy, nface, nrank, p8est->mpirank, neib_face);

                new_data = get_boundary_data_by_face(p8est, qn, nface);
                SC_PRODUCTIONF("New data from boundary: %lf\n", new_data.dummy);

                // TODO test
                if (neib_face == 3 && data->dummy == 0) {
                    data->dummy += new_data.dummy;
                }
            } else {
                SC_PRODUCTIONF("Data: %lf, nface: %d, nrank: %d, cur rank: %d, neib_face: %d\n",
                               ndata->dummy, nface, nrank, p8est->mpirank, neib_face);

                new_data = calc_flux(p8est, q, qn, neib_face);

                // TODO test
                // перенос сверху вниз
                if (neib_face == 3) {
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
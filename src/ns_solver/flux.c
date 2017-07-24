#include <p8est_mesh.h>
#include "solver.h"
#include "../util.h"
#include "data.h"

element_data_t calc_flux(p8est_t *p8est, element_data_t *idata, element_data_t *ndata, int nface) {
    element_data_t  F;
    element_data_t  fi, fn;
    element_data_t  fsum;
    //element_data_t  c;

    context_t       *ctx    = (context_t *) p8est->user_pointer;
    double          sc      = calc_speed(idata->Z.Density, idata->Z.Pressure, ctx->Adiabatic);

    init_empty_solver(p8est, &F);
    init_empty_solver(p8est, &fsum);
    //init_solver_by_double(p8est, &c, sc);

    init_empty_solver(p8est, &fi);
    init_empty_solver(p8est, &fn);

    SC_INFOF("Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                   F.Z.Density, F.Z.u1, F.Z.u2, F.Z.u3, F.Z.Pressure);

    // Z + c
    //sumZ(p8est, &c, idata);

    /* направление face от соседа(-n) к текущей ячейке(-i) */
    switch (nface) {
        case 0:                      /* -x side */
        case 1:                      /* +x side */
            SC_LDEBUGF("face %d for -x and +x side\n", nface);
            fi.Z.Density    = idata->Z.Density * idata->Z.u1;
            fi.Z.u1         = idata->Z.Density * idata->Z.u1 * idata->Z.u1 + idata->Z.Pressure;
            fi.Z.u2         = idata->Z.Density * idata->Z.u1 * idata->Z.u2;
            fi.Z.u3         = idata->Z.Density * idata->Z.u1 * idata->Z.u3;
            fi.Z.Pressure   = idata->Z.Density * idata->Z.u1 * (idata->Z.E + idata->Z.Pressure / idata->Z.Density);

            fn.Z.Density    = ndata->Z.Density * ndata->Z.u1;
            fn.Z.u1         = ndata->Z.Density * ndata->Z.u1 * ndata->Z.u1 + ndata->Z.Pressure;
            fn.Z.u2         = ndata->Z.Density * ndata->Z.u1 * ndata->Z.u2;
            fn.Z.u3         = ndata->Z.Density * ndata->Z.u1 * ndata->Z.u3;
            fn.Z.Pressure   = ndata->Z.Density * ndata->Z.u1 * (ndata->Z.E + ndata->Z.Pressure / ndata->Z.Density);
            break;
        case 2:                      /* -y side */
        case 3:                      /* +y side */
            SC_LDEBUGF("face %d for -y and +y side\n", nface);
            fi.Z.Density    = idata->Z.Density * idata->Z.u2;
            fi.Z.u1         = idata->Z.Density * idata->Z.u2 * idata->Z.u1;
            fi.Z.u2         = idata->Z.Density * idata->Z.u2 * idata->Z.u2 + idata->Z.Pressure;
            fi.Z.u3         = idata->Z.Density * idata->Z.u2 * idata->Z.u3;
            fi.Z.Pressure   = idata->Z.Density * idata->Z.u2 * (idata->Z.E + idata->Z.Pressure / idata->Z.Density);

            fn.Z.Density    = ndata->Z.Density * ndata->Z.u2;
            fn.Z.u1         = ndata->Z.Density * ndata->Z.u2 * ndata->Z.u1;
            fn.Z.u2         = ndata->Z.Density * ndata->Z.u2 * ndata->Z.u2 + ndata->Z.Pressure;
            fn.Z.u3         = ndata->Z.Density * ndata->Z.u2 * ndata->Z.u3;
            fn.Z.Pressure   = ndata->Z.Density * ndata->Z.u2 * (ndata->Z.E + ndata->Z.Pressure / ndata->Z.Density);
            break;
        case 4:                      /* -z side */
        case 5:                      /* +z side */
            SC_LDEBUGF("face %d for -z and +z side\n", nface);
            fi.Z.Density    = idata->Z.Density * idata->Z.u3;
            fi.Z.u1         = idata->Z.Density * idata->Z.u3 * idata->Z.u1;
            fi.Z.u2         = idata->Z.Density * idata->Z.u3 * idata->Z.u2;
            fi.Z.u3         = idata->Z.Density * idata->Z.u3 * idata->Z.u3 + idata->Z.Pressure;
            fi.Z.Pressure   = idata->Z.Density * idata->Z.u3 * (idata->Z.E + idata->Z.Pressure / idata->Z.Density);

            fn.Z.Density    = ndata->Z.Density * ndata->Z.u3;
            fn.Z.u1         = ndata->Z.Density * ndata->Z.u3 * ndata->Z.u1;
            fn.Z.u2         = ndata->Z.Density * ndata->Z.u3 * ndata->Z.u2;
            fn.Z.u3         = ndata->Z.Density * ndata->Z.u3 * ndata->Z.u3 + ndata->Z.Pressure;
            fn.Z.Pressure   = ndata->Z.Density * ndata->Z.u3 * (ndata->Z.E + ndata->Z.Pressure / ndata->Z.Density);
            break;
        default:
            SC_ABORT("Wrong face");
    }

    SC_INFOF("[fi] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                   fi.Z.Density, fi.Z.u1, fi.Z.u2, fi.Z.u3, fi.Z.Pressure);

    SC_INFOF("[fn] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                   fn.Z.Density, fn.Z.u1, fn.Z.u2, fn.Z.u3, fn.Z.Pressure);

    fsum = sumZ(p8est, &fi, &fn);

    SC_INFOF("[fsum] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                   fsum.Z.Density, fsum.Z.u1, fsum.Z.u2, fsum.Z.u3, fsum.Z.Pressure);

    F.Z.Density     = 0.5 * (fsum.Z.Density - (ndata->Q.D - ndata->Q.D));
    F.Z.u1          = 0.5 * (fsum.Z.u1 - (ndata->Q.D - ndata->Q.D));
    F.Z.u2          = 0.5 * (fsum.Z.u2 - (ndata->Q.D - ndata->Q.D));
    F.Z.u3          = 0.5 * (fsum.Z.u3 - (ndata->Q.D - ndata->Q.D));
    F.Z.Pressure    = 0.5 * (fsum.Z.Pressure - (ndata->Q.D - ndata->Q.D));

    SC_INFOF("[FLUX] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                   F.Z.Density, F.Z.u1, F.Z.u2, F.Z.u3, F.Z.Pressure);

    return F;
}

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
    element_data_t              from_neighbor_flux;
    element_data_t              sum_flux;
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

    SC_LDEBUG("*****************\n");

    // прохождение по всем локальным ячейкам
    for (qumid = 0; qumid < mesh->local_num_quadrants; ++qumid) {
        which_tree = -1;

        // обнуляем сумму потоков для каждой ячейки
        init_solver_by_double(p8est, &sum_flux, 0.);

        // получение текущей ячейки
        q = p8est_mesh_quadrant_cumulative(p8est,
                                           qumid,
                                           &which_tree,
                                           &quadrant_id);

        h = (double) P4EST_QUADRANT_LEN (q->level) / (double) P4EST_ROOT_LEN;
        V = h*h*h;

        // даные текущей ячейки
        data = (element_data_t *) q->p.user_data;

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

            SC_INFOF("[neighbor GET] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf; dummy: %lf\n",
                           ndata->Z.Density,
                           ndata->Z.u1, ndata->Z.u2, ndata->Z.u3,
                           ndata->Z.Pressure, ndata->dummy);

            P4EST_ASSERT (p8est_quadrant_is_equal (qn, &(ndata->quad)));
            P4EST_ASSERT (ndata->quad.p.which_tree == which_tree);
            P4EST_ASSERT (!isnan(ndata->Z.Pressure));

            // Граничные значения забираем в соответствии с функцией контекста
            if (is_boundary) {
                SC_INFOF("[BOUNDARY] Data: %lf, nface: %d, nrank: %d, cur rank: %d, neib_face: %d\n",
                           ndata->dummy, nface, nrank, p8est->mpirank, neib_face);

                // подсчёт потока с границы
                element_data_t boundary_data = get_boundary_data_by_face(p8est, qn, nface);
                SC_INFOF("New data from boundary: %lf\n", boundary_data.dummy);
                SC_INFOF("[boundary_data] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                               boundary_data.Z.Density,
                               boundary_data.Z.u1, boundary_data.Z.u2, boundary_data.Z.u3,
                               boundary_data.Z.Pressure);

                // подсчёт потока с соседом
                from_neighbor_flux = calc_flux(p8est, data, &boundary_data, nface);

                // TODO test
                if (neib_face == 3 && data->dummy == 0) {
                    data->dummy = boundary_data.dummy;
                }
            } else {
                SC_INFOF("Data: %lf, nface: %d, nrank: %d, cur rank: %d, neib_face: %d\n",
                           ndata->dummy, nface, nrank, p8est->mpirank, neib_face);

                // подсчёт потока с соседом
                from_neighbor_flux = calc_flux(p8est, data, ndata, neib_face);

                // TODO test
                // перенос сверху вниз
                if (neib_face == 3) {
                    data->dummy += ndata->dummy;
                }
            }

            sum_flux.Z.Pressure  += from_neighbor_flux.Z.Pressure * Sn;
            sum_flux.Z.Density   += from_neighbor_flux.Z.Density * Sn;
            sum_flux.Z.u1        += from_neighbor_flux.Z.u1 * Sn;
            sum_flux.Z.u2        += from_neighbor_flux.Z.u2 * Sn;
            sum_flux.Z.u3        += from_neighbor_flux.Z.u3 * Sn;
        }

        sum_flux.Z.Pressure       = sum_flux.Z.Pressure * ctx->dt / V;
        sum_flux.Z.Density        = sum_flux.Z.Density * ctx->dt / V;
        sum_flux.Z.u1             = sum_flux.Z.u1 * ctx->dt / V;
        sum_flux.Z.u2             = sum_flux.Z.u2 * ctx->dt / V;
        sum_flux.Z.u3             = sum_flux.Z.u3 * ctx->dt / V;

        SC_INFOF("[old z] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                       data->Z.Density, data->Z.u1, data->Z.u2, data->Z.u3, data->Z.Pressure);

        SC_INFOF("[sum flux z after] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                       sum_flux.Z.Density, sum_flux.Z.u1, sum_flux.Z.u2, sum_flux.Z.u3, sum_flux.Z.Pressure);

        //
        data->Z.Pressure          = data->Z.Pressure - sum_flux.Z.Pressure;
        data->Z.Density           = data->Z.Density - sum_flux.Z.Density;
        data->Z.u1                = data->Z.u1 - sum_flux.Z.u1;
        data->Z.u2                = data->Z.u2 - sum_flux.Z.u2;
        data->Z.u3                = data->Z.u3 - sum_flux.Z.u3;

        SC_INFOF("[new z] Density: %lf; u1: %lf; u2: %lf; u3: %lf; Pressure: %lf\n",
                       data->Z.Density, data->Z.u1, data->Z.u2, data->Z.u3, data->Z.Pressure);

        SC_LDEBUG("*****************\n");
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

            SC_PRODUCTIONF("Quad is (%d), neib (%d), value (%d)\n", q->p.which_tree, nq[P8EST_HALF].p.which_tree, data->Z.Pressure);
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

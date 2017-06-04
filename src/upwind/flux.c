#include <p8est_mesh.h>
#include "solver.h"

void calc_flux_face_iter (p8est_iter_face_info_t *info, void *user_data)
{
    int                 i, j;
    p8est_t             *p8est = info->p4est;
    context_t           *ctx = (context_t *) p8est->user_pointer;
    element_data_t      *ghost_data = (element_data_t *) user_data;
    element_data_t      *udata;
    p8est_quadrant_t    *quad;
    double              vdotn = 0.;
    double              uavg;
    double              q;
    double              h, facearea;
    int                 which_face;
    int                 upwindside;
    p8est_iter_face_side_t *side[2];
    sc_array_t         *sides = &(info->sides);

    side[0] = p8est_iter_fside_array_index_int (sides, 0);
    if (sides->elem_count == 1) {
        /* наткнулись на границу, копируем текущую сторону */
        side[1] = p8est_iter_fside_array_index_int (sides, 0);
    } else {
        side[1] = p8est_iter_fside_array_index_int (sides, 1);
    }

    /* which of the quadrant's faces the interface touches */
    which_face = side[0]->face;

    switch (which_face) {
        case 0:                      /* -x side */
            vdotn = -ctx->v[0];
            break;
        case 1:                      /* +x side */
            vdotn = ctx->v[0];
            break;
        case 2:                      /* -y side */
            vdotn = -ctx->v[1];
            break;
        case 3:                      /* +y side */
            vdotn = ctx->v[1];
            break;
        case 4:                      /* -z side */
            vdotn = -ctx->v[2];
            break;
        case 5:                      /* +z side */
            vdotn = ctx->v[2];
            break;
    }
    upwindside = vdotn >= 0. ? 0 : 1;

    /* Because we have non-conforming boundaries, one side of an interface can
     * either have one large ("full") quadrant or 2^(d-1) small ("hanging")
     * quadrants: we have to compute the average differently in each case.  The
     * info populated by p4est_iterate() gives us the context we need to
     * proceed. */
    uavg = 0;
    if (side[upwindside]->is_hanging) {
        /* there are 2^(d-1) (P4EST_HALF) subfaces */
        for (j = 0; j < P8EST_HALF; j++) {
            if (side[upwindside]->is.hanging.is_ghost[j]) {
                udata = (element_data_t *) &ghost_data[side[upwindside]->is.hanging.quadid[j]];
            } else {
                udata = (element_data_t *) side[upwindside]->is.hanging.quad[j]->p.user_data;
            }
            uavg += udata->u;
        }
        uavg /= P8EST_HALF;
    }
    else {
        if (side[upwindside]->is.full.is_ghost) {
            udata = (element_data_t *) & ghost_data[side[upwindside]->is.full.quadid];
        } else {
            udata = (element_data_t *) side[upwindside]->is.full.quad->p.user_data;
        }
        uavg = udata->u;
    }

    /* flux from side 0 to side 1 */
    q = vdotn * uavg;
    for (i = 0; i < 2; i++) {
        if (side[i]->is_hanging) {
            /* there are 2^(d-1) (P4EST_HALF) subfaces */
            for (j = 0; j < P8EST_HALF; j++) {
                quad = side[i]->is.hanging.quad[j];
                h = (double) P8EST_QUADRANT_LEN (quad->level) / (double) P8EST_ROOT_LEN;

                facearea = h * h;

                if (!side[i]->is.hanging.is_ghost[j]) {
                    udata = (element_data_t *) quad->p.user_data;

                    if (i == upwindside) {
                        udata->dudt += vdotn * udata->u * facearea * (i ? 1. : -1.);
                    } else {
                        udata->dudt += q * facearea * (i ? 1. : -1.);
                    }
                }
            }
        } else {
            quad = side[i]->is.full.quad;
            h = (double) P8EST_QUADRANT_LEN (quad->level) / (double) P8EST_ROOT_LEN;

            facearea = h * h;

            if (!side[i]->is.full.is_ghost) {
                udata = (element_data_t *) quad->p.user_data;
                udata->dudt += q * facearea * (i ? 1. : -1.);
            }
        }
    }
}
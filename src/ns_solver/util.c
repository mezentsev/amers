#include <p8est_bits.h>
#include "util.h"

void get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3]) {
    p4est_qcoord_t half_length = P8EST_QUADRANT_LEN(q->level) / 2;

    p8est_qcoord_to_vertex (p8est->connectivity,
                            tree,
                            q->x + half_length,
                            q->y + half_length,
                            q->z + half_length,
                            xyz);
}

int is_quadrant_on_face_boundary (p8est_t * p4est,
                                  p4est_topidx_t treeid,
                                  int face,
                                  const p8est_quadrant_t * q) {
    p4est_qcoord_t      dh, xyz;
    p8est_connectivity_t *conn = p4est->connectivity;

    P4EST_ASSERT (0 <= face && face < P8EST_FACES);
    P4EST_ASSERT (p8est_quadrant_is_valid (q));

    if (conn->tree_to_tree[P8EST_FACES * treeid + face] != treeid ||
        (int) conn->tree_to_face[P8EST_FACES * treeid + face] != face) {
        return 0;
    }

    dh = P8EST_LAST_OFFSET (q->level);
    switch (face / 2) {
        case 0:
            xyz = q->x;
            break;
        case 1:
            xyz = q->y;
            break;
        case 2:
            xyz = q->z;
            break;
        default:
            SC_ABORT_NOT_REACHED ();
            break;
    }
    return xyz == ((face & 0x01) ? dh : 0);
}

int ipow(int base, int exp) {
    int result = 1;
    while (exp) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}
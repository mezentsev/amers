#include "util.h"

void get_midpoint(p8est_t *p8est, p4est_topidx_t tree, p8est_quadrant_t *q, double xyz[3])
{
    p4est_qcoord_t half_length = P8EST_QUADRANT_LEN(q->level) / 2;

    p8est_qcoord_to_vertex (p8est->connectivity,
                            tree,
                            q->x + half_length,
                            q->y + half_length,
                            q->z + half_length,
                            xyz);
}

double calc_speed_sound(double density, double pressure, double adiabatic) {
    return sqrt(adiabatic * pressure / density);
}
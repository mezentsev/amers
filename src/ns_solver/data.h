#ifndef AMR_DATA_H
#define AMR_DATA_H

#include <p4est_to_p8est.h>
#include <p8est.h>

/**
 * Data associated with a quadrant is accessible by p.user_data
 */
typedef struct data {
    double Density;
    double Pressure;
    double u1, u2, u3;
    double E;

    struct {
        double D;
        double Du1;
        double Du2;
        double Du3;
        double PE;
    } Q;

} data_t;


/**
 * Context data associated with a forest is accessible by user_pointer
 */
typedef struct context {
    double center[P4EST_DIM];  // coordinates of the center
    double width;
    double level;

    double Adiabatic;
} context_t;

#endif //AMR_DATA_H

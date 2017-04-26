#ifndef AMR_DATA_H
#define AMR_DATA_H

/**
 * Data associated with a quadrant is accessible by p.user_data
 */
typedef struct data {

} data_t;

/**
 * Context data associated with a forest is accessible by user_pointer
 */
typedef struct ctx {
    double center[P4EST_DIM];  // coordinates of the center
    double width;
    double level;
} ctx_t;


#endif //AMR_DATA_H

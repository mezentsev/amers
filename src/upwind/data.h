#ifndef AMR_DATA_H
#define AMR_DATA_H

#include <p4est_to_p8est.h>
#include <p8est.h>

/**
 * Данные каждой ячейки, которые можно получить через p.user_data
 */
typedef struct data {
    double u;
    double dudt;            /* производная по времени */
    double dux, duy, duz;      /* значение скорости по x, y, z */
    double dummy;

    p8est_quadrant_t quad;
} element_data_t;


/**
 * Данные контекста, которые можно получить через user_pointer
 */
typedef struct context {
    double  center[P4EST_DIM];    /* центр сетки */
    int     steps;                /* количество шагов */
    double  width;

    double  dt;                   /* временной шаг */

    /**
     * Прототип функции для получения данных от границ сетки
     *
     * @param p8est сетка
     * @param q текущая просматриваемая ячейка
     * @param face номер относительно соседа
     * @return
     */
    element_data_t (*get_boundary_data_by_face)(p8est_t *p8est, p8est_quadrant_t *q, int face);
} context_t;

#endif //AMR_DATA_H

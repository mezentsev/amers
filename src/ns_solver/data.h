#ifndef AMR_DATA_H
#define AMR_DATA_H

#include <p4est_to_p8est.h>
#include <p8est.h>

/**
 * Данные каждой ячейки, которые можно получить через p.user_data
 */
typedef struct data {
    double Density;         /* плотность */
    double Pressure;        /* давление */
    double u1, u2, u3;      /* значение скорости по x, y, z */
    double E;               /* полная энергия */

    struct {
        double D;
        double Du1;
        double Du2;
        double Du3;
        double PE;
    } Q;                    /* значение Q */

    int dummy;

} element_data_t;


/**
 * Данные контекста, которые можно получить через user_pointer
 */
typedef struct context {
    double  center[P4EST_DIM];    /* центр сетки */
    int     level;                /* уровень адаптации TODO */
    double  width;

    double  dt;                   /* временной шаг */
    double  Adiabatic;            /* показатель адиабаты */

    /**
     * Прототип функции для получения данных от границ сетки
     *
     * @param quad текущая просматриваемая ячейка
     * @param face номер относительно соседа
     * @return
     */
    element_data_t (*get_boundary_data_by_face)(p8est_quadrant_t *quad, int face);
} context_t;

#endif //AMR_DATA_H

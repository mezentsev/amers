#ifndef AMR_SOLVER_H
#define AMR_SOLVER_H

#include "data.h"

/**
 * Инициализация солвера
 *
 * @param q параметр
 * @param ctx контекст приложения
 */
void init_solver(data_t *q, context_t *ctx);

/**
 * Вычисление временного шага (CFL) для ячейки со стороной h
 *
 * @param data
 * @param ctx контекст приложения
 * @param h сторона ячейки
 */
void cflq(data_t *data, context_t *ctx, double h);

/**
 * Вычисление потока
 *
 * @param data
 * @param nx
 * @param ny
 * @param nz
 * @return
 */
double flow(data_t *data, double nx, double ny, double nz);

/**
 * Преобразование из Z в Q
 *
 * @param data
 */
void setq(data_t *data);


/**
 * Вычисление скорости света
 *
 * @param density значение плотности
 * @param pressure значение давления
 * @param adiabatic параметр адиабаты
 * @return
 */
double calc_speed_sound(double density, double pressure, double adiabatic);

#endif //AMR_SOLVER_H

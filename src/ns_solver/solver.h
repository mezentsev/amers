#ifndef NS_SOLVER_H
#define NS_SOLVER_H

#include "data.h"
#include "ns.h"

/**
 * Решить шаг N
 *
 * @param p8est
 * @param step шаг
 */
void solver_step(p8est_t *p8est,
                 int step);

/**
 * Возвращает данные с границы
 *
 *        e-------f
 *       /|      /|
 *      / |     / |
 *     a--|----b  |
 *     |  g----|--h
 *     | /     | /
 *     c-------d
 *
 *           3   5
 *           +  +
 *           | /
 *           |/
 *  1 - -----.---- + 0
 *          /|
 *         / |
 *        -  -
 *       4   2
 *
 * @param quad граничная ячейка
 * @param face номер относительно границы
 * @return заполненные данные
 */
element_data_t get_boundary_data_by_face(p8est_t *p8est, p8est_quadrant_t *quad, int face);

/**
 * Обход всех ячеек для вычисления временного шага
 *
 * @param p8est_iter_volume_info_t данные ячейки
 * @param user_data
 */
void calc_cfl_timestep_iter(p8est_iter_volume_info_t *info,
                            void *user_data);

/**
 * Обход по сетке для вычисления потока
 *
 * @param mesh сетка
 */
void calc_flux_mesh_iter(p8est_t *p8est,
                         p8est_mesh_t *mesh,
                         p8est_ghost_t *ghost,
                         void *ghost_data);
/**
 * Инициализация солвера
 *
 * @param q параметр
 */
void init_solver(p8est_t *p8est, element_data_t *q);

/**
 * Инициализация солвера пустыми векторами
 */
void init_empty_solver(p8est_t *p8est, element_data_t *data);

/**
 * Инициализация солвера числом
 */
void init_solver_by_double(p8est_t *p8est, element_data_t *data, double val);

/**
 * Вычисление временного шага (CFL) для ячейки со стороной h
 *
 * @param data
 * @param ctx контекст приложения
 * @param h сторона ячейки
 */
void cflq(element_data_t *data, context_t *ctx, double h);

/**
 * Вычисление потока между двумя соседними ячейками
 *
 * @param p8est
 * @param idata данные текущей ячейки
 * @param nadata данные соседа
 * @param nface направление фейса
 * @return значение вектора потока между двумя ячейками
 */
element_data_t calc_flux(p8est_t *p8est, element_data_t *idata, element_data_t *ndata, int nface);

/**
 * Преобразование из Z в Q
 *
 * @param p8est
 * @param data
 */
void updateQ(p8est_t *p8est, element_data_t *data);

/**
 * Суммирует два Z
 *
 * @param p8est
 * @param z1
 * @param z2
 * @return
 */
element_data_t sumZ(p8est_t *p8est, element_data_t *z1, element_data_t *z2);

/**
 * Вычисление скорости света
 *
 * @param density значение плотности
 * @param pressure значение давления
 * @param adiabatic параметр адиабаты
 * @return
 */
double calc_speed(double density, double pressure, double adiabatic);

#endif

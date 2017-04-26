#include <math.h>
#include <sc.h>
#include "solver.h"
#include "data.h"

void calc_q(data_t *data){
    data->Q.D   = data->Density;
    data->Q.Du1 = data->Density * data->u1;
    data->Q.Du2 = data->Density * data->u2;
    data->Q.Du3 = data->Density * data->u3;
    data->Q.PE  = data->Pressure * data->E;
}

void init_solver(data_t *data, context_t *ctx){
    double e;

    data->Density = 1;
    data->Pressure = 1;
    data->u1 = 1;
    data->u2 = 1;
    data->u3 = 1;

    e = data->Pressure / (ctx->Adiabatic  - 1) * data->Density;
    data->E = e + (pow(data->u1, 2) + pow(data->u2, 2) + pow(data->u3, 2))/2;

    calc_q(data);
}

double cfl(data_t *data, p8est_quadrant_t *q){
    return 0.;
}
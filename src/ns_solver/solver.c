#include <math.h>
#include <sc.h>
#include "solver.h"

void init_solver(data_t *data, context_t *ctx){
    double e;

    data->Density = 1;
    data->Pressure = 1;
    data->u1 = 1;
    data->u2 = 1;
    data->u3 = 1;

    e = data->Pressure / (ctx->Adiabatic  - 1) * data->Density;
    data->E = e + (pow(data->u1, 2) + pow(data->u2, 2) + pow(data->u3, 2))/2;
}
#include <math.h>
#include <sc.h>
#include <assert.h>
#include "solver.h"
#include "data.h"
#include "util.h"
#include "math.h"

void setq(data_t *data){
    data->Q.D   = data->Density;
    data->Q.Du1 = data->Density * data->u1;
    data->Q.Du2 = data->Density * data->u2;
    data->Q.Du3 = data->Density * data->u3;
    data->Q.PE  = data->Pressure * data->E;
}

void newq(data_t *data, context_t *ctx, double cfl, double length, double neighbors_sum){
    data_t new_data;
    double V = pow(length, 3);

    // TODO implementation
}

double flow(data_t *data, double nx, double ny, double nz) {
    // TODO implementation
    return 0;
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

    setq(data);
}

void cflq(data_t *data, context_t *ctx, double length) {
    assert(length > 0);
    assert(data->Pressure > 0);
    assert(data->Density > 0);
    assert(data->u1 != 0);
    assert(data->u2 != 0);
    assert(data->u3 != 0);

    double t1;
    double t2;
    double t3;
    double dt;
    double speed_sound = calc_speed_sound(data->Density, data->Pressure, ctx->Adiabatic);

    t1 = length/fabs(data->u1) + speed_sound;
    t2 = length/fabs(data->u2) + speed_sound;
    t3 = length/fabs(data->u3) + speed_sound;

    dt = 1/(1/t1 + 1/t2 + 1/t3);

    ctx->dt = (ctx->dt == 0) ? dt : MIN(ctx->dt, dt);
}


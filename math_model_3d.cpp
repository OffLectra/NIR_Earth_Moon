#include "math_model_3d.h"

//вспомогательные функции
double gx(double x, double r) {
    return -mu*x/pow(r, 3);
}

double gy(double y, double r) {
    return -mu*y/pow(r, 3);
}

double gz(double z, double r) {
    return -mu*z/pow(r, 3);
}


//правые части
double dVx_dt(double gx) {
    return gx;
}

double dVy_dt(double gy) {
    return gy;
}

double dVz_dt(double gz) {
    return gz;
}

double dx_dt(double Vx) {
    return Vx;
}

double dy_dt(double Vy) {
    return Vy;
}

double dz_dt(double Vz) {
    return Vz;
}

double dm_dt(double beta) {
    return -beta;
}

double hight(double r, double R_Earth) {
    return r-R_Earth;
}

#ifndef MATH_MODEL_3D_H
#define MATH_MODEL_3D_H

#include "help_function.h"

//double mu_Earth = 398600.4481;

//функции правых частей
double dVx_dt(double gx);
double dVy_dt(double gy);
double dVz_dt(double gz);
double dx_dt(double Vx);
double dy_dt(double Vy);
double dz_dt(double Vz);
double dm_dt(double beta);

//вспомогательные функции
double gx(double x, double r);
double gy(double y, double r);
double gz(double z, double r);

double hight(double r, double R_Earth);

#endif // MATH_MODEL_3D_H

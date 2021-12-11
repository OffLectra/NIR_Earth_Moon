#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "help_function.h"

//double mu_Earth = 398600.4481;

//функции правых частей
double f_d_Vx_dt(double gx, double P, double m_tek, double gamma, double Fax = 0.0);
double f_d_Vy_dt(double gy, double P, double m_tek, double gamma, double Fay = 0.0);
double f_d_x_dt(double Vx);
double f_d_y_dt(double Vy);
double f_d_m_dt(double beta);

//вспомогательные функции
double f_gx(double x, double r);
double f_gy(double y, double r);
double f_r(double x, double y);
double f_gamma(double gamma_0, double d_gamma_dt, double t);

double f_V(double Vx, double Vy);
double f_TETA(double Vx, double Vy);
double f_h(double a);
double f_h_izm(double V, double r);
double f_hight(double r, double R_Earth);

double f_Fax (double V, double Vx, double rho, double Sbalxbezm, double m_tec); //Sbalx = Cx * Sm/2
double f_Fay (double V, double Vy, double rho, double Sbalxbezm, double m_tec); //Sbalx = Cx * Sm/2


#endif // FUNCTIONS_H

#include "functions.h"
#include "cmath"


//вспомогательные функции
double f_gx(double x, double r)
{
    return -mu*x/pow(r, 3);
}

double f_gy(double y, double r)
{
    return -mu*y/pow(r, 3);
}

double f_r(double x, double y)
{
    return sqrt(pow(x,2) + pow(y,2));
}

double f_gamma(double gamma_0, double d_gamma_dt, double t)
{
    return gamma_0 + d_gamma_dt*t;
}



//правые части
double f_d_Vx_dt(double gx, double P, double m_tek, double gamma, double Fax)
{
    return gx + ((P/1000.0)/m_tek)*cos(gamma) + Fax/m_tek;
}

double f_d_Vy_dt(double gy, double P, double m_tek, double gamma, double Fay)
{
    return gy + ((P/1000.0)/m_tek)*sin(gamma) + Fay/m_tek;
}

double f_d_x_dt(double Vx)
{
    return Vx;
}

double f_d_y_dt(double Vy)
{
    return Vy;
}

double f_d_m_dt(double beta)
{
    return -beta;
}

double f_V(double Vx, double Vy)
{
    return sqrt(Vx*Vx + Vy*Vy);
}

double f_TETA(double Vx, double Vy)
{
    return atan2(Vy,Vx);
}

double f_h(double a)
{
    return -mu/a;
}

double f_h_izm(double V, double r)
{
    return V*V - ((2*mu)/r);
}

double f_Fax(double V, double Vx, double rho, double Sbalxbezm, double m_tec)
{
    return  -Sbalxbezm/m_tec*V*Vx*rho*1000.0;                                            ////Проверь Размерности!!!!!
}

double f_Fay(double V, double Vy, double rho, double Sbalxbezm, double m_tec)
{
    return  -Sbalxbezm/m_tec*V*Vy*rho*1000.0;
}

double f_hight(double r, double R_Earth)
{
    return r-R_Earth;
}

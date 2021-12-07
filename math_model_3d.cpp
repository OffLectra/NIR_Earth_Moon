#include "math_model_3d.h"

//вспомогательные функции
double f_gx(double x, double r)
{
    return -mu*x/pow(r, 3);
}

double f_gy(double y, double r)
{
    return -mu*y/pow(r, 3);
}

double f_gz(double z, double r)
{
    return -mu*z/pow(r, 3);
}

//правые части
double f_d_Vx_dt(double gx)
{
    return gx;
}

double f_d_Vy_dt(double gy)
{
    return gy;
}

double f_d_Vz_dt(double gz)
{
    return gz;
}

double f_d_x_dt(double Vx)
{
    return Vx;
}

double f_d_y_dt(double Vy)
{
    return Vy;
}

double f_d_z_dt(double Vz)
{
    return Vz;
}

double f_d_m_dt(double beta)
{
    return -beta;
}

double f_hight(double r, double R_Earth)
{
    return r-R_Earth;
}

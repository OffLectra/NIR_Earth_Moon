#include "optima.h"

/*
Vector grad_RK4_upr(Vector upr, Vector dupr, NU_RK4 nu, Settings_RK4 settings) {
    int n = upr.length();

    Vector grad = Vector(n).fill(0.0);

    double f_neizm = integr_RK4_upr(nu,settings,upr);
    Vector upr_t = upr;
    for (int i = 0; i < n; ++i) {
        upr_t = upr;
        upr_t[i] += dupr[i];
        grad[i] = (integr_RK4_upr(nu, settings, upr_t) - f_neizm)/dupr[i];
    }
    return grad;
}
*/

Vector gradDescent::grad(modeling_flight_2D flight, Vector U)
{
    // flight - переменная с функцией расчёта искомого параметра (t - время)
    // U = [u1,u2, ...]

    int n = U.length();// Количество управляющих параметров

    Vector grad = Vector(n).fill(0.0); // пустой вектор градиента

    double f_0 = flight.propUpr(U); // Значение искомого параметра в начальной точке (для которой и ищем градиент)

    for (int i = 0; i < n; ++i) {
        Vector Utemp = U;
        Utemp[i] = Utemp[i] + dU[i]; //[u1+du1, u2, ...] --i+1-> [u1, u2+du2, ...]
        grad[i]  = (flight.propUpr(Utemp) - f_0)/dU[i];
    }

    return grad;
}

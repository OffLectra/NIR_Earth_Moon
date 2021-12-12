#include "optima.h"

/*
Vector grad_RK4_upr(Vector upr, Vector dupr, NU_RK4 nu, Settings_RK4 settings) {
    int n = upr.length();

    Vector grad = Vector(n).fill(0.0);
    Vector f_dx = Vector(n).fill(0.0);

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




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

Vector gradDescent::gradD(modeling_flight_2D flight, Vector upr0)
{
    double t_cur = 0, t_temp = 10;
    Vector upr_cur = upr0;

    int counter = 0;

    print(QString("gamma_0") + "\t" + "d_gamma_dt" + "\n");
    print(QString::number(upr_cur[0]*180/M_PI) + "\t" + QString::number(upr_cur[1]*180/M_PI) + "\n");

    while (fabs(t_cur-t_temp)>eps) {
        t_temp = t_cur;
        t_cur = flight.propUpr(upr_cur);
        if (isToPrint) {
            print(QString("t_cur = %1 с").arg(t_cur));
        }
        Vector grad_cur = grad(flight,upr_cur);
        Vector addU = mp_VxV(grad_cur,lam);
        upr_cur = min_VxV(upr_cur,addU);

        qApp->processEvents();

        if (counter % 1000 == 0){
            lam = mp_cxV(2,lam);
        }
        counter++;

        if (isToPrint) {
            QString out = QString::number(counter) + " iteration\n";
            out += QString("grad[0]") + "\t" + "t_grad[1]" + "\n";
            out += QString("%1\t%2\n").arg(grad_cur[0]).arg(grad_cur[1]);
            out += QString("Dgamma_0") + "\t" + "Dd_gamma_dt" + "\n";
            out += QString("%1\t%2\n").arg(addU[0]).arg(addU[1]);
            out += QString("gamma_0") + "\t" + "d_gamma_dt" + "\n";
            out += QString("%1\t%2\n").arg(upr_cur[0]*toDeg).arg(upr_cur[1]*toDeg);
            print(out);
        }

    }
    print(QString("gamma_0") + "\t" + "d_gamma_dt" + "\n");
    print(QString::number(upr_cur[0]*180/M_PI) + "\t" + QString::number(upr_cur[1]*180/M_PI) + "\n");
    return upr_cur;
}

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
        qApp->processEvents();
    }

    return grad;
}

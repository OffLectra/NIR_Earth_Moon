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

Vector gradDescent::gradD(modeling_flight_2D flight, Vector upr0) {
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
        upr_cur = min_VxV(upr_cur, addU);

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

Vector DFP::calcDFP(modeling_flight_2D flight, Vector u0)
{
    double t_u_k = 0, t_u_k_1=10; // t - потому что функционал
    Vector u_k_1 = u0;
    Vector u_k = Vector(2).fill(0.0);
    Vector grad_u_k_1 = Vector(2).fill(2.0);
    Vector grad_u_k = Vector(2).fill(0.0);
    Vector dw = {0.0,0.0};
    Vector du = {0.0,0.0};
    Vector dut = {0.0,0.0};
    Matr A = {{5E-3,0.0},{0.0,5E-10}};//{{1E-3,0.0},{0.0,2E-9}};
    int counter = 1;
    gradDescent gr(dU,vivod);
    gr.isToPrint = isToPrint;
    while (fabs(grad_u_k[u_gamma]-grad_u_k_1[u_gamma])>eps ||
           fabs(grad_u_k[u_dgdt] -grad_u_k_1[u_dgdt])>eps ||
           fabs(t_u_k-t_u_k_1)>eps ||
           fabs(u_k[u_gamma]-u_k_1[u_gamma])>eps ||
           fabs(u_k[u_dgdt]-u_k_1[u_dgdt])>eps) {
        // начало итерации поиска минимума функции
//        t_u_k = t_u_k_1;
//        t_u_k_1 = flight.propUpr(u_k_1); //integr_RK4_upr(nu_cur, settings, u_k_1); //1.
//        if (isToPrint) {
//            print(QString("t_u_k_1 = %1 с").arg(t_u_k_1));
//        }
        grad_u_k_1 = gr.grad(flight, u_k_1);//grad_RK4_upr(u_k_1,dU,nu_cur,settings); // 1. Вычисляем Градиент в точке X(к-1), начальной точке итерации

        Vector p_k_1 = mp_MxV(A, grad_u_k_1);                    // 2. Умножение матрицы А на вектор градиента в точке X(к-1) -> Вектор изменения наших параметров
        //Vector lp_k_1 = {-lambda1*p_k_1[u_gamma],-lambda2*p_k_1[u_dgdt]}; // 3. Расчет вектора приращения для поиска новой точки (минус потому что ищим минимум, антиградиент)
        Vector addU = mp_VxV(p_k_1,lam);
        u_k = min_VxV(u_k_1,addU); // 4. Находим новую (конечную) точку (итерации) X(к)

        grad_u_k = gr.grad(flight, u_k);//grad_RK4_upr(x_k,deltax,nu_cur,settings); // 5. Вычисляем Градиент в точке X(к), конечной точке итерации

        Matr dA = A;
        if (counter % 1 == 0) {
            dw = min_VxV(grad_u_k_1,grad_u_k); // 7. Вычисляем
            du = min_VxV(u_k,u_k_1);            // 8. Вычисляем
            dut = add_VxV(du,mp_MxV(A,dw));     // 9. Вычисляем  Здесь А еще не изменилась!!!
            dA = mp_cxM(1.0/smp_VxV(dw,dut),mp_VxVT(dut)); // 10. Вычисляем
            A = min_MxM(A,dA); // 11. Вычисляем новую матрицу
        }

        if (isnan(A[0][0]) || isinf(A[0][0])) { break; }
        if (isToPrint) {
            QString out = QString::number(counter) + " iteration\n";
            out += QString("A") + "\n";
            out += QString::number(A[0][0]) + "\t" + QString::number(A[0][1]) + "\n";
            out += QString::number(A[1][0]) + "\t" + QString::number(A[1][1]) + "\n";
            out += QString("gamma") + "\t" + "dgdt" + "\n";
            out += QString::number(u_k_1[u_gamma]*180/M_PI) + "\t" + QString::number(u_k_1[u_dgdt ]*180/M_PI) + "\n";
            out += QString("fx") + "\t\n";
            out += QString::number(t_u_k_1) +"\n";
            out += QString("grad[u_gamma]") + "\t" + "grad[u_dgdt]" + "\n";
            out += QString::number(grad_u_k_1[u_gamma]) + "\t" + QString::number(grad_u_k_1[u_dgdt]) + "\n";
            out += QString("p[u_gamma]") + "\t" + "p[u_dgdt]" + "\n";
            out += QString::number(p_k_1[u_gamma],'f',10) + "\t" + QString::number(p_k_1[u_dgdt],'f',10) + "\n";
            out += QString("lp[u_gamma]") + "\t" + "lp[u_dgdt]" + "\n";
            out += QString::number(addU[u_gamma],'f',10) + "\t" + QString::number(addU[u_dgdt],'f',10) + "\n";
            out += QString("gamma") + "\t" + "dgdt" + "\n";
            out += QString::number(u_k[u_gamma]*180/M_PI,'f',10) + "\t" + QString::number(u_k[u_dgdt ]*180/M_PI,'f',10) + "\n";
            out += QString("grad[u_gamma]k") + "\t" + "grad[u_dgdt]k" + "\n";
            out += QString::number(grad_u_k[u_gamma]) + "\t" + QString::number(grad_u_k[u_dgdt]) + "\n";
            out += QString("dw[u_gamma]k") + "\t" + "dw[u_dgdt]k" + "\n";
            out += QString::number(dw[u_gamma],'f',10.0) + "\t" + QString::number(dw[u_dgdt],'f',10) + "\n";
            out += QString("dx[u_gamma]") + "\t" + "dx[u_dgdt]" + "\n";
            out += QString::number(du[u_gamma],'f',10) + "\t" + QString::number(du[u_dgdt],'f',10) + "\n";
            out += QString("dxt[u_gamma]") + "\t" + "dxt[u_dgdt]" + "\n";
            out += QString::number(dut[u_gamma],'f',10.0) + "\t" + QString::number(dut[u_dgdt],'f',10) + "\n";
            out += QString("dA") + "\n";
            out += QString::number(dA[0][0],'f',10) + "\t" + QString::number(dA[0][1],'f',10) + "\n";
            out += QString::number(dA[1][0],'f',10) + "\t" + QString::number(dA[1][1],'f',10) + "\n";
            print(out);
        }
        t_u_k = flight.propUpr(u_k); //integr_RK4_upr(nu_cur, settings, u_k_1); //1.
        if (fabs(t_u_k-t_u_k_1)>20) {
            return u_k_1;
        }
        t_u_k_1 = t_u_k;

        if (isToPrint) {
            print(QString("t_u_k_1 = %1 с").arg(t_u_k_1));
        }
        u_k_1 = u_k; // Назначаем старой точкой новуюю для следующей итерации
        qApp->processEvents();
        counter++;
        // конец итерации поиска минимума

    }

    return u_k;
}

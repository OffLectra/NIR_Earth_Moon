#include "RK4_integrator.h"
#include "functions.h"
#include "atm_gost_.h"

using namespace std;

//Мат.модель без атмосферы
double integr_RK4(NU_RK4 nu, Settings_RK4 settings) {
    bool print_to_file = settings.is_export_data_to_file;
    double matr_K[4][d_LAST] = {0.0};
    cASK vec_ASK = nu.v_ASK, vec_ASK_temp = vec_ASK;    //1) явл.пер для хран век. сост. в момент оконч расч шага(y_i, y_i+1, ...);
                                                                 //2) явл врем перем для хран век сост на подшагах инт ((y, y1, y2, y3)_i, (y, y1, y2, y3)_i+1, ...);
    cKE vec_KE = nu.v_KE;
    ofstream Vivod_File(settings.file_name.toStdString(), ios_base::trunc);
    Vivod_File << fixed;
    Vivod_File.precision(16);

    QString razd = ";";
    if (print_to_file) {
      Vivod_File << "t"     << razd.toStdString() << "x"    << razd.toStdString() << "y"    << razd.toStdString() << "r"     << razd.toStdString() <<
                    "Vx"    << razd.toStdString() << "Vy"   << razd.toStdString() << "V"    << razd.toStdString() << "TETA"  << razd.toStdString() <<
                    "a"     << razd.toStdString() << "e"    << razd.toStdString() << "i"    << razd.toStdString() << "RAAN"  << razd.toStdString() <<
                    "omega" << razd.toStdString() << "u"    << razd.toStdString() << "h"    << razd.toStdString() <<
                    "m"     << razd.toStdString() << "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha" << razd.toStdString() << "ro" << razd.toStdString() <<
                    "hight" << razd.toStdString() << endl;
    }
    int shag = 0;
    double h = 0.0;

    double h_izm = f_h_izm(f_V(vec_ASK.Vx, vec_ASK.Vy), f_r(vec_ASK.x, vec_ASK.y));//
    double V = f_V(vec_ASK.Vx, vec_ASK.Vy);
    double r = f_r(vec_ASK.x, vec_ASK.y);
    double alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy) - f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK.Time))*toDeg;

//    t = 0.0;           //с
//    m0 = 8250.0;       //кг
//    P_ud = 326.0;      //с  (удельный импульс)
//    P = 20000.0;       //Н  (тяга)
//    mu_earth = 398600.4481;
//    H0 = 200.0;        //км
//    a_f = 220000.0;    //км (большая полуось лунной эллиптич. орбиты)
//    R_earth = 6371.3;  //км (в 1 приближении, сфера)

    bool EoR = false;
    bool EoS = true;

    double h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);

    while (!EoR) {    //1 итер - 1 шаг инт

        if (print_to_file && EoS) {
          Vivod_File <<
                        vec_ASK.Time   << razd.toStdString() <<
                        vec_ASK.x      << razd.toStdString() <<
                        vec_ASK.y      << razd.toStdString() <<
                        r              << razd.toStdString() <<
                        vec_ASK.Vx     << razd.toStdString() <<
                        vec_ASK.Vy     << razd.toStdString() <<
                        V              << razd.toStdString() <<
                        f_TETA(vec_ASK.Vx, vec_ASK.Vy)*toDeg << razd.toStdString() <<
                        vec_KE.a       << razd.toStdString() <<
                        vec_KE.e       << razd.toStdString() <<
                        vec_KE.i*toDeg << razd.toStdString() <<
                        vec_KE.RAAN*toDeg << razd.toStdString() <<
                        vec_KE.om*toDeg << razd.toStdString() <<
                        vec_KE.u*toDeg  << razd.toStdString() <<
                        h_izm           << razd.toStdString() <<
                        vec_ASK.m       << razd.toStdString() <<
                        nu.beta         << razd.toStdString() <<
                        f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK.Time)*toDeg << razd.toStdString() <<
                        alpha           << razd.toStdString() <<
                        0               << razd.toStdString() <<
                        h_tec           << razd.toStdString() <<
                        endl;
        }
// 1)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);
        matr_K[0][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time) );
        matr_K[0][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time));
        matr_K[0][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[0][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[0][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[0][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[0][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x + settings.dt*matr_K[0][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y + settings.dt*matr_K[0][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m + settings.dt*matr_K[0][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt * 0.5;

// 2)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);
        matr_K[1][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time) );
        matr_K[1][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time));
        matr_K[1][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[1][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[1][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt * matr_K[1][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt * matr_K[1][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x + settings.dt * matr_K[1][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y + settings.dt * matr_K[1][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m + settings.dt * matr_K[1][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt * 0.5;

// 3)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);
        matr_K[2][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time) );
        matr_K[2][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time));
        matr_K[2][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[2][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[2][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt * matr_K[2][d_Vx_dt];
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt * matr_K[2][d_Vy_dt];
        vec_ASK_temp.x  = vec_ASK.x + settings.dt * matr_K[2][d_x_dt];
        vec_ASK_temp.y  = vec_ASK.y + settings.dt * matr_K[2][d_y_dt];
        vec_ASK_temp.m  = vec_ASK.m + settings.dt * matr_K[2][d_m_dt];
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt;

// 4)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);
        matr_K[3][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time) );
        matr_K[3][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time));
        matr_K[3][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[3][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[3][d_m_dt ] = f_d_m_dt(nu.beta);

//---------------------------
        vec_ASK_temp.Vx   = vec_ASK.Vx   + settings.dt * (1.0/6.0)*(matr_K[0][d_Vx_dt] + 2*matr_K[1][d_Vx_dt] + 2*matr_K[2][d_Vx_dt] + matr_K[3][d_Vx_dt]);
        vec_ASK_temp.Vy   = vec_ASK.Vy   + settings.dt * (1.0/6.0)*(matr_K[0][d_Vy_dt] + 2*matr_K[1][d_Vy_dt] + 2*matr_K[2][d_Vy_dt] + matr_K[3][d_Vy_dt]);
        vec_ASK_temp.x    = vec_ASK.x    + settings.dt * (1.0/6.0)*(matr_K[0][d_x_dt ] + 2*matr_K[1][d_x_dt ] + 2*matr_K[2][d_x_dt ] + matr_K[3][d_x_dt ]);
        vec_ASK_temp.y    = vec_ASK.y    + settings.dt * (1.0/6.0)*(matr_K[0][d_y_dt ] + 2*matr_K[1][d_y_dt ] + 2*matr_K[2][d_y_dt ] + matr_K[3][d_y_dt ]);
        vec_ASK_temp.m    = vec_ASK.m    + settings.dt * (1.0/6.0)*(matr_K[0][d_m_dt ] + 2*matr_K[1][d_m_dt ] + 2*matr_K[2][d_m_dt ] + matr_K[3][d_m_dt ]);
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt;

        h_izm = f_h_izm(f_V(vec_ASK_temp.Vx, vec_ASK_temp.Vy), f_r(vec_ASK_temp.x, vec_ASK_temp.y));

        if (h_izm>=settings.hf ) {
            if (fabs(h_izm-settings.hf)>settings.eps) {
                EoS = false;
                settings.dt /= 10;
            } else {
                EoR = true;
            }

        } else {
            EoS = true;
        }


        if (shag>=1E5) {
            EoR = true;
        }
        if (EoS) {
            vec_ASK = vec_ASK_temp;
            vec_KE = AGESK_to_KE(vec_ASK);
            V = f_V(vec_ASK.Vx, vec_ASK.Vy);
            r = f_r(vec_ASK.x, vec_ASK.y);
            alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK.Time))*toDeg;

            shag++;
        } else {
            vec_ASK_temp = vec_ASK;
        }
    }

    if (print_to_file) {
      Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "r" << razd.toStdString() <<
                    "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "V" << razd.toStdString() << "TETA" << razd.toStdString() <<
                    "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                    "omega" <<razd.toStdString() << "u" << razd.toStdString() << "h" << razd.toStdString() <<
                    "m" << razd.toStdString() << "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha"  << razd.toStdString() << "ro" << razd.toStdString() <<
                    "hight" << razd.toStdString() << endl;
    }
    Vivod_File.close();

    return vec_ASK_temp.Time;
}

double integr_RK4_upr(NU_RK4 nu, Settings_RK4 settings, Vector upr) {
    bool print_to_file = settings.is_export_data_to_file;
    double matr_K[4][d_LAST] = {0.0};
    cASK vec_ASK = nu.v_ASK, vec_ASK_temp = vec_ASK;  //1)- явл.пер для хран век. сост. в момент оконч расч шага(y_i, y_i+1, ...);
                                                                 //2) явл врем перем для хран век сост на подшагах инт ((y, y1, y2, y3)_i, (y, y1, y2, y3)_i+1, ...);
    cKE vec_KE = nu.v_KE;
    ofstream Vivod_File(settings.file_name.toStdString(), ios_base::trunc);
    Vivod_File << fixed;
    Vivod_File.precision(16);

    QString razd = ";";
    if (print_to_file) {
      Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "r" << razd.toStdString() <<
                    "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "V" << razd.toStdString() << "TETA" << razd.toStdString() <<
                    "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                    "omega" <<razd.toStdString() << "u" << razd.toStdString() << "h" << razd.toStdString() <<
                    "m" << razd.toStdString() << "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha" << razd.toStdString() << "ro" << razd.toStdString() <<
                    "hight" << razd.toStdString() << endl;
    }
    int shag = 0;
    double h = 0.0;

    double h_izm = f_h_izm(f_V(vec_ASK.Vx, vec_ASK.Vy), f_r(vec_ASK.x, vec_ASK.y));//
    double V = f_V(vec_ASK.Vx, vec_ASK.Vy);
    double r = f_r(vec_ASK.x, vec_ASK.y);
    double alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK.Time))*toDeg;

//    t = 0.0;           //с
//    m0 = 8250.0;       //кг
//    P_ud = 326.0;      //с  (удельный импульс)
//    P = 20000.0;       //Н  (тяга)
//    mu_earth = 398600.4481;
//    H0 = 200.0;        //км
//    a_f = 220000.0;    //км (большая полуось лунной эллиптич. орбиты)
//    R_earth = 6371.3;  //км (в 1 приближении, сфера)

    bool EoR = false;
    bool EoS = true;

    double h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);

    while (!EoR) {    //1 итер - 1 шаг инт

        if (print_to_file && EoS) {
          Vivod_File <<
                        vec_ASK.Time  << razd.toStdString() <<
                        vec_ASK.x    << razd.toStdString() <<
                        vec_ASK.y<< razd.toStdString() <<
                        r << razd.toStdString() <<
                        vec_ASK.Vx<< razd.toStdString() <<
                        vec_ASK.Vy<< razd.toStdString() <<
                        V << razd.toStdString() <<
                        f_TETA(vec_ASK.Vx, vec_ASK.Vy)*toDeg << razd.toStdString() <<
                        vec_KE.a << razd.toStdString() <<
                        vec_KE.e << razd.toStdString() <<
                        vec_KE.i*toDeg << razd.toStdString() <<
                        vec_KE.RAAN*toDeg << razd.toStdString() <<
                        vec_KE.om*toDeg << razd.toStdString() <<
                        vec_KE.u*toDeg << razd.toStdString() <<
                        h_izm << razd.toStdString() <<
                        vec_ASK.m<< razd.toStdString() <<
                        nu.beta << razd.toStdString() <<
                        f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK.Time)*toDeg << razd.toStdString() <<
                        alpha << razd.toStdString() << 0 << razd.toStdString() << h_tec << razd.toStdString() <<
                        endl;
        }
// 1)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);
        matr_K[0][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[0][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[0][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[0][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[0][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[0][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[0][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x + settings.dt*matr_K[0][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y + settings.dt*matr_K[0][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m + settings.dt*matr_K[0][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt*0.5;

// 2)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);
        matr_K[1][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[1][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[1][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[1][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[1][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[1][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[1][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x + settings.dt* matr_K[1][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y + settings.dt* matr_K[1][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m + settings.dt* matr_K[1][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt*0.5;

// 3)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);
        matr_K[2][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[2][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[2][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[2][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[2][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[2][d_Vx_dt];
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[2][d_Vy_dt];
        vec_ASK_temp.x  = vec_ASK.x + settings.dt* matr_K[2][d_x_dt];
        vec_ASK_temp.y  = vec_ASK.y + settings.dt* matr_K[2][d_y_dt];
        vec_ASK_temp.m  = vec_ASK.m + settings.dt* matr_K[2][d_m_dt];
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt;

// 4)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);
        matr_K[3][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[3][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[3][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[3][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[3][d_m_dt ] = f_d_m_dt(nu.beta);

//---------------------------
        vec_ASK_temp.Vx   = vec_ASK.Vx   + settings.dt*(1.0/6.0)*(matr_K[0][d_Vx_dt] + 2*matr_K[1][d_Vx_dt] + 2*matr_K[2][d_Vx_dt] + matr_K[3][d_Vx_dt]);
        vec_ASK_temp.Vy   = vec_ASK.Vy   + settings.dt*(1.0/6.0)*(matr_K[0][d_Vy_dt] + 2*matr_K[1][d_Vy_dt] + 2*matr_K[2][d_Vy_dt] + matr_K[3][d_Vy_dt]);
        vec_ASK_temp.x    = vec_ASK.x    + settings.dt*(1.0/6.0)*(matr_K[0][d_x_dt ] + 2*matr_K[1][d_x_dt ] + 2*matr_K[2][d_x_dt ] + matr_K[3][d_x_dt ]);
        vec_ASK_temp.y    = vec_ASK.y    + settings.dt*(1.0/6.0)*(matr_K[0][d_y_dt ] + 2*matr_K[1][d_y_dt ] + 2*matr_K[2][d_y_dt ] + matr_K[3][d_y_dt ]);
        vec_ASK_temp.m    = vec_ASK.m    + settings.dt*(1.0/6.0)*(matr_K[0][d_m_dt ] + 2*matr_K[1][d_m_dt ] + 2*matr_K[2][d_m_dt ] + matr_K[3][d_m_dt ]);
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt;


        h_izm = f_h_izm(f_V(vec_ASK_temp.Vx, vec_ASK_temp.Vy), f_r(vec_ASK_temp.x, vec_ASK_temp.y));

        if (h_izm>=settings.hf ) {
            if (fabs(h_izm-settings.hf)>settings.eps)
            {
                EoS = false;
                settings.dt /= 10;
            } else {
                EoR = true;
            }

        } else {
            EoS = true;
        }


        if (shag>=1E5) {
            EoR = true;
        }
        if (EoS) {
            vec_ASK = vec_ASK_temp;
            vec_KE = AGESK_to_KE(vec_ASK);
            V = f_V(vec_ASK.Vx, vec_ASK.Vy);
            r = f_r(vec_ASK.x, vec_ASK.y);
            alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK.Time))*toDeg;

            shag++;
        } else {
            vec_ASK_temp = vec_ASK;
        }
    }

    if (print_to_file) {
      Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "r" << razd.toStdString() <<
                    "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "V" << razd.toStdString() << "TETA" << razd.toStdString() <<
                    "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                    "omega" <<razd.toStdString() << "u" << razd.toStdString() << "h" << razd.toStdString() <<
                    "m" << razd.toStdString() << "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha"  << razd.toStdString() << "ro" << razd.toStdString() <<
                    "hight" << razd.toStdString() << endl;
    }
    Vivod_File.close();

    return vec_ASK_temp.Time;
}

//Мат.модель с учетом атмосферы
void integr_RK4_atm(NU_RK4 nu, Settings_RK4 settings) {
    bool print_to_file = settings.is_export_data_to_file;
    double matr_K[4][d_LAST] = {0.0};
    cASK vec_ASK = nu.v_ASK, vec_ASK_temp = vec_ASK;  //1)- явл.пер для хран век. сост. в момент оконч расч шага(y_i, y_i+1, ...);
                                                                 //2) явл врем перем для хран век сост на подшагах инт ((y, y1, y2, y3)_i, (y, y1, y2, y3)_i+1, ...);
    cKE vec_KE = nu.v_KE;
    ofstream Vivod_File(settings.file_name.toStdString(), ios_base::trunc);
    Vivod_File << fixed;
    Vivod_File.precision(16);

    QString razd = ";";
    if (print_to_file) {
      Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "r" << razd.toStdString() <<
                    "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "V" << razd.toStdString() << "TETA" << razd.toStdString() <<
                    "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                    "omega" <<razd.toStdString() << "u" << razd.toStdString() << "h" << razd.toStdString() <<
                    "m" << razd.toStdString() << "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha"<< razd.toStdString() << "ro" << razd.toStdString() <<
                    "hight" << razd.toStdString() << "Fax" << razd.toStdString() << "Fay" << razd.toStdString() << endl;
    }
    int shag = 0;
    double h = 0.0;
    double h_tec = f_hight(f_r(vec_ASK.x, vec_ASK.y),Earth_R);

    double h_izm = f_h_izm(f_V(vec_ASK.Vx, vec_ASK.Vy), f_r(vec_ASK.x, vec_ASK.y));//
    double V = f_V(vec_ASK.Vx, vec_ASK.Vy);
    double r = f_r(vec_ASK.x, vec_ASK.y);
    double alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK.Time))*toDeg;


    /*
    for (int i = 0; i < 100; i++) {
        double Kp = 3.0;
        double h = rand_f(120.0, 1500.0);
        double F = rand_f(65.0, 260.0);
        double X0[3] = { 0.0 };
        X0[0] = h;
        std::cout << "Parameter F = " << F << "\t";
        double rho = atmosGOST_R_25645_166_2004(h, F, Kp, F, 0.0, X0, 0.0, 0.0, 0.0, 0.0);
        std::cout << std::endl;
    }*/


    double X0[3] = { 0.0 };
    X0[0] = nu.v_ASK.x;
    X0[1] = nu.v_ASK.y;

    double F0 = 250.0;

    double ro = atmosGOST_R_25645_166_2004(h_izm, F0, 3.0, F0, 0.0, X0, 0.0, 0.0, 0.0, 0.0);//atmosGOST_R_25645_166_2004(h_izm, F0_arr[3], 3.0, 0.0, 0.0, X0[3], 0.0, 0.0, 0.0, 0.0);

//    t = 0.0;           //с
//    m0 = 8250.0;       //кг
//    P_ud = 326.0;      //с  (удельный импульс)
//    P = 20000.0;       //Н  (тяга)
//    mu_earth = 398600.4481;
//    H0 = 200.0;        //км
//    a_f = 220000.0;    //км (большая полуось лунной эллиптич. орбиты)
//    R_earth = 6371.3;  //км (в 1 приближении, сфера)

    bool EoR = false;
    bool EoS = true;


    h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);
    while (!EoR) {    //1 итер - 1 шаг инт

        if (print_to_file && EoS) {
          Vivod_File <<
                        vec_ASK.Time  << razd.toStdString() <<
                        vec_ASK.x    << razd.toStdString() <<
                        vec_ASK.y<< razd.toStdString() <<
                        r << razd.toStdString() <<
                        vec_ASK.Vx<< razd.toStdString() <<
                        vec_ASK.Vy<< razd.toStdString() <<
                        V << razd.toStdString() <<
                        f_TETA(vec_ASK.Vx, vec_ASK.Vy)*toDeg << razd.toStdString() <<
                        vec_KE.a << razd.toStdString() <<
                        vec_KE.e << razd.toStdString() <<
                        vec_KE.i*toDeg << razd.toStdString() <<
                        vec_KE.RAAN*toDeg << razd.toStdString() <<
                        vec_KE.om*toDeg << razd.toStdString() <<
                        vec_KE.u*toDeg << razd.toStdString() <<
                        h_izm << razd.toStdString() <<
                        vec_ASK.m<< razd.toStdString() <<
                        nu.beta << razd.toStdString() <<
                        f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK.Time)*toDeg << razd.toStdString() <<
                        alpha << razd.toStdString() << ro << razd.toStdString() << h_tec << razd.toStdString() <<
                        f_Fax(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vx,ro,nu.Sbalxbezm,vec_ASK_temp.m) << razd.toStdString() <<
                        f_Fay(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vy,ro,nu.Sbalxbezm,vec_ASK_temp.m) << razd.toStdString() <<
                        endl;
        }
// 1)------------------------

        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);
        X0[0] = vec_ASK_temp.x;
        X0[1] = vec_ASK_temp.y;

        ro = atmosGOST_R_25645_166_2004(h_tec, F0, 3.0, F0, 0.0, X0, 0.0, 0.0, 0.0, 0.0);//atmosGOST_R_25645_166_2004(h_izm, F0_arr[3], 3.0, 0.0, 0.0, X0[3], 0.0, 0.0, 0.0, 0.0);

        matr_K[0][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fax(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vx,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[0][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fay(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vy,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[0][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[0][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[0][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[0][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[0][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x + settings.dt*matr_K[0][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y + settings.dt*matr_K[0][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m + settings.dt*matr_K[0][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt*0.5;

// 2)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);

        X0[0] = vec_ASK_temp.x;
        X0[1] = vec_ASK_temp.y;

        ro = atmosGOST_R_25645_166_2004(h_tec, F0, 3.0, F0, 0.0, X0, 0.0, 0.0, 0.0, 0.0);//atmosGOST_R_25645_166_2004(h_izm, F0_arr[3], 3.0, 0.0, 0.0, X0[3], 0.0, 0.0, 0.0, 0.0);

        matr_K[1][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fax(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vx,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[1][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fay(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vy,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[1][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[1][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[1][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[1][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[1][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x + settings.dt* matr_K[1][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y + settings.dt* matr_K[1][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m + settings.dt* matr_K[1][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt*0.5;

// 3)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);

        X0[0] = vec_ASK_temp.x;
        X0[1] = vec_ASK_temp.y;

        ro = atmosGOST_R_25645_166_2004(h_tec, F0, 3.0, F0, 0.0, X0, 0.0, 0.0, 0.0, 0.0);//atmosGOST_R_25645_166_2004(h_izm, F0_arr[3], 3.0, 0.0, 0.0, X0[3], 0.0, 0.0, 0.0, 0.0);

        matr_K[2][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fax(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vx,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[2][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fay(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vy,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[2][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[2][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[2][d_m_dt ] = f_d_m_dt(nu.beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + settings.dt*matr_K[2][d_Vx_dt];
        vec_ASK_temp.Vy = vec_ASK.Vy + settings.dt*matr_K[2][d_Vy_dt];
        vec_ASK_temp.x  = vec_ASK.x + settings.dt* matr_K[2][d_x_dt];
        vec_ASK_temp.y  = vec_ASK.y + settings.dt* matr_K[2][d_y_dt];
        vec_ASK_temp.m  = vec_ASK.m + settings.dt* matr_K[2][d_m_dt];
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt;

// 4)------------------------
        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y), Earth_R);

        X0[0] = vec_ASK_temp.x;
        X0[1] = vec_ASK_temp.y;

        ro = atmosGOST_R_25645_166_2004(h_tec, F0, 3.0, F0, 0.0, X0, 0.0, 0.0, 0.0, 0.0);//atmosGOST_R_25645_166_2004(h_izm, F0_arr[3], 3.0, 0.0, 0.0, X0[3], 0.0, 0.0, 0.0, 0.0);

        matr_K[3][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fax(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vx,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[3][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), nu.P, vec_ASK_temp.m, f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK_temp.Time), f_Fay(f_V(vec_ASK_temp.Vx,vec_ASK_temp.Vy),vec_ASK_temp.Vy,ro,nu.Sbalxbezm,vec_ASK_temp.m));
        matr_K[3][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[3][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[3][d_m_dt ] = f_d_m_dt(nu.beta);

//---------------------------
        vec_ASK_temp.Vx   = vec_ASK.Vx   + settings.dt*(1.0/6.0)*(matr_K[0][d_Vx_dt] + 2*matr_K[1][d_Vx_dt] + 2*matr_K[2][d_Vx_dt] + matr_K[3][d_Vx_dt]);
        vec_ASK_temp.Vy   = vec_ASK.Vy   + settings.dt*(1.0/6.0)*(matr_K[0][d_Vy_dt] + 2*matr_K[1][d_Vy_dt] + 2*matr_K[2][d_Vy_dt] + matr_K[3][d_Vy_dt]);
        vec_ASK_temp.x    = vec_ASK.x    + settings.dt*(1.0/6.0)*(matr_K[0][d_x_dt ] + 2*matr_K[1][d_x_dt ] + 2*matr_K[2][d_x_dt ] + matr_K[3][d_x_dt ]);
        vec_ASK_temp.y    = vec_ASK.y    + settings.dt*(1.0/6.0)*(matr_K[0][d_y_dt ] + 2*matr_K[1][d_y_dt ] + 2*matr_K[2][d_y_dt ] + matr_K[3][d_y_dt ]);
        vec_ASK_temp.m    = vec_ASK.m    + settings.dt*(1.0/6.0)*(matr_K[0][d_m_dt ] + 2*matr_K[1][d_m_dt ] + 2*matr_K[2][d_m_dt ] + matr_K[3][d_m_dt ]);
        vec_ASK_temp.Time = vec_ASK.Time + settings.dt;



        h_izm = f_h_izm(f_V(vec_ASK_temp.Vx, vec_ASK_temp.Vy), f_r(vec_ASK_temp.x, vec_ASK_temp.y));



        if (h_izm>=settings.hf ) {
            if (fabs(h_izm-settings.hf)>settings.eps) {
                EoS = false;
                settings.dt /= 10;
            } else {
                EoR = true;
            }

        } else {
            EoS = true;
        }

        h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);

        if (h_tec<=120.0 ) {
            EoR = true;
        }

        if (shag>=1E5) {
            EoR = true;
        }

        if (EoS) {
            vec_ASK = vec_ASK_temp;
            vec_KE = AGESK_to_KE(vec_ASK);
            V = f_V(vec_ASK.Vx, vec_ASK.Vy);
            r = f_r(vec_ASK.x, vec_ASK.y);
            alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(nu.gamma_0, nu.d_gamma_dt, vec_ASK.Time))*toDeg;

            shag++;
        } else {
            vec_ASK_temp = vec_ASK;
        }
    }

    if (print_to_file) {
      Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "r" << razd.toStdString() <<
                    "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "V" << razd.toStdString() << "TETA" << razd.toStdString() <<
                    "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                    "omega" <<razd.toStdString() << "u" << razd.toStdString() << "h" << razd.toStdString() <<
                    "m" << razd.toStdString() << "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha" << razd.toStdString() << "ro" << razd.toStdString() <<
                    "hight" << razd.toStdString() << "Fax" << razd.toStdString() << "Fay" << razd.toStdString() << endl;
    }
    Vivod_File.close();
}

//        for (int i = 0; i < 4; i++ ) //итератор
//        {
//            switch (i) {
//            case 0:
//            case 1:
//            case 2:
//            {
////                matr_K[i][d_Vx_dt] = f_d_Vx_dt(vec_ASK.Time, f_gx(vec_ASK.x, f_r(vec_ASK.x, vec_ASK.y)), P, );
////                matr_K[i][d_Vy_dt] = f_d_Vy_dt();
////                matr_K[i][d_x_dt ] = f_d_x_dt();
////                matr_K[i][d_y_dt ] = f_d_y_dt();
////                matr_K[i][d_m_dt ] = f_d_m_dt();
//            }

//            }

//        }

//            i += 1; //присвоение
//            i = 1;  //присвоение
//            i = i + 1; //+ арифметический оператор для int
//            i == 1; //логический
//            i++;    //присвоение
//            i || f  //логический
// == || && <= >= > < ! - лог. опер.
// = += -= *= /= ++ -- - оп.присв.



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

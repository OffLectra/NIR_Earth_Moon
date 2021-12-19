#ifndef RK4_INTEGRATOR_H
#define RK4_INTEGRATOR_H

#include <QString>
#include <QDateTime>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "perevod_.h"
#include "functions.h"

struct NU_RK4{
    QDateTime DT_0;
    cASK v_ASK;
    cKE v_KE;
    double P, P_ud, Wist;
    double gamma_0;
    double d_gamma_dt;
    double beta;
    double Sbalxbezm;

    NU_RK4 upd_with_KE(cKE newKE){
        v_KE  = newKE;
        v_ASK = KE_to_AGESK(newKE);
        return *this;
    }

    Vector U(){ return {gamma_0, d_gamma_dt};}
} ; //название типа структуры

struct Settings_RK4{
    double dt;
    double af; // Большая полуось конечного эллипса
    double hf; // Интеграл энергии конечный (целевой)
    double eps = 1E-6;
    bool is_export_data_to_file = true; //состояние "выводить ли данные в файл"
    QString file_name = "result.csv";   //имя файла

    Settings_RK4 upd_af(double af_new){
        af = af_new;
        hf = f_h(af_new);
        return *this;
    }
};

double integr_RK4(NU_RK4 nu, Settings_RK4 settings);
void integr_RK4_atm(NU_RK4 nu, Settings_RK4 settings);
double integr_RK4_upr(NU_RK4 nu, Settings_RK4 settings, Vector upr);

Vector grad_RK4_upr(Vector upr,Vector dupr, NU_RK4 nu, Settings_RK4 settings);

enum d_dt_K {
    d_Vx_dt = 0,
    d_Vy_dt,
    d_Vz_dt,
    d_x_dt,
    d_y_dt,
    d_z_dt,
    d_m_dt,
    d_LAST
};

#endif // RK4_INTEGRATOR_H

#ifndef PEREVOD__H
#define PEREVOD__H

#include <QObject>
#include <QDateTime>
#include <math.h>
#include <string>
#include <iostream>
#include <string.h>
#include "help_function.h"

typedef struct {
    QString NameToC = "Kep"; //Type of coordinate name
    double a = 0,p = 0,r = 0,e = 0,i = 0,
    RAAN = 0,om = 0,f = 0,u = 0,i_d = 0,
    RAAN_d = 0,om_d = 0,f_d = 0,u_d = 0,
    Vn = 0, Vr = 0, V = 0, Time = 0, m = 0;         //default all angles in degrees
    QString name[17] = {"a","p","r","e","i","RAAN","om","u","f","i_d",
                        "RAAN_d","om_d","u_d","f_d","Vn","Vr","V"};
    enum KE_Params {
        ke_a = 0,
        ke_p,
        ke_r,
        ke_e,
        ke_i,
        ke_RAAN,
        ke_om,
        ke_f,
        ke_u,
        ke_i_d,
        ke_RAAN_d,
        ke_om_d,
        ke_f_d,
        ke_u_d,
        ke_Vn,
        ke_Vr,
        ke_V,
        ke_LAST
    };

    int kol_params  = 6;
  //  QDateTime DateTime;   //переменная времени
    int rad_or_deg = 1;
    void to_Deg() {
        this->u_d    = (this->u)   *toDeg;
        this->f_d    = (this->f)   *toDeg;
        this->i_d    = (this->i)   *toDeg;
        this->om_d   = (this->om)  *toDeg;
        this->RAAN_d = (this->RAAN)*toDeg;
    }
    void to_Rad() {
        this->u    = (this->u_d)   *toRad;
        this->f    = (this->f_d)   *toRad;
        this->i    = (this->i_d)   *toRad;
        this->om   = (this->om_d)  *toRad;
        this->RAAN = (this->RAAN_d)*toRad;
    }
} Kep_param_vec;

// декларация, далее объявляется тип
typedef struct {
    QString NameToC = "ASK"; //Type of coordinate name
    double x = 0, y = 0, z = 0, Vx = 0, Vy = 0, Vz = 0, Time = 0, m = 0;
    QString name[7] = {"x","y","z","Vx","Vy","Vz","m"};
    int kol_params  = 7;
    //QDateTime DateTime;

    double r() {
        return f_r(x, y, z);
    }

    double V() {
        return f_V(Vx, Vy, Vz);
    }

} ASK_param_vec; //название типа структуры

Kep_param_vec AGESK_to_KE(ASK_param_vec v_AGESK); //функция перевода из агэск в кэ
ASK_param_vec KE_to_AGESK(Kep_param_vec v_KE);  //функция перевода из ке в агэск

double  check_f_rad     (double f);
double  check_u_rad     (double u);
double  check_inc_rad   (double inc);
double  check_omega_rad (double omega);
double  check_RAAN_rad  (double RAAN);


#endif // PEREVOD__H

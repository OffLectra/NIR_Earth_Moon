#ifndef HELP_FUNCTION_H
#define HELP_FUNCTION_H

#include <math.h>
#include <QObject>
#include <cmath>

////солянка из функций и классов
const double e_   = M_E;
const double pi   = M_PI;
const double pi_2 = M_PI_2;// число пи пополам
const double pi_4 = M_PI_4;
const double bm   = pow(10,-16);

const double toRad = pi/180.0;
const double toDeg = 180.0/pi;

const double mu = 398600.4481; // гравитационный параметр для Земли
const double Earth_R = 6371.14;
const double g0 = 9.80665;

typedef QVector<double> Vector;
typedef QVector<Vector> Matr;

enum upr_param {
    u_gamma = 0,
    u_dgdt,
    u_LAST
};

double f_r(double x, double y, double z);
double f_V(double Vx, double Vy, double Vz);
bool angle_compare(double u1, double u2, double tol);

// inline - функция описывается в том же месте, где объявлена (в хедере)

Vector mp_MxV(Matr M, Vector V);
Vector mp_cxV(double num, Vector V);
Matr mp_cxM(double num, Matr M);
double smp_VxV(Vector V1, Vector V2);
Vector add_VxV(Vector V1, Vector V2);
Matr add_MxM(Matr M1, Matr M2);
Vector min_VxV(Vector V1, Vector V2);
Matr min_MxM(Matr M1, Matr M2);
Matr mp_VxVT(Vector V);
Matr trans(Matr M);

class vXYZ {
public:
    vXYZ();
    vXYZ(double x0, double y0, double z0);
    double x = 0.0,y = 0.0,z = 0.0;
    enum coordinates {
        c_x = 0,
        c_y,
        c_z,
        LAST_coord
    };
    double len ();
    vXYZ normalize ();
    int kol_params = 3;
};

//void foo() {
//    vXYZ kozulia = vXYZ(1,2,3); // 1
//    vXYZ kozulia1(1,2,3);       // 2

//    vXYZ k = kozulia;           // 3

//    vXYZ k1;                    // 4
//    k1.x = 8;
//    k1.y = 9;
//    k1.z = 2;
//}

#endif // HELP_FUNCTION_H

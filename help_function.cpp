#include "help_function.h"


double f_r(double x, double y, double z) {
    return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}


double f_V(double Vx, double Vy, double Vz) {
    return sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
}

bool angle_compare(double u1, double u2, double tol) {
    if ((u1 > pi && u2 > pi) || (u1 < pi && u2 < pi)) {
        return (fabs(u1-u2)<tol)?true:false;
    } else {
        if (u1 > 3.0/2.0*pi && u2 < 1.0/2.0*pi) {
            return (fabs(u1-2*pi-u2)<tol)?true:false;
        } else if (u2 > 3.0/2.0*pi && u1 < 1.0/2.0*pi) {
            return (fabs(u1+2.0*pi-u2)<tol)?true:false;
        } else {
            return (fabs(u1-u2)<tol)?true:false;
        }
    }
}



















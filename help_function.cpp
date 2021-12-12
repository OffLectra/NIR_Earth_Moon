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


vXYZ::vXYZ() {} // конструктор по умолчанию

vXYZ::vXYZ(double x0, double y0, double z0) { // еще один конструктор
    x=x0;   y=y0;   z=z0;
}

double vXYZ::len() {
    return sqrt(x*x+y*y+z*z);
}

vXYZ vXYZ::normalize() {
    vXYZ out = *this;
    double l = len();
    if (l != 0.0) {
        out.x /= l;
        out.y /= l;
        out.z /= l;
    }
    return out;
}

Vector mp_MxV(Matr M, Vector V) {
    int rows = V.length(), columns = M.at(0).length();

    Vector out = Vector(V.length());

    for (int ix = 0; ix < rows; ix++) {
        out[ix] = 0.0;
        for (int jx = 0; jx < columns; jx++)
            out[ix] += M[ix][jx] * V[jx];
    }
    return out;
}

Vector mp_cxV(double num, Vector V) {

    for (int ix = 0; ix < V.length(); ix++)
    {
        V[ix] *= num;
    }
    return V;
}

Matr mp_cxM(double num, Matr M) {
    int rows = M.length(), columns = M.at(0).length();
    for(int row = 0; row < rows; ++row) {
        for(int col = 0; col < columns; ++col) {
            M[row][col] *= num;
        }
    }
    return M;
}

Vector add_VxV(Vector V1, Vector V2) {
    for (int ix = 0; ix < V1.length(); ix++) {
        V1[ix] += V2[ix];
    }
    return V1;
}

double smp_VxV(Vector V1, Vector V2) {
    double out = 0.0;
    for (int ix = 0; ix < V1.length(); ix++) {
        out += V1[ix]*V2[ix];
    }
    return out;
}

Matr add_MxM(Matr M1, Matr M2) {
    int rows = M1.length(), columns = M1.at(0).length();
    for(int row = 0; row < rows; ++row) {
        for(int col = 0; col < columns; ++col) {
            M1[row][col] += M2[row][col];
        }
    }
    return M1;
}

Vector min_VxV(Vector V1, Vector V2) {
    for (int ix = 0; ix < V1.length(); ix++) {
        V1[ix] -= V2[ix];
    }
    return V1;
}

Matr min_MxM(Matr M1, Matr M2) {
    int rows = M1.length(), columns = M1.at(0).length();
    for(int row = 0; row < rows; ++row) {
        for(int col = 0; col < columns; ++col) {
            M1[row][col] -= M2[row][col];
        }
    }
    return M1;
}

Matr mp_VxVT(Vector V) {
    int rows = V.length(), columns = V.length();

    Matr out = Matr(rows).fill(Vector(columns));

    for (int ix = 0; ix < rows; ix++) {
        for (int jx = 0; jx < columns; jx++)
            out[ix][jx] = V[ix] * V[jx];
    }
    return out;
}

Matr trans(Matr M) {                             // транспонирование
    int rows = M.length(), columns = M.at(0).length();
    int new_rows = columns, new_columns = rows;
    Matr MT = Matr(new_rows).fill(Vector(new_columns).fill(0.0));
    for(int row = 0; row < rows; ++row) {
        for(int col = 0; col < columns; ++col) {
            MT[col][row] = M[row][col];
        }
    }
    return MT;
}

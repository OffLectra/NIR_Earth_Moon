#ifndef HELP_FUNCTION_H
#define HELP_FUNCTION_H
#include "math.h"
#include "QObject"

////солянка из функций и классов
const double e_         = M_E;
const double pi         = M_PI;
const double pi_2       = M_PI_2;// пополам
const double pi_4       = M_PI_4;
const double bm = pow(10,-16);

const double toRad = pi/180.0;
const double toDeg = 180.0/pi;

const double mu = 398600.4481; //грав параметр для Земли
const double Earth_R = 6371.14;
const double g0 = 9.80665;

typedef QVector<double> Vector;
typedef QVector<Vector> Matr;

enum upr_param{
    u_gamma = 0,
    u_dgdt,
    u_LAST
};

inline Vector mp_MxV(Matr M, Vector V){
    int rows = V.length(), columns = M.at(0).length();

    Vector out = Vector(V.length());

    for (int ix = 0; ix < rows; ix++)
    {
        out[ix] = 0.0;
        for (int jx = 0; jx < columns; jx++)
            out[ix] += M[ix][jx] * V[jx];
    }
    return out;
}

inline Vector mp_cxV(double num, Vector V){

    for (int ix = 0; ix < V.length(); ix++)
    {
        V[ix] *= num;
    }
    return V;
}

inline Matr mp_cxM(double num, Matr M){
    int rows = M.length(), columns = M.at(0).length();
    for(int row = 0; row < rows; ++row)
    {
        for(int col = 0; col < columns; ++col)
        {
            M[row][col] *= num;
        }
    }
    return M;
}

inline double smp_VxV(Vector V1, Vector V2){
    double out = 0.0;
    for (int ix = 0; ix < V1.length(); ix++)
    {
        out += V1[ix]*V2[ix];
    }
    return out;
}

inline Vector add_VxV(Vector V1, Vector V2){

    for (int ix = 0; ix < V1.length(); ix++)
    {
        V1[ix] += V2[ix];
    }
    return V1;
}

inline Matr add_MxM(Matr M1, Matr M2){
    int rows = M1.length(), columns = M1.at(0).length();
    for(int row = 0; row < rows; ++row)
    {
        for(int col = 0; col < columns; ++col)
        {
            M1[row][col] += M2[row][col];
        }
    }
    return M1;
}

inline Vector min_VxV(Vector V1, Vector V2){

    for (int ix = 0; ix < V1.length(); ix++)
    {
        V1[ix] -= V2[ix];
    }
    return V1;
}

inline Matr min_MxM(Matr M1, Matr M2){
    int rows = M1.length(), columns = M1.at(0).length();
    for(int row = 0; row < rows; ++row)
    {
        for(int col = 0; col < columns; ++col)
        {
            M1[row][col] -= M2[row][col];
        }
    }
    return M1;
}

inline Matr mp_VxVT(Vector V){
    int rows = V.length(), columns = V.length();

    Matr out = Matr(rows).fill(Vector(columns));

    for (int ix = 0; ix < rows; ix++)
    {
        for (int jx = 0; jx < columns; jx++)
            out[ix][jx] = V[ix] * V[jx];
    }
    return out;
}

inline Matr trans(Matr M) //транспонирование
{
    int rows = M.length(), columns = M.at(0).length();
    int new_rows = columns, new_columns = rows;
    Matr MT = Matr(new_rows).fill(Vector(new_columns).fill(0.0));
    for(int row = 0; row < rows; ++row)
    {
        for(int col = 0; col < columns; ++col)
        {
            MT[col][row] = M[row][col];
        }
    }
    return MT;
}

class vXYZ
{
public:
    vXYZ() {}
    vXYZ(double x0,double y0,double z0, bool normalize = false) {
        x=x0;   y=y0;   z=z0;
    }

    double x = 0.0,y = 0.0,z = 0.0;
    bool norm = false;
    double l = 0.0;
    enum coordinates{
        c_x = 0,
        c_y,
        c_z,
        LAST_coord
    };

    double len (){
        if (norm){
            return l;
        }else {
            return sqrt(x*x+y*y+z*z);
        }
    }

    void normalize (){
        l = len();
        if (l != 0.0){
            x /= len();
            y /= len();
            z /= len();
            norm = true;
        }
    }

    int kol_params = 3;
};


#endif // HELP_FUNCTION_H

#ifndef MODELING_FLIGHT_RK4_H
#define MODELING_FLIGHT_RK4_H

#include <iomanip>
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include "math_model_3d.h"
//#include "help_function.h"
#include "RK4_integrator.h"

using namespace std;




class StepData {
public:
    cASK vASK;
    cKE vKE;
    QDateTime DT;
    double globalTime;
    double alpha, gamma, gt, TETA, beta; //2D

    StepData() {}
    StepData(cASK ask, cKE ke) : vASK(ask), vKE(ke) {}

    QString Str_All_Header(QString del);

    QString toStr_All(QString del);
};







class modeling_flight_3D {
    NU_RK4 data;
    Settings_RK4 config;
    QVector<StepData> calcData;

public:
    modeling_flight_3D(NU_RK4 nu, Settings_RK4 settings) {
        data = nu;
        config = settings;
    };

    QVector<StepData> getCalcData() {
        return calcData;
    }

    void printCalcDataToFile(QString filename);
    void propagate();
};





class modeling_flight_2D {
    NU_RK4 data;
    Settings_RK4 config;
    QVector<StepData> resultData;
    bool isSaveOldData = false;
    Vector upr;

public:
    modeling_flight_2D(NU_RK4 nu, Settings_RK4 settings) {
        data = nu;
        config = settings;
    };

    QVector<StepData> getCalcData() {
        return resultData;
    }

    Vector get_vU0(){return data.U();}

    void upd(){

    }


    void printCalcDataToFile(QString filename);
    void propagate();

    double propUpr(Vector U);
};





//    void propagate_with_atm();
//void propagate_with_upr(Vector upr);
#endif // MODELING_FLIGHT_RK4_H

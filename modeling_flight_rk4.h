#ifndef MODELING_FLIGHT_RK4_H
#define MODELING_FLIGHT_RK4_H
#include "RK4_integrator.h"
#include "iostream"
#include "math.h"
#include <iomanip>
#include <cmath>
#include "fstream"
#include "math_model_3d.h"
#include "help_function.h"

using namespace std;

class StepData{
public:
    ASK_param_vec vASK;
    Kep_param_vec vKE;
    QDateTime DT;
    double time;
    double alpha, gamma, gt, TETA, beta;

    StepData(){}
    StepData(ASK_param_vec ask, Kep_param_vec ke) : vASK(ask),vKE(ke){}
    QString toStr_All(QString del){
        return  QString::number(vASK.Time)      + del +
                QString::number(vASK.x)         + del +
                QString::number(vASK.y)         + del +
                QString::number(vASK.z)         + del +
                QString::number(vASK.r())       + del +
                QString::number(vASK.Vx)        + del +
                QString::number(vASK.Vy)        + del +
                QString::number(vASK.Vz)        + del +
                QString::number(vASK.V())       + del +
                QString::number(TETA)           + del +
                QString::number(alpha)          + del +
                QString::number(gamma)          + del +
                QString::number(gt)             + del +
                QString::number(vKE.a)          + del +
                QString::number(vKE.e)          + del +
                QString::number(vKE.i*toDeg)    + del +
                QString::number(vKE.RAAN*toDeg) + del +
                QString::number(vKE.om*toDeg)   + del +
                QString::number(vKE.u*toDeg)    + del +
                QString::number(vASK.m)         + del +
                QString::number(beta);

    }
};

class modeling_flight_RK4
{
    NU_RK4 data;
    Settings_RK4 config;
    QVector<StepData> calcData;

public:
    modeling_flight_RK4(NU_RK4 nu, Settings_RK4 settings) {
        data = nu;
        config = settings;
    };

    QVector<StepData> getCalcData(){
        return calcData;
    }

    void printCalcDataToFile(QString filename) {
        bool print_to_file = config.is_export_data_to_file;
        double matr_K[4][d_LAST] = {0.0};
        ofstream Vivod_File(config.file_name.toStdString(), ios_base::trunc);
        Vivod_File << fixed;
        Vivod_File.precision(16);

        QString razd = ";";
        if () {
          Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "z" << razd.toStdString() << "r" << razd.toStdString() <<
                        "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "Vz" << razd.toStdString() << "V" << razd.toStdString() <<
                        "TETA" << razd.toStdString() << "alpha" << razd.toStdString() << "gamma" << razd.toStdString() << "gt" << razd.toStdString() <<
                        "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                        "omega" <<razd.toStdString() << "u" << razd.toStdString() << "m" << razd.toStdString() <<

                        "h" << razd.toStdString() <<
                         "beta" << razd.toStdString() << "gamma"<< razd.toStdString() << "alpha" << razd.toStdString() << "ro" << razd.toStdString() <<
                        "hight" << razd.toStdString() << endl;
        }

        for (StepData step : calcData) {
            Vivod_File << step.toStr_All(razd).toStdString() << endl;
        }

        if (print_to_file) {
          Vivod_File << "t"    << razd.toStdString() << "x"  << razd.toStdString() << "y" << razd.toStdString() << "z" << razd.toStdString() << "r" << razd.toStdString() <<
                        "Vx"   << razd.toStdString() << "Vy" << razd.toStdString() << "Vz" << razd.toStdString() << "V" << razd.toStdString() <<
                        "TETA" << razd.toStdString() << "alpha" << razd.toStdString() << "gamma" << razd.toStdString() << "gt" << razd.toStdString() <<
                        "a" << razd.toStdString() << "e" << razd.toStdString() << "i" << razd.toStdString() << "RAAN" << razd.toStdString() <<
                        "omega" <<razd.toStdString() << "u" << razd.toStdString() << "m" << razd.toStdString() <<


                        endl;
        }
        Vivod_File.close();
    }

    void propagate_unperturbed();

    void propagate_with_atm();
    void propagate_with_upr(Vector upr);
};

#endif // MODELING_FLIGHT_RK4_H

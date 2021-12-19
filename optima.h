#ifndef OPTIMA_H
#define OPTIMA_H

#include <QObject>
#include <QApplication>
#include <QTextBrowser>
#include "functions.h"


#include <math.h>
#include "modeling_flight_rk4.h"

typedef function<double(Vector)> f_x; // Функция f(x), где x - вектор параметров x1,x2,...

class gradDescent {

public:
    Vector dU  = {1E-4,1E-7};
    double eps = 1e-6;
    Vector lam = {5E-4,2E-9};
    gradDescent (QTextBrowser *TB_vivod) {
        vivod = TB_vivod;
    }

    gradDescent (Vector _dU,
                 QTextBrowser *TB_vivod = nullptr) {
        vivod = TB_vivod;
        dU = _dU;
        isToPrint = (vivod!=nullptr)?true:false;
    }

    gradDescent (
            Vector _dU,
            double _eps,
            Vector _lam,
            QTextBrowser *TB_vivod = nullptr
            ){
        dU  = _dU;
        eps = _eps;
        lam = _lam;
        vivod = TB_vivod;
        isToPrint = (vivod!=nullptr)?true:false;
    }


    Vector grad(f_x func, Vector upr);
    Vector gradD(f_x func, Vector upr);

    Vector gradD(modeling_flight_2D flight, Vector upr);
    Vector grad(modeling_flight_2D flight,Vector upr);



    void print(QString msg){vivod->append(msg);}
    bool isToPrint = true;
    QTextBrowser *vivod = nullptr;
};

class DFP {

public:

    DFP () {
        isToPrint = false;
    }

    DFP (QTextBrowser *TB_vivod) {
        vivod = TB_vivod;
    }



    double eps = 1E-6; // точность проверки условия выхода их расчета (нашли минимум)
    Vector dU = {1E-4,1E-7};
    Vector lam = {1,1};
    Matr   mA = {{5E-3,0.0},{0.0,5E-10}};

    double dt_max = 20;

    bool checkVector(Vector x,Vector x_1, double n, double e){
        for (int i = 0; i < n; ++i) {
            if((fabs(x[i]-x_1[i])>e)){
                return true;
            }
        }
        return false;
    }




    Vector calcDFP(f_x func, Vector upr);

    Vector calcDFP(modeling_flight_2D flight, Vector upr);



    void print(QString msg){vivod->append(msg);}
    bool isToPrint = true;
    QTextBrowser *vivod = nullptr;
};

struct data_from_optima_calc
{
    data_from_optima_calc() {}
};

class OptimaFlight2D
{
public:
    NU_RK4 NU0;
    Settings_RK4 settings0;
    Vector da, a0, ak;

    OptimaFlight2D(
            NU_RK4 _NU0,
            Settings_RK4 _settings0,
            Vector _da,
            Vector _a0,
            Vector _ak,
            QTextBrowser *TB_vivod) :
        NU0(_NU0),
        settings0(_settings0),
        da(_da),
        a0(_a0),
        ak(_ak),
        vivod(TB_vivod)
    {}

    QStringLL solve(int n);


//    data_from_optima_calc iter(NU_RK4 NU,
//                               Settings_RK4 s,
//                               int n){
//        if (n != 0) {
//            DFP opt(vivod);
//            modeling_flight_2D raschet1(NU, s);
//            Vector U = opt.calcDFP(raschet,raschet.get_vU0());
//            double t = raschet.propUpr(U);
////            ui->TB_vivod->append(QString("Полет до а = %1:\t t3 = %2").arg(ak).arg(t3));
//        } else {
//            for (int a = aFrom; a < ak.at(n-1); ++a) {
//                DFP opt(vivod);
//                Vector U = opt.calcDFP(raschet,nu.U());
//                double t = raschet.propUpr(U);



//                modeling_flight_2D newRaschet;
//                newRaschet.

//                iter();
//            }
//        }
//    }


    double task2v(NU_RK4 NU,
                  Settings_RK4 settings,
                  Vector apk){
        double ap = apk[0]; // Запоминаем величину большой полуоси промежуточного эллипса
        double ak = apk[1]; // Запоминаем величину большой полуоси конечного эллипса

        Vector Uparam {NU.gamma_0, NU.d_gamma_dt};

        modeling_flight_2D raschet1(NU, settings.upd_af(ap));

        DFP opt(vivod);
        opt.isToPrint = false;
        Vector Uisk1 = opt.calcDFP(raschet1,Uparam);
        double t1 = raschet1.propUpr(Uisk1);

        modeling_flight_2D raschet2(
                    NU.upd_with_KE(raschet1.getCalcData().last().vKE.getKE_with_aem()),
                    settings.upd_af(ak));
        Vector Uisk2 = opt.calcDFP(raschet2,Uisk1);
        double t2 = raschet2.propUpr(Uisk2);

        if (isToPrint) {
            print(QString("Result for ap = %1:\n\t t = %2").
                                 arg(ap).arg(t1+t2));
        }

        return t1+t2;
    }


    void print(QString msg){vivod->append(msg);}
    bool isToPrint = true;
    QTextBrowser *vivod = nullptr;

};


#endif // OPTIMA_H

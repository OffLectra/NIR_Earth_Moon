#ifndef OPTIMA_H
#define OPTIMA_H

#include <QObject>
#include <QApplication>
#include <QTextBrowser>


#include <math.h>
#include "modeling_flight_rk4.h"

class gradDescent {

public:
    Vector dU  = {1E-4,1E-7};
    double eps = 1e-6;
    Vector lam = {5E-4,2E-9};
    gradDescent (QTextBrowser *TB_vivod) {
        vivod = TB_vivod;
    }



    Vector gradD(modeling_flight_2D flight, Vector upr);
    Vector grad(modeling_flight_2D flight,Vector upr);



    void print(QString msg){vivod->append(msg);}
    bool isToPrint = true;
    QTextBrowser *vivod = nullptr;
};


#endif // OPTIMA_H

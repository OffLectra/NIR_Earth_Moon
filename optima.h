#ifndef OPTIMA_H
#define OPTIMA_H

#include <QObject>
#include <QWidget>

#include <math.h>
#include "modeling_flight_rk4.h"

class gradDescent {

public:
    Vector dU = {1E-4,1E-7};

    gradDescent () {}


    Vector grad(modeling_flight_2D flight,Vector upr);



};


#endif // OPTIMA_H

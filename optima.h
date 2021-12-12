#ifndef OPTIMA_H
#define OPTIMA_H

#include <QObject>
#include <QWidget>

#include <math.h>
#include "math_model_3d.h"
#include "RK4_integrator.h"

class gradDescent {

public:
    gradDescent () {}

    Vector grad(Vector upr);



};


#endif // OPTIMA_H

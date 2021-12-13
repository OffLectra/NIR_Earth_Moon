#include "perevod_.h"


ASK_param_vec KE_to_AGESK(Kep_param_vec v_KE) {
    int ToO = 0;
    vXYZ r0,n0;

    ASK_param_vec v_AGESK;
//__перевод в радианы из градусов___________________________________________
//    v_KE.u   =   deg_to_rad(v_KE.u_d);
//    v_KE.i   =   deg_to_rad(v_KE.i_d);
//    v_KE.om  =   deg_to_rad(v_KE.om_d);
//    v_KE.RAAN =  deg_to_rad(v_KE.RAAN_d);     //добавил перевод из градусов в радианы
//_ToO_______________________________________________________________
//    v_AGESK.DateTime = v_KE.DateTime;

    if (v_KE.e == 0.0) { // надо изменить условие (наверное)
        ToO = 0;
    }
    if ((v_KE.e>0)and(v_KE.e<1)) {
        ToO=1;
    }
    if (v_KE.e==1.0) {
        ToO=2;
    }
    if (v_KE.e>1) {
        ToO=3;
    }
//___p____________________________________________________________________
    switch (ToO) {
        case 0: {
               v_KE.p=v_KE.a;
        }break;

        case 1: {
               v_KE.p=v_KE.a*(1-pow(v_KE.e,2));
        }break;

        case 2: {
               v_KE.p=0;
        }break;

        case 3: {
               v_KE.p=v_KE.a*(pow(v_KE.e,2)-1);
        }break;
    }



//__teta____________________________________________________________________
    v_KE.f = v_KE.u-v_KE.om;
    v_KE.f = check_f_rad(v_KE.f);

//__r____________________________________________________________________

    v_KE.r=v_KE.p/(1+v_KE.e*cos(v_KE.f));

//__r0,n0____________________________________________________________________

    r0.x=cos(v_KE.RAAN)*cos(v_KE.u)-sin(v_KE.RAAN)*sin(v_KE.u)*cos(v_KE.i);
    r0.y=sin(v_KE.RAAN)*cos(v_KE.u)+cos(v_KE.RAAN)*sin(v_KE.u)*cos(v_KE.i);
    r0.z=sin(v_KE.u)*sin(v_KE.i);
    n0.x=-cos(v_KE.RAAN)*sin(v_KE.u)-sin(v_KE.RAAN)*cos(v_KE.u)*cos(v_KE.i);
    n0.y=-sin(v_KE.RAAN)*sin(v_KE.u)+cos(v_KE.RAAN)*cos(v_KE.u)*cos(v_KE.i);
    n0.z=cos(v_KE.u)*sin(v_KE.i);

//__x,y,z_________________________________________________________________

    v_AGESK.x = v_KE.r*r0.x;
    v_AGESK.y = v_KE.r*r0.y;
    v_AGESK.z = v_KE.r*r0.z;

//--------------------------------------------------------------------------

//__vr,vn______________________________________________________________________
    v_KE.Vr=sqrt(mu/v_KE.p)*v_KE.e*sin(v_KE.f);      //подсвечиваются глобальные переменные
    v_KE.Vn=sqrt(mu/v_KE.p)*(1+v_KE.e*cos(v_KE.f));  //vn=c/r;

//__v__________________________________________________________________________
    v_AGESK.Vx=v_KE.Vr*r0.x+v_KE.Vn*n0.x;
    v_AGESK.Vy=v_KE.Vr*r0.y+v_KE.Vn*n0.y;
    v_AGESK.Vz=v_KE.Vr*r0.z+v_KE.Vn*n0.z;
    return v_AGESK;
}

Kep_param_vec AGESK_to_KE(ASK_param_vec v_AGESK) {
    Kep_param_vec v_KE;
    double h;
    double c1,c2,c3,c,f1,f2,f3,lap;
    int ToO = 0;

    //__Проверка______________________________________________________________
    v_KE.r = sqrt(v_AGESK.x*v_AGESK.x + v_AGESK.y*v_AGESK.y + v_AGESK.z*v_AGESK.z);
    v_KE.V = sqrt(v_AGESK.Vx*v_AGESK.Vx + v_AGESK.Vy*v_AGESK.Vy + v_AGESK.Vz*v_AGESK.Vz);
    //__h____________________________________________________________________
    h = pow(v_KE.V,2)-2*mu/v_KE.r;
    if (h<-bm) {
        ToO = 0;
    }
    if (fabs(h) < bm) {
        ToO = 1;
    }
    if (h>bm) {
        ToO = 2;
    }

    //__c____________________________________________________________________
    c1  =  v_AGESK.y * v_AGESK.Vz - v_AGESK.z * v_AGESK.Vy;
    c2  =  v_AGESK.z * v_AGESK.Vx - v_AGESK.x * v_AGESK.Vz;
    c3  =  v_AGESK.x * v_AGESK.Vy - v_AGESK.y * v_AGESK.Vx;
    c  =  sqrt(c1*c1 + c2*c2 + c3*c3);
    if (fabs(c) < bm) {
        c =  bm;
    }
    //__lap____________________________________________________________________
    f1  =  -mu * v_AGESK.x/v_KE.r - v_AGESK.Vz * c2 + v_AGESK.Vy * c3;
    f2  =  -mu * v_AGESK.y/v_KE.r - v_AGESK.Vx * c3 + v_AGESK.Vz * c1;
    f3  =  -mu * v_AGESK.z/v_KE.r - v_AGESK.Vy * c1 + v_AGESK.Vx * c2;
    lap  =  sqrt(f1*f1 + f2*f2 + f3*f3);

    //__a____________________________________________________________________
    v_KE.a = -mu/h;

    //__p_  ___________________________________________________________________
            //       v_KE.p = c*c/mu;
    switch (ToO) {
        case 0: {
               v_KE.p = v_KE.a*(1-pow(v_KE.e,2));
        }break;

        case 1: {
               v_KE.p = v_KE.r;
        }break;

        case 2: {
               v_KE.p = v_KE.a*(pow(v_KE.e,2)-1);
        }break;
    }

    //__e____________________________________________________________________
    v_KE.e = sqrt(1.0+pow(c/mu,2)*h);
    //__i____________________________________________________________________
    v_KE.i = acos(c3/c);
    v_KE.i = check_inc_rad(v_KE.i); //check_inc_rad(v_KE.i);
    if (!(fabs(v_KE.i) < bm)) {
        //__RAAN_________________________________________________________
        v_KE.RAAN =  atan2(c1/sqrt(c1*c1+c2*c2),
                     -c2/sqrt(c1*c1+c2*c2));
        v_KE.RAAN = check_RAAN_rad(v_KE.RAAN);
        if(fabs(v_KE.RAAN-2*pi)<bm){
            v_KE.RAAN = 0.0;
        }
        //__omega________________________________________________________
        if (fabs(v_KE.e) > 1E-6){
            v_KE.om  =  atan2(f3/(lap*sin(v_KE.i)),
                         (cos(v_KE.RAAN)*f1+sin(v_KE.RAAN)*f2)/(mu*v_KE.e));
        }else{
            v_KE.om = 0.0;
        }
        v_KE.om  =  check_omega_rad(v_KE.om);
        //__u____________________________________________________________
        v_KE.u  =  atan2(v_AGESK.z/(v_AGESK.r() * sin(v_KE.i)),
          (v_AGESK.x*cos(v_KE.RAAN)+v_AGESK.y*sin(v_KE.RAAN))/v_AGESK.r());
        v_KE.u  =  check_u_rad(v_KE.u);
    } else {
        v_KE.RAAN = 0.0;
        v_KE.om = atan2(f2,f1);
        v_KE.u = atan2(v_AGESK.y/v_AGESK.r(),v_AGESK.x/v_AGESK.r()) - v_KE.RAAN;
        v_KE.om  =  check_omega_rad(v_KE.om);
        v_KE.u  =  check_u_rad(v_KE.u);
    }

    //__teta____________________________________________________________________
    v_KE.f = v_KE.u-v_KE.om;
    v_KE.f = check_f_rad(v_KE.f);
    //__vr____________________________________________________________________
    v_KE.Vr = sqrt(mu/v_KE.p)*v_KE.e*sin(v_KE.f);
    //__vn____________________________________________________________________
          //vn = c/r;
    v_KE.Vn = sqrt(mu/v_KE.p)*(1+v_KE.e*cos(v_KE.f));
    //__перевод в градусы из радиан___________________________________________
    v_KE.u_d     =  (v_KE.u)   *toDeg;
    v_KE.f_d     =  (v_KE.f)   *toDeg;
    v_KE.i_d     =  (v_KE.i)   *toDeg;
    v_KE.om_d    =  (v_KE.om)  *toDeg;
    v_KE.RAAN_d  =  (v_KE.RAAN)*toDeg;     //добавил перевод из радиан в градусы

    return v_KE;
}

double  check_RAAN_rad   (double RAAN) {
    if (RAAN < 0.0){
        RAAN += 2.0*pi;
        if (fabs(2.0*pi-RAAN)<bm) {
            RAAN = 2.0*pi-bm;
        }
    }
    if (RAAN >= 2.0*pi) {
        RAAN -= 2.0*pi;
        if (fabs(RAAN)<bm) {
            RAAN = bm;
        }
    }
    if ((RAAN <= bm)|| (fabs(RAAN-2.0*pi)<bm)) {
            RAAN = bm;
    }
    return RAAN;
}

double  check_omega_rad  (double omega) {
    if (omega < 0.0) {
        omega += 2.0*pi;
    }
    if (omega > 2.0*pi) {
        omega -= 2.0*pi;
    }
    return omega;
}

double  check_inc_rad    (double inc) {
    if (inc < 0.0) {
        inc += pi;
    }
    if (inc > pi) {
        inc -= pi;
    }
    return inc;
}
double  check_u_rad      (double u) {
    if (u < 0.0) {
        u += 2.0*pi;
    }
    if (u > 2.0*pi) {
        u -= 2.0*pi;
    }
    return u;
}
double  check_f_rad      (double f) {
    if (f < 0.0) {
        f += 2.0*pi;
    }
    if (f>2*pi) {
        f -= 2.0*pi;
    }
    return f;
}

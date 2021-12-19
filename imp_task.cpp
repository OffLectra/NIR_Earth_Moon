#include "imp_task.h"

Vector ImpTask::dV_xap(cKE vKE0, cKE &vKEk, double beta, double Wist, double m0, double mu)
{
    double rp0 = vKE0.rp();
    double Vp0 = vKE0.Vp(mu);
    double Vpk = sqrt(mu*(2.0/rp0-1.0/vKEk.a));
    double dVid= Vpk-Vp0;
    vKEk.e = rp0*Vpk*Vpk/mu-1.0;
    double mf = m0*pow(e_,dVid*1000/(-Wist));
    double dm = mf-m0;
    double dt = fabs(dm/beta);
    double dVgrav= 0.0189*mu/pow(vKE0.a,3)*dVid*dt;

    return {dVid + dVgrav, mf,dt};
}

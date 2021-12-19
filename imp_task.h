#ifndef IMP_TASK_H
#define IMP_TASK_H
#include "perevod_.h"
#include <QApplication>
#include <QTextBrowser>



class ImpTask
{
    cKE vKE0, vKEk;
    double beta, Wist, m0, mu;
public:
    enum res_vals{
        dV,
        mf,
        dt
    };
    ImpTask(cKE _vKE0,
            cKE _vKEk,
            double _beta,
            double _Wist,
            double _m0,
            double _mu,
            QTextBrowser *TB_vivod = nullptr) :
        vKE0(_vKE0), vKEk(_vKEk),
        beta(_beta), Wist(_Wist),
        m0(_m0), mu(_mu){
        vivod = TB_vivod;
        isToPrint = (vivod!=nullptr)?true:false;

    }


    double findFulldVxap2V(double ap){
        cKE vKEp(ap,0.0,0.0,0.0,0.0,0.0);
        Vector R1 = dV_xap(vKE0,vKEp,beta,Wist,m0,mu);
        Vector R2 = dV_xap(vKEp,vKEk,beta,Wist,R1[mf],mu);
        return R1[dV]+R2[dV];
    }

    double findFulldVxap(Vector ap){
        double result = 0;
        cKE vKEp(ap.first(),0.0,0.0,0.0,0.0,0.0), vKEp1;
        Vector r = dV_xap(vKE0,vKEp,beta,Wist,m0,mu);
        result += r[dV];
        for (int i = 1; i < ap.length(); ++i) {
            vKEp1 = vKEp;
            vKEp = cKE(ap[i],0.0,0.0,0.0,0.0,0.0);
            r = dV_xap(vKEp1,vKEp,beta,Wist,r[mf],mu);
            result += r[dV];
        }
        r = dV_xap(vKEp,vKEk,beta,Wist,r[mf],mu);
        result += r[dV];
        return result;
    }

    Vector findFullVals(Vector ap){
        cKE vKEp(ap.first(),0.0,0.0,0.0,0.0,0.0), vKEp1;
        Vector r = dV_xap(vKE0,vKEp,beta,Wist,m0,mu);
        Vector result = r;
        double dVsum = r[dV];
        double dtsum = r[dt];
        for (int i = 1; i < ap.length(); ++i) {
            vKEp1 = vKEp;
            vKEp = cKE(ap[i],0.0,0.0,0.0,0.0,0.0);
            r = dV_xap(vKEp1,vKEp,beta,Wist,r[mf],mu);
            result.append(r);
            dVsum += r[dV];
            dtsum += r[dt];
        }
        r = dV_xap(vKEp,vKEk,beta,Wist,r[mf],mu);
        result.append(r);
        dVsum += r[dV];
        dtsum += r[dt];
        result.append({dVsum,dtsum});
        result = ap + result;
        return result;
    }


    Vector dV_xap(cKE vKE0, cKE &vKEk, double beta, double Wist, double m0, double mu);


    void print(QString msg){vivod->append(msg);}
    bool isToPrint = true;
    QTextBrowser *vivod = nullptr;

};



#endif // IMP_TASK_H

#include "modeling_flight_rk4.h"


void modeling_flight_3D::printCalcDataToFile(QString filename) {
    ofstream Vivod_File(filename.toStdString(), ios_base::trunc);  // путь.имя файла с его расширением, ios - настройка для работы с файлом
    Vivod_File << fixed; // настройка для вывода вещественных чисел с определенной точностью
    Vivod_File.precision(16);

    QString razd = ";";
    Vivod_File << calcData.first().Str_All_Header(razd).toStdString() << endl;

    for (StepData step : calcData) {
        Vivod_File << step.toStr_All(razd).toStdString() << endl;
    }
    // - ниже аналогичная конструкция
    //        for (int i = 0; i < calcData.length(); i++) {
    //            StepData step = calcData[i];
    //            Vivod_File << step.toStr_All(razd).toStdString() << endl;
    //        }

    Vivod_File << calcData[0].Str_All_Header(razd).toStdString() << endl;
    Vivod_File.close();
}

void modeling_flight_3D::propagate() {
    double m_K[4][d_LAST] = {0.0};
    cASK v_ASK = data.v_ASK, v_ASK_t = v_ASK;  //1)- явл.пер для хран век. сост. в момент оконч расч шага(y_i, y_i+1, ...);
    double dt = config.dt;                                                    //2) явл врем перем для хран век сост на подшагах инт ((y, y1, y2, y3)_i, (y, y1, y2, y3)_i+1, ...);
    cKE v_KE = data.v_KE;
    double beta = 0.0;

    calcData.append(StepData(v_ASK, v_KE));

    int shag = 0;

    bool EoR = false;
    bool EoS = true;


    while (!EoR) {    //1 итер - 1 шаг инт
        // 1)------------------------
        m_K[0][d_Vx_dt] = dVx_dt(gx(v_ASK.x, v_ASK.r()));
        m_K[0][d_Vy_dt] = dVy_dt(gy(v_ASK.y, v_ASK.r()));
        m_K[0][d_Vz_dt] = dVz_dt(gz(v_ASK.z, v_ASK.r()));
        m_K[0][d_x_dt ] = dx_dt(v_ASK.Vx);
        m_K[0][d_y_dt ] = dy_dt(v_ASK.Vy);
        m_K[0][d_z_dt ] = dz_dt(v_ASK.Vz);
        m_K[0][d_m_dt ] = dm_dt(beta);

        // 2)------------------------
        v_ASK_t.Vx = v_ASK.Vx + dt*m_K[0][d_Vx_dt]*0.5;
        v_ASK_t.Vy = v_ASK.Vy + dt*m_K[0][d_Vy_dt]*0.5;
        v_ASK_t.Vz = v_ASK.Vz + dt*m_K[0][d_Vz_dt]*0.5;
        v_ASK_t.x  = v_ASK.x  + dt*m_K[0][d_x_dt]*0.5;
        v_ASK_t.y  = v_ASK.y  + dt*m_K[0][d_y_dt]*0.5;
        v_ASK_t.z  = v_ASK.z  + dt*m_K[0][d_z_dt]*0.5;
        v_ASK_t.m  = v_ASK.m  + dt*m_K[0][d_m_dt]*0.5;
        v_ASK_t.Time = v_ASK.Time + dt * 0.5;

        m_K[1][d_Vx_dt] = dVx_dt(gx(v_ASK_t.x, v_ASK_t.r()));
        m_K[1][d_Vy_dt] = dVy_dt(gy(v_ASK_t.y, v_ASK_t.r()));
        m_K[1][d_Vz_dt] = dVz_dt(gz(v_ASK_t.z, v_ASK_t.r()));
        m_K[1][d_x_dt ] = dx_dt(v_ASK_t.Vx);
        m_K[1][d_y_dt ] = dy_dt(v_ASK_t.Vy);
        m_K[1][d_z_dt ] = dz_dt(v_ASK_t.Vz);
        m_K[1][d_m_dt ] = dm_dt(beta);

        // 3)------------------------
        v_ASK_t.Vx = v_ASK.Vx + dt * m_K[1][d_Vx_dt]*0.5;
        v_ASK_t.Vy = v_ASK.Vy + dt * m_K[1][d_Vy_dt]*0.5;
        v_ASK_t.Vz = v_ASK.Vz + dt * m_K[1][d_Vz_dt]*0.5;
        v_ASK_t.x  = v_ASK.x  + dt * m_K[1][d_x_dt]*0.5;
        v_ASK_t.y  = v_ASK.y  + dt * m_K[1][d_y_dt]*0.5;
        v_ASK_t.z  = v_ASK.z  + dt * m_K[1][d_z_dt]*0.5;
        v_ASK_t.m  = v_ASK.m  + dt * m_K[1][d_m_dt]*0.5;
        v_ASK_t.Time = v_ASK.Time + dt * 0.5;

        m_K[2][d_Vx_dt] = dVx_dt(gx(v_ASK_t.x, v_ASK_t.r()));
        m_K[2][d_Vy_dt] = dVy_dt(gy(v_ASK_t.y, v_ASK_t.r()));
        m_K[2][d_Vz_dt] = dVy_dt(gz(v_ASK_t.z, v_ASK_t.r()));
        m_K[2][d_x_dt ] = dx_dt(v_ASK_t.Vx);
        m_K[2][d_y_dt ] = dy_dt(v_ASK_t.Vy);
        m_K[2][d_z_dt ] = dz_dt(v_ASK_t.Vz);
        m_K[2][d_m_dt ] = dm_dt(beta);

        // 4)------------------------
        v_ASK_t.Vx = v_ASK.Vx + dt * m_K[2][d_Vx_dt];
        v_ASK_t.Vy = v_ASK.Vy + dt * m_K[2][d_Vy_dt];
        v_ASK_t.Vz = v_ASK.Vz + dt * m_K[2][d_Vz_dt];
        v_ASK_t.x  = v_ASK.x  + dt * m_K[2][d_x_dt];
        v_ASK_t.y  = v_ASK.y  + dt * m_K[2][d_y_dt];
        v_ASK_t.z  = v_ASK.z  + dt * m_K[2][d_z_dt];
        v_ASK_t.m  = v_ASK.m  + dt * m_K[2][d_m_dt];
        v_ASK_t.Time = v_ASK.Time + dt;

        m_K[3][d_Vx_dt] = dVx_dt(gx(v_ASK_t.x, v_ASK_t.r()));
        m_K[3][d_Vy_dt] = dVy_dt(gy(v_ASK_t.y, v_ASK_t.r()));
        m_K[3][d_Vz_dt] = dVz_dt(gz(v_ASK_t.z, v_ASK_t.r()));
        m_K[3][d_x_dt ] = dx_dt(v_ASK_t.Vx);
        m_K[3][d_y_dt ] = dy_dt(v_ASK_t.Vy);
        m_K[3][d_z_dt ] = dz_dt(v_ASK_t.Vz);
        m_K[3][d_m_dt ] = dm_dt(beta);

        //---------------------------
        v_ASK_t.Vx   = v_ASK.Vx   + dt * (1.0/6.0)*(m_K[0][d_Vx_dt] + 2*m_K[1][d_Vx_dt] + 2*m_K[2][d_Vx_dt] + m_K[3][d_Vx_dt]);
        v_ASK_t.Vy   = v_ASK.Vy   + dt * (1.0/6.0)*(m_K[0][d_Vy_dt] + 2*m_K[1][d_Vy_dt] + 2*m_K[2][d_Vy_dt] + m_K[3][d_Vy_dt]);
        v_ASK_t.Vz   = v_ASK.Vz   + dt * (1.0/6.0)*(m_K[0][d_Vz_dt] + 2*m_K[1][d_Vz_dt] + 2*m_K[2][d_Vz_dt] + m_K[3][d_Vz_dt]);
        v_ASK_t.x    = v_ASK.x    + dt * (1.0/6.0)*(m_K[0][d_x_dt ] + 2*m_K[1][d_x_dt ] + 2*m_K[2][d_x_dt ] + m_K[3][d_x_dt ]);
        v_ASK_t.y    = v_ASK.y    + dt * (1.0/6.0)*(m_K[0][d_y_dt ] + 2*m_K[1][d_y_dt ] + 2*m_K[2][d_y_dt ] + m_K[3][d_y_dt ]);
        v_ASK_t.z    = v_ASK.z    + dt * (1.0/6.0)*(m_K[0][d_z_dt ] + 2*m_K[1][d_z_dt ] + 2*m_K[2][d_z_dt ] + m_K[3][d_z_dt ]);
        v_ASK_t.m    = v_ASK.m    + dt * (1.0/6.0)*(m_K[0][d_m_dt ] + 2*m_K[1][d_m_dt ] + 2*m_K[2][d_m_dt ] + m_K[3][d_m_dt ]);
        v_ASK_t.Time = v_ASK.Time + dt;

// метод эйлера
//        v_ASK_t.Vx   = v_ASK.Vx   + config.dt * dVx_dt(gx(v_ASK_t.x, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
//        v_ASK_t.Vy   = v_ASK.Vy   + config.dt * dVy_dt(gy(v_ASK_t.y, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
//        v_ASK_t.Vz   = v_ASK.Vz   + config.dt * dVz_dt(gz(v_ASK_t.z, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
//        v_ASK_t.x    = v_ASK.x    + config.dt * dx_dt(v_ASK_t.Vx);
//        v_ASK_t.y    = v_ASK.y    + config.dt * dy_dt(v_ASK_t.Vy);
//        v_ASK_t.z    = v_ASK.z    + config.dt * dz_dt(v_ASK_t.Vz);
//        v_ASK_t.m    = v_ASK.m    + config.dt * dm_dt(data.beta);
//        v_ASK_t.Time = v_ASK.Time + config.dt;

        v_KE = AGESK_to_KE(v_ASK_t); //функция перевода из агэск в кэ

        //проверка выхода по положению на орбите
        bool isEqual = angle_compare(v_KE.om, v_KE.u, 3*toRad);
        if (isEqual && v_KE.u>v_KE.om) {
            if (!angle_compare(v_KE.om, v_KE.u, 0.01*toRad)) {
                EoS = false;
                dt /= 10.0;
            } else {
                EoR = true;
            }

        } else {
            EoS = true;
        }
        //Проверка, что программа не работает в холостую
        if (shag>=1E5) {
            EoR = true;
        }
        //Если данные удовлетворяют мы их сохраняем, иначе откат
        if (EoS) {
            v_ASK = v_ASK_t;
            v_KE = AGESK_to_KE(v_ASK);
            shag++;

            calcData.append(StepData(v_ASK,v_KE));

        } else {
            v_ASK_t = v_ASK;
        }
    }
}

void modeling_flight_2D::printCalcDataToFile(QString filename) {
    ofstream Vivod_File(filename.toStdString(), ios_base::trunc);  // путь.имя файла с его расширением, ios - настройка для работы с файлом
    Vivod_File << fixed;
    Vivod_File.precision(16);

    QString razd = ";";
    Vivod_File << resultData[0].Str_All_Header(razd).toStdString() << endl;

    for (StepData step : resultData) {
        Vivod_File << step.toStr_All(razd).toStdString() << endl;
    }
    // - ниже аналогичная конструкция
    //        for (int i = 0; i < calcData.length(); i++) {
    //            StepData step = calcData[i];
    //            Vivod_File << step.toStr_All(razd).toStdString() << endl;
    //        }

    Vivod_File << resultData[0].Str_All_Header(razd).toStdString() << endl;

    Vivod_File.close();
}

void modeling_flight_2D::propagate() {

    double matr_K[4][d_LAST] = {0.0};
    cASK vec_ASK = data.v_ASK, vec_ASK_temp = vec_ASK;  //1)- явл.пер для хран век. сост. в момент оконч расч шага(y_i, y_i+1, ...);
    //2) явл врем перем для хран век сост на подшагах инт ((y, y1, y2, y3)_i, (y, y1, y2, y3)_i+1, ...);
    cKE vec_KE = data.v_KE;
    double dt = config.dt;
    double beta = data.beta;

    if (!isSaveOldData) {
        resultData.clear();
    }
    resultData.append(StepData(vec_ASK, vec_KE));

    int shag = 0;

    double h_izm = f_h_izm(f_V(vec_ASK.Vx, vec_ASK.Vy), f_r(vec_ASK.x, vec_ASK.y));//
    double alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK.Time))*toDeg;

    bool EoR = false;
    bool EoS = true;

    double h_tec = f_hight(f_r(vec_ASK_temp.x, vec_ASK_temp.y),Earth_R);

    while (!EoR) {    //1 итер - 1 шаг инт

// 1)------------------------

        matr_K[0][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[0][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[0][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[0][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[0][d_m_dt ] = f_d_m_dt(beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + dt * matr_K[0][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + dt * matr_K[0][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x  + dt * matr_K[0][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y  + dt * matr_K[0][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m  + dt * matr_K[0][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + dt*0.5;

// 2)------------------------

        matr_K[1][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[1][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[1][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[1][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[1][d_m_dt ] = f_d_m_dt(beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + dt * matr_K[1][d_Vx_dt]*0.5;
        vec_ASK_temp.Vy = vec_ASK.Vy + dt * matr_K[1][d_Vy_dt]*0.5;
        vec_ASK_temp.x  = vec_ASK.x  + dt * matr_K[1][d_x_dt]*0.5;
        vec_ASK_temp.y  = vec_ASK.y  + dt * matr_K[1][d_y_dt]*0.5;
        vec_ASK_temp.m  = vec_ASK.m  + dt * matr_K[1][d_m_dt]*0.5;
        vec_ASK_temp.Time = vec_ASK.Time + dt*0.5;

// 3)------------------------

        matr_K[2][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[2][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[2][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[2][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[2][d_m_dt ] = f_d_m_dt(beta);

        vec_ASK_temp.Vx = vec_ASK.Vx + dt * matr_K[2][d_Vx_dt];
        vec_ASK_temp.Vy = vec_ASK.Vy + dt * matr_K[2][d_Vy_dt];
        vec_ASK_temp.x  = vec_ASK.x  + dt * matr_K[2][d_x_dt];
        vec_ASK_temp.y  = vec_ASK.y  + dt * matr_K[2][d_y_dt];
        vec_ASK_temp.m  = vec_ASK.m  + dt * matr_K[2][d_m_dt];
        vec_ASK_temp.Time = vec_ASK.Time + dt;

// 4)------------------------

        matr_K[3][d_Vx_dt] = f_d_Vx_dt(f_gx(vec_ASK_temp.x, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time) );
        matr_K[3][d_Vy_dt] = f_d_Vy_dt(f_gy(vec_ASK_temp.y, f_r(vec_ASK_temp.x, vec_ASK_temp.y)), data.P, vec_ASK_temp.m, f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK_temp.Time));
        matr_K[3][d_x_dt ] = f_d_x_dt(vec_ASK_temp.Vx);
        matr_K[3][d_y_dt ] = f_d_y_dt(vec_ASK_temp.Vy);
        matr_K[3][d_m_dt ] = f_d_m_dt(beta);

//---------------------------
        vec_ASK_temp.Vx   = vec_ASK.Vx   + dt*(1.0/6.0)*(matr_K[0][d_Vx_dt] + 2*matr_K[1][d_Vx_dt] + 2*matr_K[2][d_Vx_dt] + matr_K[3][d_Vx_dt]);
        vec_ASK_temp.Vy   = vec_ASK.Vy   + dt*(1.0/6.0)*(matr_K[0][d_Vy_dt] + 2*matr_K[1][d_Vy_dt] + 2*matr_K[2][d_Vy_dt] + matr_K[3][d_Vy_dt]);
        vec_ASK_temp.x    = vec_ASK.x    + dt*(1.0/6.0)*(matr_K[0][d_x_dt ] + 2*matr_K[1][d_x_dt ] + 2*matr_K[2][d_x_dt ] + matr_K[3][d_x_dt ]);
        vec_ASK_temp.y    = vec_ASK.y    + dt*(1.0/6.0)*(matr_K[0][d_y_dt ] + 2*matr_K[1][d_y_dt ] + 2*matr_K[2][d_y_dt ] + matr_K[3][d_y_dt ]);
        vec_ASK_temp.m    = vec_ASK.m    + dt*(1.0/6.0)*(matr_K[0][d_m_dt ] + 2*matr_K[1][d_m_dt ] + 2*matr_K[2][d_m_dt ] + matr_K[3][d_m_dt ]);
        vec_ASK_temp.Time = vec_ASK.Time + dt;


        h_izm = f_h_izm(f_V(vec_ASK_temp.Vx, vec_ASK_temp.Vy), f_r(vec_ASK_temp.x, vec_ASK_temp.y));

        if (h_izm>=config.hf ) {
            if (fabs(h_izm-config.hf)>config.eps)
            {
                EoS = false;
                dt /= 10;
            } else {
                EoR = true;
            }

        } else {
            EoS = true;
        }


        if (shag>=1E5) {
            EoR = true;
        }
        if (EoS) {
            vec_ASK = vec_ASK_temp;
            vec_KE = AGESK_to_KE(vec_ASK);


            StepData result;
            result.alpha = (f_TETA(vec_ASK.Vx, vec_ASK.Vy)-f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK.Time)); //*toDeg
            result.TETA = f_TETA(vec_ASK.Vx, vec_ASK.Vy);
            result.gamma = f_gamma(upr[u_gamma], upr[u_dgdt], vec_ASK.Time);
            result.gt = upr[u_dgdt];
            result.globalTime = vec_ASK.Time;
            result.vASK = vec_ASK;
            result.vKE = vec_KE;
          //  result.DT = data.DT_0.addMSecs(1000*result.time);

            resultData.append(result);


            shag++;
        } else {
            vec_ASK_temp = vec_ASK;
        }
    }
}

double modeling_flight_2D::propUpr(Vector U) {
    upr = U; // вектор управляющих параметров
    propagate();
    return resultData.last().vASK.Time;
}

QString StepData::Str_All_Header(QString del) {
    return  "t"     + del +
            "x"     + del +
            "y"     + del +
            "z"     + del +
            "r"     + del +
            "Vx"    + del +
            "Vy"    + del +
            "Vz"    + del +
            "V"     + del +
            "TETA"  + del +
            "alpha" + del +
            "gamma" + del +
            "gt"    + del +
            "a"     + del +
            "e"     + del +
            "i"     + del +
            "RAAN"  + del +
            "omega" + del +
            "u"     + del +
            "m"     + del +
            "beta";
}

QString StepData::toStr_All(QString del) {
    return  QString::number(vASK.Time)      + del +
            QString::number(vASK.x)         + del +
            QString::number(vASK.y)         + del +
            QString::number(vASK.z)         + del +
            QString::number(vASK.r())       + del +
            QString::number(vASK.Vx)        + del +
            QString::number(vASK.Vy)        + del +
            QString::number(vASK.Vz)        + del +
            QString::number(vASK.V())       + del +
            QString::number(TETA*toDeg)     + del +
            QString::number(alpha*toDeg)    + del +
            QString::number(gamma*toDeg)    + del +
            QString::number(gt*toDeg)       + del +
            QString::number(vKE.a)          + del +
            QString::number(vKE.e)          + del +
            QString::number(vKE.i*toDeg)    + del +
            QString::number(vKE.RAAN*toDeg) + del +
            QString::number(vKE.om*toDeg)   + del +
            QString::number(vKE.u*toDeg)    + del +
            QString::number(vASK.m)         + del +
            QString::number(beta);
}

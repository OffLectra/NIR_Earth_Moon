#include "modeling_flight_rk4.h"



void modeling_flight_RK4::propagate_unperturbed(){
    double m_K[4][d_LAST] = {0.0};
    ASK_param_vec v_ASK = data.v_ASK, v_ASK_t = v_ASK;  //1)- явл.пер для хран век. сост. в момент оконч расч шага(y_i, y_i+1, ...);
    //2) явл врем перем для хран век сост на подшагах инт ((y, y1, y2, y3)_i, (y, y1, y2, y3)_i+1, ...);
    Kep_param_vec v_KE = data.v_KE;

    calcData.append(StepData(v_ASK,v_KE));

    int shag = 0;

    bool EoR = false;
    bool EoS = true;


    while (!EoR) {    //1 итер - 1 шаг инт

        // 1)------------------------
        m_K[0][d_Vx_dt] = f_d_Vx_dt(f_gx(v_ASK_t.x, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[0][d_Vy_dt] = f_d_Vy_dt(f_gy(v_ASK_t.y, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[0][d_Vz_dt] = f_d_Vz_dt(f_gz(v_ASK_t.z, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[0][d_x_dt ] = f_d_x_dt(v_ASK_t.Vx);
        m_K[0][d_y_dt ] = f_d_y_dt(v_ASK_t.Vy);
        m_K[0][d_z_dt ] = f_d_z_dt(v_ASK_t.Vz);
        m_K[0][d_m_dt ] = f_d_m_dt(data.beta);

        v_ASK_t.Vx = v_ASK.Vx + config.dt*m_K[0][d_Vx_dt]*0.5;
        v_ASK_t.Vy = v_ASK.Vy + config.dt*m_K[0][d_Vy_dt]*0.5;
        v_ASK_t.Vz = v_ASK.Vz + config.dt*m_K[0][d_Vz_dt]*0.5;
        v_ASK_t.x  = v_ASK.x  + config.dt*m_K[0][d_x_dt]*0.5;
        v_ASK_t.y  = v_ASK.y  + config.dt*m_K[0][d_y_dt]*0.5;
        v_ASK_t.z  = v_ASK.z  + config.dt*m_K[0][d_z_dt]*0.5;
        v_ASK_t.m  = v_ASK.m  + config.dt*m_K[0][d_m_dt]*0.5;
        v_ASK_t.Time = v_ASK.Time + config.dt * 0.5;

        // 2)------------------------
        m_K[1][d_Vx_dt] = f_d_Vx_dt(f_gx(v_ASK_t.x, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[1][d_Vy_dt] = f_d_Vy_dt(f_gy(v_ASK_t.y, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[1][d_Vz_dt] = f_d_Vz_dt(f_gz(v_ASK_t.z, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[1][d_x_dt ] = f_d_x_dt(v_ASK_t.Vx);
        m_K[1][d_y_dt ] = f_d_y_dt(v_ASK_t.Vy);
        m_K[1][d_z_dt ] = f_d_z_dt(v_ASK_t.Vz);
        m_K[1][d_m_dt ] = f_d_m_dt(data.beta);

        v_ASK_t.Vx = v_ASK.Vx + config.dt * m_K[1][d_Vx_dt]*0.5;
        v_ASK_t.Vy = v_ASK.Vy + config.dt * m_K[1][d_Vy_dt]*0.5;
        v_ASK_t.Vz = v_ASK.Vz + config.dt * m_K[1][d_Vz_dt]*0.5;
        v_ASK_t.x  = v_ASK.x  + config.dt * m_K[1][d_x_dt]*0.5;
        v_ASK_t.y  = v_ASK.y  + config.dt * m_K[1][d_y_dt]*0.5;
        v_ASK_t.z  = v_ASK.z  + config.dt * m_K[1][d_z_dt]*0.5;
        v_ASK_t.m  = v_ASK.m  + config.dt * m_K[1][d_m_dt]*0.5;
        v_ASK_t.Time = v_ASK.Time + config.dt * 0.5;

        // 3)------------------------
        m_K[2][d_Vx_dt] = f_d_Vx_dt(f_gx(v_ASK_t.x, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[2][d_Vy_dt] = f_d_Vy_dt(f_gy(v_ASK_t.y, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[2][d_Vz_dt] = f_d_Vy_dt(f_gz(v_ASK_t.z, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[2][d_x_dt ] = f_d_x_dt(v_ASK_t.Vx);
        m_K[2][d_y_dt ] = f_d_y_dt(v_ASK_t.Vy);
        m_K[2][d_z_dt ] = f_d_z_dt(v_ASK_t.Vz);
        m_K[2][d_m_dt ] = f_d_m_dt(data.beta);

        v_ASK_t.Vx = v_ASK.Vx + config.dt * m_K[2][d_Vx_dt];
        v_ASK_t.Vy = v_ASK.Vy + config.dt * m_K[2][d_Vy_dt];
        v_ASK_t.Vz = v_ASK.Vz + config.dt * m_K[2][d_Vz_dt];
        v_ASK_t.x  = v_ASK.x  + config.dt * m_K[2][d_x_dt];
        v_ASK_t.y  = v_ASK.y  + config.dt * m_K[2][d_y_dt];
        v_ASK_t.z  = v_ASK.z  + config.dt * m_K[2][d_z_dt];
        v_ASK_t.m  = v_ASK.m  + config.dt * m_K[2][d_m_dt];
        v_ASK_t.Time = v_ASK.Time + config.dt;

        // 4)------------------------
        m_K[3][d_Vx_dt] = f_d_Vx_dt(f_gx(v_ASK_t.x, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[3][d_Vy_dt] = f_d_Vy_dt(f_gy(v_ASK_t.y, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[3][d_Vz_dt] = f_d_Vz_dt(f_gz(v_ASK_t.z, f_r(v_ASK_t.x, v_ASK_t.y, v_ASK_t.z)));
        m_K[3][d_x_dt ] = f_d_x_dt(v_ASK_t.Vx);
        m_K[3][d_y_dt ] = f_d_y_dt(v_ASK_t.Vy);
        m_K[3][d_z_dt ] = f_d_z_dt(v_ASK_t.Vz);
        m_K[3][d_m_dt ] = f_d_m_dt(data.beta);

        //---------------------------
        v_ASK_t.Vx   = v_ASK.Vx   + config.dt * (1.0/6.0)*(m_K[0][d_Vx_dt] + 2*m_K[1][d_Vx_dt] + 2*m_K[2][d_Vx_dt] + m_K[3][d_Vx_dt]);
        v_ASK_t.Vy   = v_ASK.Vy   + config.dt * (1.0/6.0)*(m_K[0][d_Vy_dt] + 2*m_K[1][d_Vy_dt] + 2*m_K[2][d_Vy_dt] + m_K[3][d_Vy_dt]);
        v_ASK_t.Vz   = v_ASK.Vz   + config.dt * (1.0/6.0)*(m_K[0][d_Vz_dt] + 2*m_K[1][d_Vz_dt] + 2*m_K[2][d_Vz_dt] + m_K[3][d_Vz_dt]);
        v_ASK_t.x    = v_ASK.x    + config.dt * (1.0/6.0)*(m_K[0][d_x_dt ] + 2*m_K[1][d_x_dt ] + 2*m_K[2][d_x_dt ] + m_K[3][d_x_dt ]);
        v_ASK_t.y    = v_ASK.y    + config.dt * (1.0/6.0)*(m_K[0][d_y_dt ] + 2*m_K[1][d_y_dt ] + 2*m_K[2][d_y_dt ] + m_K[3][d_y_dt ]);
        v_ASK_t.z    = v_ASK.z    + config.dt * (1.0/6.0)*(m_K[0][d_z_dt ] + 2*m_K[1][d_z_dt ] + 2*m_K[2][d_z_dt ] + m_K[3][d_z_dt ]);
        v_ASK_t.m    = v_ASK.m    + config.dt * (1.0/6.0)*(m_K[0][d_m_dt ] + 2*m_K[1][d_m_dt ] + 2*m_K[2][d_m_dt ] + m_K[3][d_m_dt ]);
        v_ASK_t.Time = v_ASK.Time + config.dt;

        v_KE = AGESK_to_KE(v_ASK_t); //функция перевода из агэск в кэ

        //проверка выхода по положению на орбите
        bool isEqual = angle_compare(v_KE.om, v_KE.u, 3*toRad);
        if (isEqual && v_KE.u>v_KE.om) {
            if (!angle_compare(v_KE.om, v_KE.u, 0.01*toRad)) {
                EoS = false;
                config.dt /= 10.0;
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

#include "atm_gost_.h"
#include "koef.h"

using namespace std;

const double ro_0 = 1.58868e-08; // кг/ м^3, плотность ночной атмосферы на высоте 120км
const double w_z = 7.292115e-05; // рад/с, угловая скорость вращения Земли

const double dt_f107 = 1.7;		 // запаздывание по F107
const double dt_Kp = 0.6;		 // запаздывание по Kp
const double dt_kpp = 0.25;  	 // сут, запаздывание по kpp


//  ro_n = ro_0*exp(a_0 + a_1*h + a_2*pow(h,2) + a_3*pow(h,3) + a_4*pow(h,4) + a_5*pow(h,5) + a_6*pow(h,6));

//  K_4 = (e_0 + e_1*h + e_2*pow(h,2) + e_3*pow(h, 3) + e_4*pow(h,4))*(e_5 + e_6*K_p + e_7*pow(K_p,2) + e_8*pow(K_p, 3));
//  K_3 = (b_0 + b_1*h + b_2*pow(h,2) + b_3*pow(h,3) + b_4*pow(h, 4))*(F_10_7 - F_81)/(F_81 + abs(F_10_7 - F_81));
//  A_ = A_0 + A_1*d + A_2*pow(d,2) + A_3*pow(d,3) + A_4*pow(d,4) + A_5*pow(d,5) + A_6*pow(d,6) + A_7*pow(d,7) + A_8*pow(d,8);
//  K_2 = A*(d_0 + d_1*h + d_2*pow(h,2) + d_3*pow(h,3) + d_4*pow(h,4));
//  beta = alpha - S_zv - w*t + fi_1;
//  fi = acos(1/r*(z*sin(delta)+cos(delta)*(x*cos(beta)+y*sin(beta)))); ///???????
//  K_1 = (c_0 + c_1*h + c_2*pow(h,2) + c_3*pow(h,3) + c_4*pow(h,4))*pow(cos(fi/2), (n_0+n_1*h+n_2*h*h));
//  K_0 = 1 + (l_0 + l_1*h + l_2*pow(h,2) + l_3*pow(h,3) + l_4*pow(h,4))*(F_81-F_0)/F_0;
//  ro = ro_n*K_0*(1+K_1+K_2+K_3+K_4);

    double calcF81(double const F107[81]) //значения F107 начиная с 80 суток перед текущей датой до текущей даты
    {
        /* формула по расчёту F81 из ГОСТ */
        double num = 0.0; // числитель
        double denom = 0.0; // знаменатель
        for (int i = -80; i < 1; i++)
        {
            double Wi = 1.0 + (0.5 * i) / 80.0;
            num += F107[i + 80] * Wi;
            denom += Wi;
        }
        double F81 = num / denom;
        return F81;
    }

    /*
        [RU] Функция atmosGOST_R_25645_166_2004() для расчёта плотности верхней земной атмосферы
        по модели из российского стандарта ГОСТ Р 25645.166-2004. В случае успешного расчёта возвращает
        плотность в кг/куб.м.
    */
    double atmosGOST_R_25645_166_2004(
        double const h_km,	    // altitude above Earth ellipsoid 120<h_km<1500, km / высота над уровнем земного эллипсоида, от 120 до 1500, км
        double const F107,	    // F10.7 solar emission index / индекс солнечной активности
        double const Kp,	    // квазилогарифмический планетарный среднесуточный индекс геомагнитной активности, баллы
        double const F81,	    // averaged weighted F10.7 for previous 80 days + current day / усреднённый за 81 сутки (80 предыдущих + 1 текущие) и взвешенный индекс солнечной активности
        double const DoY,	    // number of day from the beginning of the year / номер суток от начала года //TODO заменить на Mjd_TT
        double const X[3],	    // x, y, z - geocentric greenwich coordinates, km / гринвичские координаты точи пространства, км
        double const t_s,	    // всемирное время, с
        double const S_rad,	    // sidereal midnight time, rad / звёздное время в гринвическую полночь, рад
        double const alpha_rad,	// right ascention of the Sun, rad / прямое восхождение Солнца, рад
        double const delta_rad)	// declination of the Sun, rad / склонение Солнца, рад
    {
//        const double a_c = -0.004206801077112602;
//        const double b_c = 3.0743437315404796;

//        const double c_c = pow(fabs(a_c),b_c);

        // 1. Предварительный расчёт степеней от h_km, DoY, Kp:
        // h_vec = {1, h, h^2, h^3, ..., h^6}
        double h_vec[7];
        h_vec[0] = 1.0;
            for (int i = 1; i < 7; i++)
            {
                h_vec[i] = h_vec[i - 1] * h_km;
            }

        // doy_vec = {1, DoY, DoY^2, DoY^3, ..., DoY^8}
        double doy_vec[9] = { 0.0 };
        doy_vec[0] = 1.0;
            for (int i = 1; i < 9; i++)
            {
                doy_vec[i] = doy_vec[i - 1] * DoY;
            }

        double Kp_vec[4] = { 0.0 };
        Kp_vec[0] = 1.0;
        for (int i = 1; i < 4; i++)
            {
                Kp_vec[i] = Kp_vec[i - 1] * Kp;
            }

        // 2. Выбор опорного значения индекса солнечной активности F0
        // и номера колонки для таблицы 2 или 3

        int n_col = 0;
        double diff1 = abs(F0_arr[n_col] - F81);
        double diff2 = diff1;
            while (n_col < 6)
            {
                n_col++;
                diff2 = abs(F0_arr[n_col] - F81);
                    if (diff1 <= diff2)
                    {
                        --n_col;
                        break;
                    }
                diff1 = diff2;
            }


        // 3. Выбор коэффициентов из таблиц, учитывая высоту (2 либо 3 таблицы ГОСТа)
        // в том порядке, в котором они приведены в таблицах 2-3
        double a[7] = { 0.0 };
        if (h_km > a_high_table[0][n_col])
        {
            for (int i = 0; i < 7; i++)
            {
                a[i] = a_high_table[i + 1][n_col];
            }
        }
        else {
            for (int i = 0; i < 7; i++)
            {
                a[i] = a_low_table[i + 1][n_col];
            }
        }

        double b[5] = { 0.0 };
        if (h_km > b_high_table[0][n_col])
        {
            for (int i = 0; i < 5; i++)
            {
                b[i] = b_high_table[i + 1][n_col];
            }
        }
        else {
            for (int i = 0; i < 5; i++)
            {
                b[i] = b_low_table[i + 1][n_col];
            }
        }

        double c[5] = { 0.0 };
        if (h_km > c_high_table[0][n_col])
        {
            for (int i = 0; i < 5; i++)
            {
                c[i] = c_high_table[i + 1][n_col];
            }
        }
        else {
            for (int i = 0; i < 5; i++)
            {
                c[i] = c_low_table[i + 1][n_col];
            }
        }

        double n[3] = { 0.0 };
        if (h_km > c_high_table[0][n_col])
        { // высота именно из таблицы коэффициента c
            for (int i = 0; i < 3; i++)
            {
                n[i] = n_high_table[i][n_col];
            }
        }
        else {
            for (int i = 0; i < 3; i++)
            {
                n[i] = n_low_table[i][n_col];
            }
        }

        // коэффициент модели, равный углу запаздывания максимума плотности по отношению к максимуму освещённости, рад
        double phi1 = 0.0;
        if (h_km > c_high_table[0][n_col])
        { // высота именно из таблицы коэффициента c
            phi1 = phi1_high_table[n_col];
        }
        else {
            phi1 = phi1_low_table[n_col];
        }

        double d[5] = { 0.0 };
        if (h_km > d_high_table[0][n_col])
        {
            for (int i = 0; i < 5; i++)
            {
                d[i] = d_high_table[i + 1][n_col];
            }
        }
        else {
            for (int i = 0; i < 5; i++)
            {
                d[i] = d_low_table[i + 1][n_col];
            }
        }

        double e[9] = { 0.0 };
        if (h_km > e_high_table[0][n_col])
        {
            for (int i = 0; i < 9; i++)
            {
                e[i] = e_high_table[i + 1][n_col];
            }
        }
        else {
            for (int i = 0; i < 9; i++)
            {
                e[i] = e_low_table[i + 1][n_col];
            }
        }

        double et[4] = { 0.0 };
        if (h_km > e_high_table[0][n_col]) { // высота именно из таблицы коэффициента e
            for (int i = 0; i < 4; i++)
            {
                et[i] = et_high_table[i][n_col];
            }
        }
        else {
            for (int i = 0; i < 4; i++)
            {
                et[i] = et_low_table[i][n_col];
            }
        }

        double l[5] = { 0.0 };
        if (h_km > l_high_table[0][n_col])
        {
            for (int i = 0; i < 5; i++)
            {
                l[i] = l_high_table[i + 1][n_col];
            }
        }
        else {
            for (int i = 0; i < 5; i++)
            {
                l[i] = l_low_table[i + 1][n_col];
            }
        }

        // 4. Расчёт коэффициентов K0-K4
        // 4.1 коэффициент K0, учитывающий изменение плотности атмосферы, связанное с отклонением
        double K0 = 0.0;
        for (int i = 0; i < 5; i++) // сборка полинома
        {
            K0 += l[i] * h_vec[i];
        }
//        test_glob_K0 = K0;

        K0 *= (F81 - F0_arr[n_col]) / F0_arr[n_col] + 1.0; // коэффициент в связи с отклонением F81 от F0

        // 4.2 коффициент K1, учитывающий суточный эффект в распределении плотности
        double K1 = 0.0; // (c' * h_vec(1:5) ) * (cos_phi ^ ( n * h_vec(1:3) )) / 2; // суточный коэффициент распределения плотности
        for (int i = 0; i < 5; i++) // сборка полинома
        {
            K1 += c[i] * h_vec[i];
        }
//        test_glob_K1 = K1;

        double cos_power = 0.0;
        for (int i = 0; i < 3; i++) // сборка полинома
        {
            cos_power += n[i] * h_vec[i];
        };

        double r = sqrt(X[0]* X[0] + X[1]*X[1] + X[2]*X[2]); // расстояние от центра гривничской СК

        // разность между долготой, для которой рассчитывают плотность атмосферы
        // и долготой с максимальным значением плотности в её суточном распределении, рад
        double beta_rad = alpha_rad - S_rad - w_z * t_s + phi1;

        double cos_phi = 1 / r * (X[2] * sin(delta_rad)
            + cos(delta_rad) * (X[0] * cos(beta_rad) + X[1] * sin(beta_rad)));
        double cos_phi_pow = pow(fabs(cos_phi), cos_power);

        K1 *= copysign(cos_phi_pow/2,cos_phi);

        // 4.3 коэффициент K2, учитывающий полугодовой эффект
        double K2 = 0.0;
        for (int i = 0; i < 5; i++) // сборка полинома
        {
            K2 += d[i] * h_vec[i];
        }
//        test_glob_K2 = K2;

        double Ad = 0.0;
        for (int i = 0; i < 9; i++) // сборка полинома
        {
            Ad += A[i] * doy_vec[i]; // полугодовой коэффициент
        }
        K2 *= Ad;

        // 4.4 коэффициент K3, учитывающий изменение плотности, связанное с отклонением F107 от F81
        double K3 = 0.0;
        for (int i = 0; i < 5; i++) // сборка полинома
        {
            K3 += b[i] * h_vec[i];
        }
//        test_glob_K3 = K3;

        K3 *= (F107 - F81) / (F81 + abs(F107 - F81)); // изменение плотности в связи с отклонением F10.7 от F81

        // 4.5 Расчёт коэффициента K4, учитывающий зависимость плотности атмосферы от геомагнитной возмущенности
        // при использовании среднесуточных коэффициентов геомагнитной активности
        // e[5]...e[8]  и  Kp, при использовании 3-х часовых - et[5]...et[8] и kpp
        double K4_1 = 0.0; // первый множитель K4
        for (int i = 0; i < 5; i++) // сборка полинома
        {
            K4_1 += e[i] * h_vec[i];
        }

        double K4_2 = 0.0; // второй множитель K4
        for (int i = 0; i < 4; i++) // сборка полинома
        {
            K4_2 += e[i+5] * Kp_vec[i];
        }
//        test_glob_K4_1 = K4_1;
//        test_glob_K4_2 = K4_2;

        double K4 = K4_1 * K4_2;

        // 4.6 финальная формула
        double polynom = 0.0;
        for (int i = 0; i < 7; i++) // сборка полинома
        {
            polynom += a[i] * h_vec[i];
        }
        double rho_night = ro_0 * exp(polynom);
//        test_glob_rho_night = rho_night;

        // Определение плотности атмосферы
        return rho_night * K0 * (1 + K1 + K2 + K3 + K4);
    }


    double F0_arr_test[7] = { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };

    double rand_f(double a, double b) {
        return ((double)rand() / RAND_MAX) * (b - a) + a;
    }

//    bool test_by_tables4_9(
//        double const tol_rho,
//        double const tol_coef,
//        bool const verbose)
//    {
//        if (verbose) {
//            std::cout << std::endl << "=== Testing by matching GOST R 25645.166-2004 Tables 4-9 ==="
//                << std::endl;
//            std::cout << "Result" << "\th,km" << "\t\tF107"
//                << "\t\trho erros"
//                << "\t\t\t\tK0-K4' erros"
//                << std::endl;
//        }
//        double max_rho_err = 0.0;
//        double max_coef_err = 0.0;
//        double max_rho_err_i = 0.0;
//        double max_coef_err_i = 0.0;
//        for (int i = 0; i < 7; i++) // ������������ �� ����� F
//        {
//            for (int j = 0; j < 70; j++) // ������������ �� ����� ������
//            {
//                max_rho_err_i = 0.0;
//                max_coef_err_i = 0.0;
//                double h = H0 + j * H_step;
//                double X0[3] = { 0.0 };
//                X0[0] = h;

//                double rho = atmosGOST_R_25645_166_2004(
//                    h,
//                    F0_arr_test[i],
//                    1.0,
//                    F0_arr_test[i], //
//                    0.0,
//                    X0,
//                    0.0,
//                    0.0,
//                    0.0,
//                    0.0);

//                double err_rho0 = abs(test_glob_rho_night - table_rho_night[j][i]);
//                max_rho_err_i = std::fmax(max_rho_err_i, err_rho0);

//                double err_rhoK0 = abs(test_glob_K0 - table_K0[j][i]);
//                double err_rhoK1 = abs(test_glob_K1 - table_K1[j][i]);
//                double err_rhoK2 = abs(test_glob_K2 - table_K2[j][i]);
//                double err_rhoK3 = abs(test_glob_K3 - table_K3[j][i]);
//                double err_rhoK4_1 = abs(test_glob_K4_1 - table_K4[j][i]);

//                max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK0);
//                max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK1);
//                max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK2);
//                max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK3);
//                max_coef_err_i = std::fmax(max_coef_err_i, err_rhoK4_1);

//                if (verbose) {
//                    std::string rslt =
//                        ((max_rho_err_i < tol_rho)&&(max_coef_err_i < tol_coef))
//                        ? ("OK")
//                        : ("**ERR");
//                    std::cout << rslt << "\th = " << h
//                        << ";\tF107 = " << F0_arr_test[i]
//                        /*<< ";\tmax_rho_error = " << max_rho_err_i
//                        << ";\t\tmax_coef_error = " << max_coef_err_i
//                        << std::endl*/
//                        << "\t" /*<< table_rho_night[j][i] << "/"*/ << test_glob_rho_night
//                        << "\t" /*<< table_K0[j][i] << "/"*/ << test_glob_K0
//                        << "\t" /*<< table_K1[j][i] << "/"*/ << test_glob_K1
//                        << "\t" /*<< table_K2[j][i] << "/"*/ << test_glob_K2
//                        << "\t" /*<< table_K3[j][i] << "/"*/ << test_glob_K3
//                        << "\t" /*<< table_K4[j][i] << "/"*/ << test_glob_K4_1
//                        << std::endl;
//                }
//                max_rho_err = std::fmax(max_rho_err, max_rho_err_i);
//                max_coef_err = std::fmax(max_coef_err, max_coef_err_i);
//            }
//        }

//        return (max_rho_err <= tol_rho)&&(max_coef_err <= tol_coef);
//    }

//    bool test_by_tables10_11(
//        double const tol_coef,
//        bool const verbose)
//    {
//        if (verbose) {
//            std::cout << std::endl << "=== Testing by matching GOST R 25645.166-2004 Table 10 ==="
//                << std::endl;
//            std::cout << "Result" << "\th,km" << "\t\tF107"
//                << "\t\tK4'' error" << std::endl;
//        }
//        double max_coef_err = 0.0, err_rhoK4_2 = 0.0;
//        for (int i = 0; i < 7; i++) // F
//        {
//            for (int j = 0; j < 22; j++) // Kp
//            {
//                double Kp = Kp0 + j * Kp_step;
//                double h = rand_f(120.0, 1500.0);
//                double X0[3] = { 0.0 };
//                X0[0] = h;

//                double rho = atmosGOST_R_25645_166_2004(
//                    h,
//                    F0_arr_test[i],
//                    Kp,
//                    F0_arr_test[i],
//                    0.0,
//                    X0,
//                    0.0,
//                    0.0,
//                    0.0,
//                    0.0);

//                err_rhoK4_2 = abs(test_glob_K4_2 - table_K4_2_24h[j][i]);

//                if (verbose) {
//                    std::string rslt =
//                        (err_rhoK4_2 < tol_coef)
//                        ? ("OK")
//                        : ("**ERR");
//                    std::cout << rslt
//                        << "\th = " << h
//                        << ";\tF107 = " << F0_arr_test[i]
//                        << ";\terr_rhoK4_2 = " << err_rhoK4_2
//                        << std::endl;
//                }
//            }
//            max_coef_err = std::fmax(max_coef_err, err_rhoK4_2);
//        }
//        if (verbose) {
//            std::cout << "\tMaximal K4'' error = " << max_coef_err << std::endl;
//        }
//        return (max_coef_err <= tol_coef);
//    }

//int main()
//{
//    // TODO ����� �� �������

//    double tol_rho = 1e-10, tol_coef = 1.0e-03;
//    test_by_tables4_9(tol_rho, tol_coef, true);
//    test_by_tables10_11(tol_coef, true);

//    std::system("pause");
//    return 0;
//}



























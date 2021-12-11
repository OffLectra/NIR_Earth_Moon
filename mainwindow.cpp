#include <ui_mainwindow.h>
#include <QMessageBox>
#include <QDebug>
#include <QVector>
#include <QApplication>
#include "mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->on_pushButton_KE_to_AGESK_clicked();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_KE_to_AGESK_clicked() {
    Kep_param_vec v_KE;

    v_KE.a = ui->lineEdit_KE_a->text().toDouble()  ;
    v_KE.e = ui->lineEdit_KE_e->text().toDouble();
    v_KE.i = ui->lineEdit_KE_i->text().toDouble()   *toRad;
    v_KE.RAAN = ui->lineEdit_KE_RAAN->text().toDouble()*toRad;
    v_KE.om = ui->lineEdit_KE_omega->text().toDouble()  *toRad;
    v_KE.u = ui->lineEdit_KE_u->text().toDouble()   *toRad;

    ASK_param_vec v_AGESK;
    v_AGESK = KE_to_AGESK(v_KE);

    ui->lineEdit_ASK_x0->setText(QString::number(v_AGESK.x));
    ui->lineEdit_ASK_y0->setText(QString::number( v_AGESK.y));
    ui->lineEdit_ASK_z0->setText(QString::number( v_AGESK.z));
    ui->lineEdit_ASK_Vx0->setText(QString::number(v_AGESK.Vx));
    ui->lineEdit_ASK_Vy0->setText(QString::number(v_AGESK.Vy));
    ui->lineEdit_ASK_Vz0->setText(QString::number(v_AGESK.Vz));

}

void MainWindow::on_pushButton_AGESK_to_KE_clicked() {
    ASK_param_vec v_AGESK;

    v_AGESK.x = ui->lineEdit_ASK_x0->  text().toDouble();
    v_AGESK.y = ui->lineEdit_ASK_y0->  text().toDouble();
    v_AGESK.z = ui->lineEdit_ASK_z0->  text().toDouble();
    v_AGESK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    v_AGESK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    v_AGESK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();

    Kep_param_vec v_KE;
    v_KE = AGESK_to_KE(v_AGESK);

    ui->lineEdit_KE_a->setText(QString::number(  v_KE.a));
    ui->lineEdit_KE_e->setText(QString::number(v_KE.e));
    ui->lineEdit_KE_i->setText(QString::number(v_KE.i_d));
    ui->lineEdit_KE_RAAN->setText(QString::number(v_KE.RAAN_d));
    ui->lineEdit_KE_omega->setText(QString::number(v_KE.om_d));
    ui->lineEdit_KE_u->setText(QString::number(v_KE.u_d));

}


//-----

//----------------------------------------------
//---------------...----------...---------------
//-------------.****.--------.****.-------------
//-----------.*******.------.*******.-----------
//---------.**_______**.---.**______**.---------
//--------.**__________**.**_________**.--------
//-------.**____________***___________**.-------
//------.**______________*_____________*.-------
//------.**____________________________**.------
//------.**________CCC___+___У__У______**.------
//-------.**_______C____+++___УУ______**.-------
//--------.**______CCC___+___У_______**.--------
//----------.**____________________**.----------
//------------.**_________________**.-----------
//--------------.**_____________**.-------------
//----------------.**_________**.---------------
//------------------.**______**.----------------
//--------------------.**_**.-------------------
//---------------------.***.--------------------
//----------------------.*.---------------------
//-----------------------.----------------------
//----------------------------------------------


void MainWindow::on_B_start_clicked() {
    NU_RK4 nu;
    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
    nu.P = ui->lineEdit_P->text().toDouble();
    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
    nu.beta = fabs(nu.P/(nu.P_ud*g0));
    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

    nu.v_KE = AGESK_to_KE(nu.v_ASK);

    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

    settings.dt = ui->lineEdit_dt->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());

    integr_RK4(nu, settings);

    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}

void MainWindow::on_B_start_Class_propagate_clicked() {
    NU_RK4 nu;
    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
    nu.P = ui->lineEdit_P->text().toDouble();
    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
    nu.beta = fabs(nu.P/(nu.P_ud*g0));
    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

    nu.v_KE = AGESK_to_KE(nu.v_ASK);

    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

    settings.dt = ui->lineEdit_dt->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


    modeling_flight_3D raschet(nu, settings);
    raschet.propagate();
    raschet.printCalcDataToFile(settings.file_name);

    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}

void MainWindow::on_B_Startatm_clicked() {
    NU_RK4 nu;
    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
    nu.P = ui->lineEdit_P->text().toDouble();
    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
    nu.beta = fabs(nu.P/(nu.P_ud*g0));
    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

    nu.v_KE = AGESK_to_KE(nu.v_ASK);

    nu.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;


    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

    settings.dt = ui->lineEdit_dt->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());

    integr_RK4_atm(nu, settings);



    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}

void MainWindow::on_B_Optima_clicked() {
    NU_RK4 nu_cur,nu_temp;
    nu_cur.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
    nu_cur.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
    nu_cur.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
    nu_cur.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    nu_cur.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    nu_cur.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
    nu_cur.v_ASK.m = ui->lineEdit_m->text().toDouble();
    nu_cur.P = ui->lineEdit_P->text().toDouble();
    nu_cur.P_ud = ui->lineEdit_P_ud->text().toDouble();
    nu_cur.beta = fabs(nu_cur.P/(nu_cur.P_ud*g0));
    nu_cur.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu_cur.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

    nu_cur.v_KE = AGESK_to_KE(nu_cur.v_ASK);

    nu_cur.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;


    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();
    settings.is_export_data_to_file = false;
    settings.dt = ui->lineEdit_dt->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());

    double eps = 1E-6;
    double t_cur = 0,t_temp=10;
    QVector<double> t_var = QVector<double>(2).fill(0.0);
    QVector<double> grad  = QVector<double>(2).fill(0.0);

    double dgamma = 1E-4;
    double ddgdt  = 1E-6;

    double lambda1 = 1E-5;
    double lambda2 = 2E-9;

    int counter = 1;
    while (fabs(t_cur-t_temp)>eps || fabs(nu_cur.gamma_0-nu_temp.gamma_0)>eps || fabs(nu_cur.d_gamma_dt-nu_temp.d_gamma_dt)>eps) {

        t_temp = t_cur;
        nu_temp = nu_cur;
        t_cur = integr_RK4(nu_cur, settings);

        nu_temp.gamma_0 = nu_temp.gamma_0 + dgamma;
        t_var[0] = integr_RK4(nu_temp, settings);

        nu_temp = nu_cur;

        nu_temp.d_gamma_dt = nu_temp.d_gamma_dt + ddgdt;
        t_var[1] = integr_RK4(nu_temp, settings);

        QString out = QString::number(counter) + " iteration\n";
        out += QString("t_cur") + "\t" + "t_var[0]" + "\t" + "t_var[1]" + "\n";
        out += QString::number(t_cur) + "\t" + QString::number(t_var[0]) + "\t" + QString::number(t_var[1]) + "\n";
        nu_temp = nu_cur;

        grad[0] = (t_var.at(0)-t_cur)/dgamma;
        grad[1] = (t_var.at(1)-t_cur)/ddgdt;
        out += QString("grad[0]") + "\t" + "t_grad[1]" + "\n";
        out += QString::number(grad[0]) + "\t" + QString::number(grad[1]) + "\n";

        nu_cur.gamma_0 -= lambda1* grad[0];
        nu_cur.d_gamma_dt -= lambda2* grad[1];
        out += QString("gamma_0") + "\t" + "d_gamma_dt" + "\n";
        out += QString::number(nu_cur.gamma_0*180/M_PI) + "\t" + QString::number(nu_cur.d_gamma_dt*180/M_PI) + "\n";
        ui->TB_vivod->append(out);

        qApp->processEvents();

        if (counter % 1000 == 0){
            lambda1 *= 2;
            lambda2 *= 2;
        }
        counter++;
    }


    settings.is_export_data_to_file = true;
    integr_RK4(nu_temp, settings);

}

//double Optima3(double a2, double a1, double a0, tExport *data)
//{
//    int i=1;
//    double x_cur = 100.0, x_pr = 0.0, da0=1E-4, da1=2E-4, da2=1E-5, step0 = 2*da0, step1 = 2*da1, step2 = 2*da2;
//    QVdouble a_cur = {a0, a1, a2};
//    QVdouble a_0 = {a0, a1, a2};
//    QVdouble grad_a = {0.0, 0.0, 0.0};
//    QVdouble x_cur_ARR = {0.0, 0.0, 0.0};
//    while (fabs(x_cur-x_pr) > 1E-6)
//    {

//        x_pr = x_cur;

//        x_cur = IntegratorDP_E1(a_cur[2], a_cur[1], a_cur[0]);
//        x_cur_ARR[0] = IntegratorDP_E1(a_cur[2], a_cur[1], a_cur[0] + da0);
//        x_cur_ARR[1] = IntegratorDP_E1(a_cur[2], a_cur[1]+da1, a_cur[0]);
//        x_cur_ARR[2] = IntegratorDP_E1(a_cur[2]+da2, a_cur[1], a_cur[0]);

//        grad_a[0]=(x_cur_ARR[0]-x_cur)/da0;
//        grad_a[1]=(x_cur_ARR[1]-x_cur)/da1;
//        grad_a[2]=(x_cur_ARR[2]-x_cur)/da2;


//        a_cur[0]+=copysign(step0,grad_a[0]);
//        a_cur[1]+=copysign(step1,grad_a[1]);
//        a_cur[2]+=copysign(step2,grad_a[2]);

//       i++;
//       if (i%10==0) {
//           step0/=2;
//           step1/=2;
//           step2/=2;
//       }
//    }
//    IntegratorDP_E1(a_cur[2], a_cur[1], a_cur[0], data);
//    return x_cur;
//}

//double Optima1(double a2, double a1, double a0, tExport *data)
//{
//    int i=1;
//    double x_cur = 100.0, x_pr = 0.0, da0=1E-3, da1=2E-4, da2=1E-5, step0 = 2*da0, step1 = 2*da1, step2 = 2*da2;
//    QVdouble a_cur = {a0, a1, a2};
//    QVdouble a_0 = {a0, a1, a2};
//    QVdouble grad_a = {0.0, 0.0, 0.0};
//    QVdouble x_cur_ARR = {0.0, 0.0, 0.0};
//    while (fabs(x_cur-x_pr) > 1E-6)
//    {

//        x_pr = x_cur;

//        x_cur = IntegratorDP_E1(a_cur[2], a_cur[1], a_cur[0]);
//        x_cur_ARR[0] = IntegratorDP_E1(a_cur[2], a_cur[1], a_cur[0] + da0);

//        grad_a[0]=(x_cur_ARR[0]-x_cur)/da0;

//        a_cur[0]+=copysign(step0,grad_a[0]);

//       i++;
//       if (i%10==0) {
//           step0/=2;

//       }
//    }
//    IntegratorDP_E1(a_cur[2], a_cur[1], a_cur[0], data);
//    return x_cur;
//}

void MainWindow::on_B_Opt_DFP_clicked() {
    NU_RK4 nu_cur;
    nu_cur.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
    nu_cur.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
    nu_cur.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
    nu_cur.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    nu_cur.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    nu_cur.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
    nu_cur.v_ASK.m = ui->lineEdit_m->text().toDouble();
    nu_cur.P = ui->lineEdit_P->text().toDouble();
    nu_cur.P_ud = ui->lineEdit_P_ud->text().toDouble();
    nu_cur.beta = fabs(nu_cur.P/(nu_cur.P_ud*g0));
    nu_cur.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu_cur.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

    nu_cur.v_KE = AGESK_to_KE(nu_cur.v_ASK);

    nu_cur.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;


    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();
    settings.is_export_data_to_file = false;
    settings.dt = ui->lineEdit_dt->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());

    double eps = 1E-6; // точность проверки условия выхода их расчета (нашли минимум)
    double f_x_k = 0, f_x_k_1=10;

    Vector x_0 = {nu_cur.gamma_0,nu_cur.d_gamma_dt}; //Вектор управляющих параметров (начальный)

    Vector x_k_1 = x_0;                 //Вектор управляющих параметров (к-1) начальная точка итерации
    Vector x_k = Vector(2).fill(0.0);   //Вектор управляющих параметров (к)   конечная  точка итерации

    Vector grad_x_k_1  = Vector(2).fill(2.0);   //Вектор градиентов в точке Х(к-1)
    Vector grad_x_k  = Vector(2).fill(0.0);     //Вектор градиентов в точке Х(к)


    double dgamma = 1E-4; // дельты необходимые для поиска градиентов в точках
    double ddgdt  = 1E-7;



    Vector deltax = {dgamma,ddgdt}; // дельты для удобства в векторной форме

    double lambda1 = 1;//1E-5;   // некие вспомогательные константы, понижающие влияние вектора градиента по направлениям
    double lambda2 = 1;//2E-10;
//    double lambda_mn = 1E-3;
    Vector dw = min_VxV(grad_x_k_1,grad_x_k); // вектор разницы градиентов в точках -grad_x_k + grad_x_k_1 (т.к. антиградиенты)
    Vector dx = min_VxV(x_k,x_k_1);           // вектор разницы значения координат точек Х(к) - Х(к-1) (т.к. антиградиенты)
    Vector dxt = Vector(2).fill(0.0); //dx~      вспомогательный вектор dx + А*dw

    Vector lam = {lambda1,lambda2};

//    Matr A = {{1.0,0.0},{0.0,1.0}};
    Matr A = {{1E-3,0.0},{0.0,2E-9}};     //Некая вспомогательная матрица метода ДФП, обычно принимают единичной, но необходимо парировать величины градиентов, поэтому немного другая

    int counter = 1;
    while (fabs(grad_x_k[u_gamma]-grad_x_k_1[u_gamma])>eps ||
           fabs(grad_x_k[u_dgdt]-grad_x_k_1[u_dgdt])>eps ||
           fabs(f_x_k-f_x_k_1)>eps ||
           fabs(x_k[u_gamma]-x_k_1[u_gamma])>eps ||
           fabs(x_k[u_dgdt]-x_k_1[u_dgdt])>eps) {
        // начало итерации поиска минимума функции
        f_x_k = f_x_k_1;

        f_x_k_1 = integr_RK4_upr(nu_cur, settings,x_k_1); //1.

//        nu_temp.gamma_0 = nu_temp.gamma_0 + dgamma;
//        f_dx[0] = integr_RK4(nu_temp, settings); //2.

//        nu_temp = nu_cur;

//        nu_temp.d_gamma_dt = nu_temp.d_gamma_dt + ddgdt; //3.
//        f_dx[1] = integr_RK4(nu_temp, settings);

        QString out = QString::number(counter) + " iteration\n";

        out += QString("A") + "\n";
        out += QString::number(A[0][0]) + "\t" + QString::number(A[0][1]) + "\n";
        out += QString::number(A[1][0]) + "\t" + QString::number(A[1][1]) + "\n";

        out += QString("gamma") + "\t" + "dgdt" + "\n";
        out += QString::number(x_k_1[u_gamma]*180/M_PI) + "\t" + QString::number(x_k_1[u_dgdt ]*180/M_PI) + "\n";

        out += QString("fx") + "\t\n";
        out += QString::number(f_x_k_1) +"\n";

        grad_x_k_1 = grad_RK4_upr(x_k_1,deltax,nu_cur,settings); // 1. Вычисляем Градиент в точке X(к-1), начальной точке итерации

        out += QString("grad[u_gamma]") + "\t" + "grad[u_dgdt]" + "\n";
        out += QString::number(grad_x_k_1[u_gamma]) + "\t" + QString::number(grad_x_k_1[u_dgdt]) + "\n";

        Vector p_k_1 = mp_MxV(A,grad_x_k_1);                    // 2. Умножение матрицы А на вектор градиента в точке X(к-1) -> Вектор изменения наших параметров
        Vector lp_k_1 = {-lambda1*p_k_1[u_gamma],-lambda2*p_k_1[u_dgdt]}; // 3. Расчет вектора приращения для поиска новой точки (минус потому что ищим минимум, антиградиент)

        out += QString("p[u_gamma]") + "\t" + "p[u_dgdt]" + "\n";
        out += QString::number(p_k_1[u_gamma],'f',10) + "\t" + QString::number(p_k_1[u_dgdt],'f',10) + "\n";
        out += QString("lp[u_gamma]") + "\t" + "lp[u_dgdt]" + "\n";
        out += QString::number(lp_k_1[u_gamma],'f',10) + "\t" + QString::number(lp_k_1[u_dgdt],'f',10) + "\n";

        x_k = add_VxV(x_k_1,lp_k_1); // 4. Находим новую (конечную) точку (итерации) X(к)
//        x_k[u_gamma] -= lam[u_gamma] * grad_x_k_1[u_gamma];
//        x_k[u_dgdt ] -= lam[u_dgdt ] * grad_x_k_1[u_dgdt ];
        out += QString("gamma") + "\t" + "dgdt" + "\n";
        out += QString::number(x_k[u_gamma]*180/M_PI,'f',10) + "\t" + QString::number(x_k[u_dgdt ]*180/M_PI,'f',10) + "\n";


        qApp->processEvents();

        grad_x_k = grad_RK4_upr(x_k,deltax,nu_cur,settings); // 5. Вычисляем Градиент в точке X(к), конечной точке итерации

        out += QString("grad[u_gamma]k") + "\t" + "grad[u_dgdt]k" + "\n";
        out += QString::number(grad_x_k[u_gamma]) + "\t" + QString::number(grad_x_k[u_dgdt]) + "\n";

        // 6. Вычисляем новую матрицу А каждые сколько-то интераций, изначально каждую

        if (counter % 1 == 0) {
            Matr dA = A;

            dw = min_VxV(grad_x_k_1,grad_x_k); // 7. Вычисляем
            out += QString("dw[u_gamma]k") + "\t" + "dw[u_dgdt]k" + "\n";
            out += QString::number(dw[u_gamma],'f',10.0) + "\t" + QString::number(dw[u_dgdt],'f',10) + "\n";

            dx = min_VxV(x_k,x_k_1);            // 8. Вычисляем


            out += QString("dx[u_gamma]") + "\t" + "dx[u_dgdt]" + "\n";
            out += QString::number(dx[u_gamma],'f',10) + "\t" + QString::number(dx[u_dgdt],'f',10) + "\n";

            dxt = add_VxV(dx,mp_MxV(A,dw));     // 9. Вычисляем  Здесь А еще не изменилась!!!

            out += QString("dxt[u_gamma]") + "\t" + "dxt[u_dgdt]" + "\n";
            out += QString::number(dxt[u_gamma],'f',10.0) + "\t" + QString::number(dxt[u_dgdt],'f',10) + "\n";

            dA = mp_cxM(1.0/smp_VxV(dw,dxt),mp_VxVT(dxt)); // 10. Вычисляем

            out += QString("dA") + "\n";
            out += QString::number(dA[0][0],'f',10) + "\t" + QString::number(dA[0][1],'f',10) + "\n";
            out += QString::number(dA[1][0],'f',10) + "\t" + QString::number(dA[1][1],'f',10) + "\n";

            A = min_MxM(A,dA); // 11. Вычисляем новую матрицу

        }


//        if (counter % 1 == 0){
//            lambda1 /= 2;
//            lambda2 /= 2;
//        }
//        if (counter % 2 == 0){
//            lambda1 = fabs(lambda_mn/grad_x_k[u_gamma]);
//            lambda2 = fabs(lambda_mn/grad_x_k[u_dgdt]);
//        }

//        if (counter % 50  == 0 && lambda_mn>1e-4){
//            lambda_mn /= 2;
//        }

        if (isnan(A[0][0]) || isinf(A[0][0])) { break; }

        x_k_1 = x_k; // Назначаем старой точкой новуюю для следующей итерации
        ui->TB_vivod->append(out);
        counter++;
        // конец итерации поиска минимума
    }


    settings.is_export_data_to_file = true;
    integr_RK4_upr(nu_cur, settings,x_k);
}

void MainWindow::on_B_Optima_new_clicked() {
    NU_RK4 nu_cur;
    nu_cur.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
    nu_cur.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
    nu_cur.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
    nu_cur.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    nu_cur.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    nu_cur.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
    nu_cur.v_ASK.m = ui->lineEdit_m->text().toDouble();
    nu_cur.P = ui->lineEdit_P->text().toDouble();
    nu_cur.P_ud = ui->lineEdit_P_ud->text().toDouble();
    nu_cur.beta = fabs(nu_cur.P/(nu_cur.P_ud*g0));
    nu_cur.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu_cur.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

    nu_cur.v_KE = AGESK_to_KE(nu_cur.v_ASK);

    nu_cur.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;


    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();
    settings.is_export_data_to_file = false;
    settings.dt = ui->lineEdit_dt->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


    double eps = 1E-6;
    double f_x_k = 0, f_x_k_1=10;

    Vector x_0 = {nu_cur.gamma_0,nu_cur.d_gamma_dt};

    Vector x_k_1 = x_0;
    Vector x_k = Vector(2).fill(0.0);



    Vector f_dx  = Vector(2).fill(0.0);
    Vector grad_x_k_1  = Vector(2).fill(0.0);
    Vector grad_x_k  = Vector(2).fill(0.0);

    double dgamma = 1E-4;
    double ddgdt  = 1E-6;

    Vector deltax = {dgamma,ddgdt};

    double lambda1 = 1E-5;
    double lambda2 = 2E-9;
    Vector dw = min_VxV(grad_x_k_1,grad_x_k);
    Vector dx = min_VxV(x_k,x_k_1);
    Vector dxt = Vector(2).fill(0.0); //dx~

    Vector lam = {lambda1,lambda2};

    Matr A = {{1.0,0.0},{0.0,1.0}};

    int counter = 1;
    while (fabs(f_x_k-f_x_k_1)>eps || fabs(x_k[u_gamma]-x_k_1[u_gamma])>eps || fabs(x_k[u_dgdt]-x_k_1[u_dgdt])>eps) {

        f_x_k = f_x_k_1;

        f_x_k_1 = integr_RK4_upr(nu_cur, settings,x_k_1); //1.

//        nu_temp.gamma_0 = nu_temp.gamma_0 + dgamma;
//        f_dx[0] = integr_RK4(nu_temp, settings); //2.

//        nu_temp = nu_cur;

//        nu_temp.d_gamma_dt = nu_temp.d_gamma_dt + ddgdt; //3.
//        f_dx[1] = integr_RK4(nu_temp, settings);


        QString out = QString::number(counter) + " iteration\n";
        out += QString("f") + "\t" + "f_dx[0]" + "\t" + "f_dx[1]" + "\n";
        out += QString::number(f_x_k) + "\t" + QString::number(f_dx[0]) + "\t" + QString::number(f_dx[1]) + "\n";

        grad_x_k_1 = grad_RK4_upr(x_k_1,deltax,nu_cur,settings);

        out += QString("grad[u_gamma]") + "\t" + "grad[u_dgdt]" + "\n";
        out += QString::number(grad_x_k_1[u_gamma]) + "\t" + QString::number(grad_x_k_1[u_dgdt]) + "\n";

        Vector p_k_1 = mp_MxV(A,grad_x_k_1); // Умножение матрицы на вектор
        Vector lp_k_1 = {-lambda1*p_k_1[u_gamma],-lambda2*p_k_1[u_dgdt]}; //приращение для поиска новой точки

        out += QString("p[u_gamma]") + "\t" + "p[u_dgdt]" + "\n";
        out += QString::number(grad_x_k_1[u_gamma]) + "\t" + QString::number(grad_x_k_1[u_dgdt]) + "\n";
        out += QString("lp[u_gamma]") + "\t" + "lp[u_dgdt]" + "\n";
        out += QString::number(grad_x_k_1[u_gamma]) + "\t" + QString::number(grad_x_k_1[u_dgdt]) + "\n";

        x_k = add_VxV(x_k_1,lp_k_1);
//        x_k[u_gamma] -= lam[u_gamma] * grad_x_k_1[u_gamma];
//        x_k[u_dgdt ] -= lam[u_dgdt ] * grad_x_k_1[u_dgdt ];
        out += QString("gamma") + "\t" + "dgdt" + "\n";
        out += QString::number(x_k[u_gamma]*180/M_PI) + "\t" + QString::number(x_k[u_dgdt ]*180/M_PI) + "\n";
        ui->TB_vivod->append(out);

        qApp->processEvents();

        if (counter % 500 == 0) {
            lambda1 *= 2;
            lambda2 *= 2;
        }

        x_k_1 = x_k;
        counter++;
    }
}



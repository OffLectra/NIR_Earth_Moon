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

    ui->Plot->addGraph();
    ui->Plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

MainWindow::~MainWindow()
{
    delete ui;
}

NU_RK4 MainWindow::prepareNUASK() {
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
    nu.Wist = fabs(nu.P_ud*g0);
    nu.beta = fabs(nu.P/(nu.P_ud*g0));
    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;
    nu.v_KE = AGESK_to_KE(nu.v_ASK);
    nu.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;
    return nu;
}

Settings_RK4 MainWindow::prepareSettings() {
    Settings_RK4 settings;
    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

    settings.dt = ui->lineEdit_dt->text().toDouble();    
    settings.af = ui->lineEdit_a_per->text().toDouble();
    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());
    return settings;
}

void MainWindow::on_pushButton_KE_to_AGESK_clicked() {
    cKE v_KE;

    v_KE.a = ui->lineEdit_KE_a->text().toDouble()  ;
    v_KE.e = ui->lineEdit_KE_e->text().toDouble();
    v_KE.i = ui->lineEdit_KE_i->text().toDouble()   *toRad;
    v_KE.RAAN = ui->lineEdit_KE_RAAN->text().toDouble()  *toRad;
    v_KE.om = ui->lineEdit_KE_omega->text().toDouble()  *toRad;
    v_KE.u = ui->lineEdit_KE_u->text().toDouble()   *toRad;

    cASK v_AGESK;
    v_AGESK = KE_to_AGESK(v_KE);

    ui->lineEdit_ASK_x0->setText(QString::number(v_AGESK.x));
    ui->lineEdit_ASK_y0->setText(QString::number( v_AGESK.y));
    ui->lineEdit_ASK_z0->setText(QString::number( v_AGESK.z));
    ui->lineEdit_ASK_Vx0->setText(QString::number(v_AGESK.Vx));
    ui->lineEdit_ASK_Vy0->setText(QString::number(v_AGESK.Vy));
    ui->lineEdit_ASK_Vz0->setText(QString::number(v_AGESK.Vz));
}

void MainWindow::on_pushButton_AGESK_to_KE_clicked() {
    cASK v_AGESK;

    v_AGESK.x = ui->lineEdit_ASK_x0->  text().toDouble();
    v_AGESK.y = ui->lineEdit_ASK_y0->  text().toDouble();
    v_AGESK.z = ui->lineEdit_ASK_z0->  text().toDouble();
    v_AGESK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
    v_AGESK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
    v_AGESK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();

    cKE v_KE;
    v_KE = AGESK_to_KE(v_AGESK);

    ui->lineEdit_KE_a->setText(QString::number(  v_KE.a));
    ui->lineEdit_KE_e->setText(QString::number(v_KE.e));
    ui->lineEdit_KE_i->setText(QString::number(v_KE.i_d));
    ui->lineEdit_KE_RAAN->setText(QString::number(v_KE.RAAN_d));
    ui->lineEdit_KE_omega->setText(QString::number(v_KE.om_d));
    ui->lineEdit_KE_u->setText(QString::number(v_KE.u_d));
}

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
//    NU_RK4 nu;
//    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu.P = ui->lineEdit_P->text().toDouble();
//    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu.beta = fabs(nu.P/(nu.P_ud*g0));
//    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

//    nu.v_KE = AGESK_to_KE(nu.v_ASK);

//    Settings_RK4 settings;
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());

    integr_RK4(prepareNUASK(), prepareSettings());

    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}

void MainWindow::on_B_start_Class_propagate_clicked() {
//    NU_RK4 nu;
//    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu.P = ui->lineEdit_P->text().toDouble();
//    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu.beta = fabs(nu.P/(nu.P_ud*g0));
//    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

//    nu.v_KE = AGESK_to_KE(nu.v_ASK);

//    Settings_RK4 settings;
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


    modeling_flight_3D raschet(prepareNUASK(), prepareSettings());
    raschet.propagate();
    raschet.printCalcDataToFile(prepareSettings().file_name);

    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}

void MainWindow::on_B_Startatm_clicked() {
//    NU_RK4 nu;
//    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu.P = ui->lineEdit_P->text().toDouble();
//    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu.beta = fabs(nu.P/(nu.P_ud*g0));
//    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

//    nu.v_KE = AGESK_to_KE(nu.v_ASK);

//    nu.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;

//    Settings_RK4 settings;
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());

    integr_RK4_atm(prepareNUASK(), prepareSettings());


    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}


// Решение градиентым спуском (В лоб)
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
    NU_RK4 nu_cur = prepareNUASK();
//    nu_cur.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu_cur.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu_cur.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu_cur.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu_cur.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu_cur.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu_cur.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu_cur.P = ui->lineEdit_P->text().toDouble();
//    nu_cur.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu_cur.beta = fabs(nu_cur.P/(nu_cur.P_ud*g0));
//    nu_cur.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu_cur.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;
//    nu_cur.v_KE = AGESK_to_KE(nu_cur.v_ASK);
//    nu_cur.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;

    Settings_RK4 settings = prepareSettings();
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();
//    settings.is_export_data_to_file = false;
//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


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
    NU_RK4 nu_cur = prepareNUASK();
//    nu_cur.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu_cur.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu_cur.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu_cur.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu_cur.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu_cur.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu_cur.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu_cur.P = ui->lineEdit_P->text().toDouble();
//    nu_cur.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu_cur.beta = fabs(nu_cur.P/(nu_cur.P_ud*g0));
//    nu_cur.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu_cur.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

//    nu_cur.v_KE = AGESK_to_KE(nu_cur.v_ASK);

//    nu_cur.Sbalxbezm = ui->lineEdit_Cx->text().toDouble()* ui->lineEdit_Sm->text().toDouble() / 2;


    Settings_RK4 settings = prepareSettings();
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();
//    settings.is_export_data_to_file = false;
//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


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



void MainWindow::on_B_start_Class_propagate_2_clicked()
{
    NU_RK4 nu = prepareNUASK();
//    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu.P = ui->lineEdit_P->text().toDouble();
//    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu.beta = fabs(nu.P/(nu.P_ud*g0));
//    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

//    nu.v_KE = AGESK_to_KE(nu.v_ASK);

    Settings_RK4 settings = prepareSettings();
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


    modeling_flight_2D raschet(nu, settings);
    Vector Uparam {nu.gamma_0, nu.d_gamma_dt};
    raschet.propUpr(Uparam);
    raschet.printCalcDataToFile(settings.file_name);

    QMessageBox::information(this, "Сообщение", "Расчет окончен");
}

void MainWindow::on_B_start_gradD_clicked() {
    NU_RK4 nu = prepareNUASK();
//    nu.v_ASK.x = ui->lineEdit_ASK_x0->text().toDouble();
//    nu.v_ASK.y = ui->lineEdit_ASK_y0->text().toDouble();
//    nu.v_ASK.z = ui->lineEdit_ASK_z0->text().toDouble();
//    nu.v_ASK.Vx = ui->lineEdit_ASK_Vx0->text().toDouble();
//    nu.v_ASK.Vy = ui->lineEdit_ASK_Vy0->text().toDouble();
//    nu.v_ASK.Vz = ui->lineEdit_ASK_Vz0->text().toDouble();
//    nu.v_ASK.m = ui->lineEdit_m->text().toDouble();
//    nu.P = ui->lineEdit_P->text().toDouble();
//    nu.P_ud = ui->lineEdit_P_ud->text().toDouble();
//    nu.beta = fabs(nu.P/(nu.P_ud*g0));
//    nu.gamma_0 = ui->lineEdit_gamma->text().toDouble()*M_PI/180;
//    nu.d_gamma_dt = ui->lineEdit_dgamma_dt->text().toDouble()*M_PI/180;

//    nu.v_KE = AGESK_to_KE(nu.v_ASK);

    Settings_RK4 settings = prepareSettings();
//    settings.file_name = ui->lineEdit_file_name->text() + ui->comboBox_file_type->currentText();

//    settings.dt = ui->lineEdit_dt->text().toDouble();
//    settings.hf = f_h(ui->lineEdit_a_per->text().toDouble());


    modeling_flight_2D raschet(nu, settings);
    Vector Uparam {nu.gamma_0, nu.d_gamma_dt};

    gradDescent opt(ui->TB_vivod);
    Vector Uisk = opt.gradD(raschet,Uparam);

    raschet.propUpr(Uisk);
    raschet.printCalcDataToFile(settings.file_name);

    QMessageBox::information(this, "Сообщение", "Расчет градиентного спуска окончен");
}


void MainWindow::on_B_start_DFP_clicked() {
    NU_RK4 nu = prepareNUASK();
    Settings_RK4 settings = prepareSettings();

    modeling_flight_2D raschet(nu, settings);
    Vector Uparam {nu.gamma_0, nu.d_gamma_dt};

    DFP opt(ui->TB_vivod);
    f_x f = [raschet](Vector U){
        modeling_flight_2D fl = modeling_flight_2D(raschet);
        return fl.propUpr(U);
    };
//    Vector Uisk = opt.calcDFP(f,Uparam);
    Vector Uisk = opt.calcDFP(raschet,Uparam);

    double t_res = raschet.propUpr(Uisk);
    raschet.printCalcDataToFile(settings.file_name);



    QMessageBox::information(this, "Сообщение", QString("Расчет градиентного спуска окончен, t = %1").arg(t_res));
}

//11772 13372 14172 17472
void MainWindow::on_B_start_Task_2v_clicked()
{
    NU_RK4 nu0 = prepareNUASK();
    Settings_RK4 settings0 = prepareSettings();
    Vector a_p,t_result;
    int i = 1;
    double ak = settings0.af; // Запоминаем величину большой полуоси конечного эллипса
    for (int a = nu0.v_KE.a+1; a < ak; a+=1000) {
//        if (a>14072) {
//            i++;
//        }
        NU_RK4 nu = nu0;      //Каждый новый рассчёт сброс параметров
        Settings_RK4 s = settings0;
        Vector Uparam {nu.gamma_0, nu.d_gamma_dt};

        s.hf = f_h(a);
        modeling_flight_2D raschet1(nu, s);

        DFP opt(ui->TB_vivod);
        opt.isToPrint = false;
        Vector Uisk1 = opt.calcDFP(raschet1,Uparam);
        double t1 = raschet1.propUpr(Uisk1);
        nu.upd_with_KE(raschet1.getCalcData().last().vKE.getKE_with_aem());

//        nu.v_ASK = KE_to_AGESK(newKE);
//        nu.v_KE = AGESK_to_KE(nu.v_ASK);
        s.hf = f_h(ak);
        modeling_flight_2D raschet2(nu,s);
        Vector Uisk2 = opt.calcDFP(raschet2,Uisk1);
        double t2 = raschet2.propUpr(Uisk2);

        a_p.append(a);
        t_result.append(t1+t2);

//        ui->TB_vivod->append(QString("Result for ap = %1:\n\t t = %2").
//                             arg(a).arg(t1+t2));
    }
    QStringList a_pSL = vec2strL(a_p,'f',0);
    QStringList t_rSL = vec2strL(t_result);
    for (int i = 0; i < t_result.length(); ++i) {
        if (t_result[i]>1000.0) {
            a_p.removeAt(i);
            t_result.removeAt(i);
        }
        ui->TB_vivod->append(QString("%1;%2").arg(a_pSL[i]).arg(t_rSL[i]));
    }



    plot_draw_graph(ui->Plot,"ap","t, с",a_p,t_result);
}

void MainWindow::on_B_start_Task_3v_clicked()
{
    NU_RK4 nu0 = prepareNUASK();
    Settings_RK4 settings0 = prepareSettings();
    Vector a_p1,a_p2,t_result;

    double da1 = 1000, da2 = 5000;
    double ak = settings0.af; // Запоминаем величину большой полуоси конечного эллипса
    for (int a1 = nu0.v_KE.a+1; a1 < ak-da1-da2-1; a1+=da1) {
        NU_RK4 nu = nu0;      //Каждый новый рассчёт сброс параметров
        Settings_RK4 s = settings0;
        Vector Uparam {nu.gamma_0, nu.d_gamma_dt};

        s.hf = f_h(a1);
        modeling_flight_2D raschet1(nu, s);

        DFP opt(ui->TB_vivod);
        opt.isToPrint = false;
        Vector Uisk1 = opt.calcDFP(raschet1,Uparam);
        double t1 = raschet1.propUpr(Uisk1);
        ui->TB_vivod->append(QString("Полет до а = %1:\t t1 = %2").arg(a1).arg(t1));
        for (int a2 = a1+da1; a2 < ak; a2+=da2) {
            cKE newKE1 = raschet1.getCalcData().last().vKE;
            newKE1.om = 0.0;
            newKE1.u  = 0.0;
            nu.v_ASK = KE_to_AGESK(newKE1);
            nu.v_KE = AGESK_to_KE(nu.v_ASK);
            s.hf = f_h(a2);
            modeling_flight_2D raschet2(nu,s);
            Vector Uisk2 = opt.calcDFP(raschet2,Uisk1);
            double t2 = raschet2.propUpr(Uisk2);
            ui->TB_vivod->append(QString("Полет до а = %1:\t t2 = %2").arg(a2).arg(t2));

            if (t2>2000.0) break;

            cKE newKE2 = raschet2.getCalcData().last().vKE;
            newKE2.om = 0.0;
            newKE2.u  = 0.0;
            nu.v_ASK = KE_to_AGESK(newKE2);
            nu.v_KE = AGESK_to_KE(nu.v_ASK);
            s.hf = f_h(ak);
            modeling_flight_2D raschet3(nu,s);
            Vector Uisk3 = opt.calcDFP(raschet3,Uisk2);
            double t3 = raschet3.propUpr(Uisk3);
            ui->TB_vivod->append(QString("Полет до а = %1:\t t3 = %2").arg(ak).arg(t3));

            a_p1.append(a1);
            a_p2.append(a2);
            t_result.append(t1+t2+t3);

            ui->TB_vivod->append(QString("Result for ap1 = %1, ap2 = %2:\t t = %3\t t1 = %4\t t2 = %5\t t3 = %6")
                                 .arg(a1).arg(a2).arg(t1+t2+t3).arg(t1).arg(t2).arg(t3));
        }
    }
    QStringList a_p1SL = vec2strL(a_p1,'f',0);
    QStringList a_p2SL = vec2strL(a_p2,'f',0);
    QStringList t_rSL = vec2strL(t_result);
    for (int i = 0; i < t_result.length(); ++i) {
        if (t_result[i]>900.0 || t_result[i]<500.0) {
            a_p1.removeAt(i);
            a_p2.removeAt(i);
            t_result.removeAt(i);
        }
        ui->TB_vivod->append(QString("%1 | %2 | %3").arg(a_p1SL[i]).arg(a_p2SL[i]).arg(t_rSL[i]));
    }

    plot_draw_graph(ui->Plot,"ap","t, с",a_p1,t_result);

}

void MainWindow::on_B_start_Task_4v_clicked()
{
    NU_RK4 nu0 = prepareNUASK();
    Settings_RK4 settings0 = prepareSettings();
    Vector a_p1,a_p2,a_p3,t_result;

    double da1 = 1000, da2 = 5000, da3 = 2000;
    double ak = settings0.af; // Запоминаем величину большой полуоси конечного эллипса
    for (int a1 = nu0.v_KE.a+1; a1 < ak-da1-da2-da3-1; a1+=da1) {
        NU_RK4 nu = nu0;      //Каждый новый рассчёт сброс параметров
        Settings_RK4 s = settings0;
        Vector Uparam {nu.gamma_0, nu.d_gamma_dt};

        s.hf = f_h(a1);
        modeling_flight_2D raschet1(nu, s);

        DFP opt;
        opt.isToPrint = false;
        Vector Uisk1 = opt.calcDFP(raschet1,Uparam);
        double t1 = raschet1.propUpr(Uisk1);
        ui->TB_vivod->append(QString("Полет до а = %1:\t t1 = %2").arg(a1).arg(t1));
        for (int a2 = a1+da1; a2 < ak-da2-da3; a2+=da2) {
            cKE newKE1 = raschet1.getCalcData().last().vKE;
            newKE1.om = 0.0;
            newKE1.u  = 0.0;
            nu.v_ASK = KE_to_AGESK(newKE1);
            nu.v_KE = AGESK_to_KE(nu.v_ASK);
            s.hf = f_h(a2);
            modeling_flight_2D raschet2(nu,s);
            Vector Uisk2 = opt.calcDFP(raschet2,Uisk1);
            double t2 = raschet2.propUpr(Uisk2);
            ui->TB_vivod->append(QString("Полет до а = %1:\t t2 = %2").arg(a2).arg(t2));

            for (int a3 = a2+da2; a3 < ak; a3+=da3) {
                cKE newKE2 = raschet2.getCalcData().last().vKE;
                newKE2.om = 0.0;
                newKE2.u  = 0.0;
                nu.v_ASK = KE_to_AGESK(newKE2);
                nu.v_KE = AGESK_to_KE(nu.v_ASK);
                s.hf = f_h(a3);
                modeling_flight_2D raschet3(nu,s);
                Vector Uisk3 = opt.calcDFP(raschet3,Uisk2);
                double t3 = raschet3.propUpr(Uisk3);
                ui->TB_vivod->append(QString("Полет до а = %1:\t t3 = %2").arg(a3).arg(t3));

                if (t2>2000.0) break;

                cKE newKE3 = raschet3.getCalcData().last().vKE;
                newKE3.om = 0.0;
                newKE3.u  = 0.0;
                nu.v_ASK = KE_to_AGESK(newKE3);
                nu.v_KE = AGESK_to_KE(nu.v_ASK);
                s.hf = f_h(ak);
                modeling_flight_2D raschet4(nu,s);
                Vector Uisk4 = opt.calcDFP(raschet4,Uisk3);
                double t4 = raschet4.propUpr(Uisk4);
                ui->TB_vivod->append(QString("Полет до а = %1:\t t4 = %2").arg(ak).arg(t4));

                a_p1.append(a1);
                a_p2.append(a2);
                a_p3.append(a3);
                t_result.append(t1+t2+t3+t4);

                ui->TB_vivod->append(QString("Result for ap1 = %1, ap2 = %2, ap3 = %3:\t t = %4\t t1 = %5\t t2 = %6\t t3 = %7\t t4 = %8")
                                     .arg(a1).arg(a2).arg(a3).arg(t1+t2+t3+t4).arg(t1).arg(t2).arg(t3).arg(t4));
            }
        }
    }
    QStringList a_p1SL = vec2strL(a_p1,'f',0);
    QStringList a_p2SL = vec2strL(a_p2,'f',0);
    QStringList a_p3SL = vec2strL(a_p3,'f',0);
    QStringList t_rSL  = vec2strL(t_result);
    for (int i = 0; i < t_result.length(); ++i) {
        if (t_result[i]>900.0 || t_result[i]<500.0) {
            a_p1.removeAt(i);
            a_p2.removeAt(i);
            a_p3.removeAt(i);
            t_result.removeAt(i);
        }
        ui->TB_vivod->append(QString("%1 | %2 | %3 | %4").arg(a_p1SL[i]).arg(a_p2SL[i]).arg(t_rSL[i]));
    }

    plot_draw_graph(ui->Plot,"ap","t, с",a_p1,t_result);
}

void MainWindow::plot_draw_graph(QCustomPlot *curPlot, QString xName, QString yName, Vector x, Vector y)
{

    curPlot->setLocale(QLocale(QLocale::Russian, QLocale::Russia));
//    ui->Plot->graph(0)->setAntialiased(false);


    int r = rand() % 255+1,
        b = rand() % 255+1,
        g = rand() % 255+1;

//    ui->Plot->addGraph();

//    QCPCurve a(ui->Plot->xAxis,ui->Plot->yAxis);
//    a.setData(x,x,y);

    curPlot->graph(0)->data().data()->clear(); // ????????
    curPlot->graph(0)->setPen(QPen(QColor(r,g,b,255)));
    curPlot->graph(0)->setData(x, y);
    ui->Plot->yAxis->setLabel(yName);
    ui->Plot->xAxis->setLabel(xName);
    curPlot->rescaleAxes(true);
    curPlot->replot();

}

void MainWindow::plot_draw_curves(QCustomPlot *curPlot, QString xName, QString yName, Vector x1, Vector x2, Vector y)
{
    curPlot->setLocale(QLocale(QLocale::Russian, QLocale::Russia));
    for (QCPCurve* curve : curves) {
//        ui->Plot->re
    }


//    ui->Plot->graph(0)->setAntialiased(false);


    int r = rand() % 255+1,
        b = rand() % 255+1,
        g = rand() % 255+1;

//    ui->Plot->addGraph();

//    QCPCurve a(ui->Plot->xAxis,ui->Plot->yAxis);
//    a.setData(x,x,y);

    curPlot->graph(0)->data().data()->clear(); // ????????
    curPlot->graph(0)->setPen(QPen(QColor(r,g,b,255)));
    curPlot->graph(0)->setData(x1, y);
    ui->Plot->yAxis->setLabel(yName);
    ui->Plot->xAxis->setLabel(xName);
    curPlot->rescaleAxes(true);
    curPlot->replot();
}






void MainWindow::on_B_start_Task_2vmin_clicked()
{
    NU_RK4 nu = prepareNUASK();
    Settings_RK4 settings = prepareSettings();

    modeling_flight_2D raschet(nu, settings);
    Vector Uparam {nu.gamma_0, nu.d_gamma_dt};


    OptimaFlight2D task(nu,settings,{},{},{},ui->TB_vivod);

    gradDescent opt({1e-2},1e-6,{10},ui->TB_vivod);

    f_x f = [&task,nu,settings](Vector U){
        return task.task2v(nu,settings,{U[0],settings.af});
    };

    Vector Uisk = opt.gradD(f,{12000});

    double t_res = task.task2v(nu,settings,{Uisk[0],settings.af});


    QMessageBox::information(this, "Сообщение", QString("Расчет градиентного спуска окончен, t_res = %1").arg(t_res));
}


void MainWindow::on_B_start_Task_Imp_clicked()
{
    NU_RK4 nu = prepareNUASK();
    Settings_RK4 s = prepareSettings();

    gradDescent opt({1e-2},1e-6,{10},ui->TB_vivod);

    QTextBrowser *TB = ui->TB_vivod;

    f_x f = [nu,s,TB](Vector U){
        ImpTask task(nu.v_KE,
                     cKE(s.af,0.0,0.0,0.0,0.0,0.0),
                     nu.beta,
                     nu.Wist,
                     nu.v_ASK.m,
                     mu,
                     TB);
        return task.findFulldVxap2V(U.first());
    };

    QStringList h1 = {"ap","dV1","mk1","dt1","dV2","mk2","dt2","dV","dt"};
    QStringList h2 = {"ap1","ap2","dV1","mk1","dt1","dV2","mk2","dt2","dV3","mk3","dt3","dV","dt"};
    QString     sep = ";",
                out = strL2str(h1,sep)+"\n";

    for (double a = nu.v_KE.a+10; a < s.af; a+=50) {
        ImpTask task(nu.v_KE,
                     cKE(s.af,0.0,0.0,0.0,0.0,0.0),
                     nu.beta,
                     nu.Wist,
                     nu.v_ASK.m,
                     mu,
                     TB);
        out += strL2str(vec2strL(task.findFullVals({a})),sep) + "\n";
    }
    out += "\n\n\n"+strL2str(h2,sep)+"\n";
    for (double a1 = nu.v_KE.a+10; a1 < s.af; a1+=500) {
        for (double a2 = a1+10; a2 < s.af; a2+=1000) {
            ImpTask task(nu.v_KE,
                         cKE(s.af,0.0,0.0,0.0,0.0,0.0),
                         nu.beta,
                         nu.Wist,
                         nu.v_ASK.m,
                         mu,
                         TB);
            out += strL2str(vec2strL(task.findFullVals({a1,a2})),sep) + "\n";
        }
    }

    TB->append(out);
//    Vector Uisk = opt.gradD(f,{12000});
//    double dV = task.findFulldVxap2V(Uisk.first());
    QMessageBox::information(this, "Сообщение",
                             QString("Расчет градиентного спуска окончен, dV = %1"));

}


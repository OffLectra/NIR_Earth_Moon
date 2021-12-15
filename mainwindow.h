#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <cmath>
#include "perevod_.h"
#include "RK4_integrator.h"
#include "modeling_flight_rk4.h"
#include "optima.h"
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    NU_RK4 prepareNUASK ();
    Settings_RK4 prepareSettings ();

private slots:

    void on_pushButton_KE_to_AGESK_clicked();

    void on_pushButton_AGESK_to_KE_clicked();

    void on_B_start_clicked();

    void on_B_Startatm_clicked();

    void on_B_Optima_clicked();

    void on_B_Opt_DFP_clicked();

    void on_B_Optima_new_clicked();

    void on_B_start_Class_propagate_clicked();

    void on_B_start_Class_propagate_2_clicked();

    void on_B_start_gradD_clicked();

    void on_B_start_DFP_clicked();

    void on_B_start_Task_2v_clicked();

    void on_B_start_Task_3v_clicked();

    void on_B_start_Task_4v_clicked();

private:
    void plot_draw_graph(QCustomPlot* curPlot,QString xName, QString yName, Vector x, Vector y);
    void plot_draw_curves(QCustomPlot* curPlot,QString xName, QString yName, Vector x1, Vector x2, Vector y);
    QVector<QCPCurve*> curves;

    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H

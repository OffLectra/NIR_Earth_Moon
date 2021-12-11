#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <cmath>
#include "perevod_.h"
#include "RK4_integrator.h"
#include "modeling_flight_rk4.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();



private slots:

    void on_pushButton_KE_to_AGESK_clicked();

    void on_pushButton_AGESK_to_KE_clicked();

    void on_B_start_clicked();

    void on_B_Startatm_clicked();

    void on_B_Optima_clicked();

    void on_B_Opt_DFP_clicked();

    void on_B_Optima_new_clicked();

    void on_B_start_Class_propagate_clicked();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H

#include "mainwindow.h"

#include <QApplication>
#include <QLabel> //библиотека, которая отвечает за надписи на форме
//здесь начинает работу прога


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show(); //вызов формы



    return a.exec();
}

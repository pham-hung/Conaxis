#ifndef TERZAGHI1D_H
#define TERZAGHI1D_H

#include <QDialog>
#include <Eigen/Dense>
#include "math.h"
#include "writetofile.h"
#include <QString>
#include <iostream>
#include <QFileDialog>
#include <QInputDialog>

using namespace std;
using namespace Eigen;

namespace Ui {
class Terzaghi1D;
}

class Terzaghi1D : public QDialog
{
    Q_OBJECT

public:
    explicit Terzaghi1D(QWidget *parent = 0);
    ~Terzaghi1D();
    void getUserData();
    void calculateSolution();

private slots:
    void on_closeButton_clicked();
    void on_okButton_clicked();

private:
    Ui::Terzaghi1D *ui;
    double k,K,G,v,dt,ns,totalTime,Cv,h,z,p0;
    MatrixXd result;
};

#endif // TERZAGHI1D_H

#ifndef MATERIALBASE_H
#define MATERIALBASE_H

#include <Eigen/Dense>
#include <iostream>
#include <QString>

using namespace std;
using namespace Eigen;
class MaterialBase
{

public:
    MaterialBase();
    int matIndex=1;
    int KFunction=0;
    MatrixXd KCurve=MatrixXd::Zero(1,1);
    QString fileNameK;
    double poission=0.2;
    int kFunction=0;
    MatrixXd kCurve=MatrixXd::Zero(1,1);
    QString fileNamek;
    double ratio=1;
    double gf=10;
    double e0=1.8;
    double Cd=1.0f;
private:

};

#endif // MATERIALBASE_H

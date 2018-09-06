#ifndef BOUNDARYCONDITIONBASE_H
#define BOUNDARYCONDITIONBASE_H
#include <iostream>
#include <QString>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class BoundaryConditionBase
{
public:
    BoundaryConditionBase();
    int boundaryIndex;
    QString boundaryName;
    int boundaryType;
    int loadType;
    double val0;
    double val1;
    QString fileName;
    MatrixXd loadCurve=MatrixXd::Zero(1,1);
};

#endif // BOUNDARYCONDITIONBASE_H

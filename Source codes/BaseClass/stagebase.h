#ifndef STAGEBASE_H
#define STAGEBASE_H
#include <iostream>
#include <QString>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


class StageBase
{
public:
    StageBase();
    int stageIndex;
    QString stageName;
    int stageType;
    double dt;
    double t0;
    double t1;
    int subStep;
    int gravityLoad;
    int timeStepType;
    QString fileName;
    MatrixXd timeStep;
};

#endif // STAGEBASE_H

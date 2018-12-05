#ifndef STAGEBOUNDARYBASE_H
#define STAGEBOUNDARYBASE_H
#include <iostream>
#include <QString>
#include "boundaryconditionbase.h"
#include "stagebase.h"
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
class StageBoundaryBase:private BoundaryConditionBase,private StageBase
{
public:
    StageBoundaryBase();
    int stageIndex;
    int stageName;
    int stageBoundaryIndex;

    vector<int> v_stageBoundaryIndex;
    vector<MatrixXd> v_nodeList;
    vector<QString> v_nodeFileName;
    vector<int> v_assignType;
    vector<int> v_boundaryIndex;

    vector<double> v_x0;
    vector<double> v_x1;
    vector<double> v_y0;
    vector<double> v_y1;
    vector<bool> v_gradientBool;
    vector<double> v_aFactor;
    vector<double> v_bFactor;
    vector<double> v_cFactor;
};

#endif // STAGEBOUNDARYBASE_H

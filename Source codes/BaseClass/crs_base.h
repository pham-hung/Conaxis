#ifndef CRS_BASE_H
#define CRS_BASE_H

#include <QObject>
#include <Eigen/Dense>

using namespace Eigen;

class CRS_Base
{
public:
    CRS_Base();
    double H0=0.0254;
    double R=0.03175;
    double e0=1.84;
    double possionRatio=0.4;
    MatrixXd crsData;
    int analysisType=0;
    bool crsFlag=false;
};

#endif // CRS_BASE_H

#ifndef GEOMETRYBASE_H
#define GEOMETRYBASE_H

#include <QObject>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class GeometryBase
{
public:
    GeometryBase();
    double re=0.565;
    double rw=0.05;
    double rs=0.15;
    double length=1.0;
    int numberOfLayer=1;
    double qw=1e-5;
    int numberOfElementSmear=1;
    int numberOfElementSoil=3;
    double surfaceElevation=1;
    int analysisType=1;
    MatrixXd layerInfo=MatrixXd::Zero(0,0);
    int defaultSubLayer=5;
};

#endif // GEOMETRYBASE_H

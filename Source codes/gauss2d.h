#ifndef GAUSS2D_H
#define GAUSS2D_H
#include <Eigen/Dense>
#include <iostream>
#include "math.h"

using namespace std;
using namespace Eigen;
class Gauss2D
{
public:
    Gauss2D();
    MatrixXd gauss3,gauss4;
    void createGauss3(); //3 gauss-Point
    void createGauss4(); //4 gauss-Point
};

#endif // GAUSS2D_H

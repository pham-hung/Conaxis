#include "gauss2d.h"

Gauss2D::Gauss2D()
{
    createGauss3();
    createGauss4();
}

void Gauss2D::createGauss3()
{
    gauss3=MatrixXd::Zero(3,3);
    gauss3<<1.0f/6.0f,1.0f/6.0f,1.0f/6.0f,
            2.0f/3.0f,1.0f/6.0f,1.0f/6.0f,
            1.0f/6.0f,2.0f/3.0f,1.0f/6.0f;
}

void Gauss2D::createGauss4()
{
    //2 x 2 gauss-points
    gauss4=MatrixXd::Zero(4,3);
    double tempVal=1.0f/sqrt(3.0f);
    gauss4<<tempVal,tempVal,1.0f,
            tempVal,-tempVal,1.0f,
            -tempVal,tempVal,1.0f,
            -tempVal,-tempVal,1.0f;
}

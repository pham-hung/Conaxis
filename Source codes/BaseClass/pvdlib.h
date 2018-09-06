#ifndef PVDLIB_H
#define PVDLIB_H
#include <math.h>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>
#include <QObject>
#include <QDebug>

using namespace std;
using namespace boost::math;

class PvdLib : public QObject
{
    Q_OBJECT
public:
    explicit PvdLib(QObject *parent = nullptr);
    //math function
    double factorial (int n);
    double besselj(double u,double x);
    double bessely(double u,double x);

    //PVD no smear function
    double PvdNoSmearBessel(double x,double n);
    double rootFunctionNoSmearBiSection(double initial,double sizeStep,double n);
    void findRootNoSmear(vector<double> &alpha,double n);
    double analyticalNoSmear(double r,double time);
    void setParameterNoSmear(double re,double rw,double Cvr,double p0);
    void printInforNoSmear();

    //PVD smear zone function
    double PvdSmearBessel(double x,double kh,double ks,double s,double n);
    double rootFunctionSmearBiSection(double initial,double sizeStep,double kh, double ks, double s, double n);
    void findRootSmear(vector<double> &alpha, double kh, double ks, double s, double n);
    double analyticalSmear(double r,double time);
    void setParameterSmear(double re, double rw, double rs, double Cvr,double kh, double ks,double p0);
    void printInforSmear();

    //support function
    void pauseSystem();

signals:

public slots:

private:
    //No Smear
    double re=0.565;
    double rw=0.0264;
    double n=re/rw; //Not depend on r
    double dt=86400;
    double r=0.565;
    double p0=100;
    double Cvr=1e-7;

    double kh=1e-9;
    double ks=0.5*kh;
    double rs=0.102;
    double s=rs/rw;

    vector<double> alpha;
    bool rootNoSmearCheck=false;
    bool rootSmearCheck=false;
};

#endif // PVDLIB_H

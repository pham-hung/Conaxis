#include "pvdlib.h"

PvdLib::PvdLib(QObject *parent) : QObject(parent)
{

}
//-------------------------------------------------
//Math functions
double PvdLib::factorial(int n)
{
    double f=1;
    if (n==0){f=1;}
    for (int j=1;j<=n;j++)
    {
        f=f*j;
    }
    return f;
}

double PvdLib::besselj(double u, double x)
{
    return cyl_bessel_j(u,x);
}

double PvdLib::bessely(double u, double x)
{
    return cyl_neumann(u,x);
}

//-------------------------------------------------
//No Smear Zone
double PvdLib::PvdNoSmearBessel(double x, double n)
{
    double fx;
    fx=besselj(1,x*n)*bessely(0,x)-bessely(1,x*n)*besselj(0,x);
    return fx;
}

double PvdLib::rootFunctionNoSmearBiSection(double initial, double sizeStep, double n)
{
    double x0,x1,x;
    double fx0,fx1,fx;
    double tol=1e-6;
    x0=initial;
    x1=x0+sizeStep;
    fx0=PvdNoSmearBessel(x0,n);
    fx1=PvdNoSmearBessel(x1,n);

    if(fx0>0)
    {
        int jj=0;
        while(fx1>0)
        {
            x1=x0+jj*sizeStep;
            fx1=PvdNoSmearBessel(x1,n);
            jj++;
        }
    }

    if(fx0<0)
    {
        int jj=0;
        while(fx1<0)
        {
            x1=x0+jj*sizeStep;
            fx1=PvdNoSmearBessel(x1,n);
            jj++;
        }
    }
    x=0.5*(x0+x1);
    fx=PvdNoSmearBessel(x,n);

    while(abs(fx)>tol)
    {
        double xOld=x;

        if((fx0*fx)>0)
        {
            x0=x;
        }
        else if ((fx1*fx)>0)
        {
            x1=x;
        }

        fx0=PvdNoSmearBessel(x0,n);
        fx1=PvdNoSmearBessel(x1,n);
        x=0.5*(x0+x1);
        fx=PvdNoSmearBessel(x,n);

        if(abs(x-xOld)<tol)
        {

            break;
            return x;
        }
    }
    return x;
}

void PvdLib::findRootNoSmear(vector<double> &alpha, double n)
{
    double initial=rootFunctionNoSmearBiSection(0.001,0.001,n);
    int j=0;
    alpha.clear();
    alpha.resize(10);
    alpha[j]=initial;
    while(j<(alpha.size()-1))
    {
        double val=rootFunctionNoSmearBiSection(alpha[j]+initial,initial,n);
        if(abs(val-alpha[j])>initial)
        {
            j=j+1;
            alpha[j]=val;
        }
    }
    rootNoSmearCheck=true;
    rootSmearCheck=false;
}

double PvdLib::analyticalNoSmear(double r,double time)
{
    if(rootNoSmearCheck==false)
    {
        findRootNoSmear(alpha,n);
    }
    if(time==0)
    {
        return p0;
    }
    else
    {
        if(r<=rw){return 0;}
        else{
            double de=2*re;
            double Tr=Cvr*time/de/de;
            double pore=0;
            for (int i=0;i<alpha.size();i++)
            {
                double alpha0=alpha[i];
                double u1a=besselj(1,alpha0)*bessely(0,alpha0)-bessely(1,alpha0)*besselj(0,alpha0);
                double u0an=besselj(0,alpha0*n)*bessely(0,alpha0)-bessely(0,alpha0*n)*besselj(0,alpha0);
                double u0ar=besselj(0,alpha0*r/rw)*bessely(0,alpha0)-bessely(0,alpha0*r/rw)*besselj(0,alpha0);
                double ea=exp(-4*alpha0*alpha0*n*n*Tr);
                double ualpha=-2*u1a*u0ar*ea/alpha0/(n*n*u0an*u0an-u1a*u1a);
                pore=pore+ualpha;
            }
            return pore*p0;
        }
    }
}

void PvdLib::setParameterNoSmear(double re, double rw, double Cvr,double p0)
{
    this->re=re;
    this->rw=rw;
    this->Cvr=Cvr;
    this->p0=p0;
    this->n=re/rw;
    rootNoSmearCheck=false;
    rootSmearCheck=false;
    printInforNoSmear();
}

void PvdLib::printInforNoSmear()
{
    cout<<"------------------------------------"<<endl;
    cout<<"INPUT PARAMTERS- NO SMEAR ZONE"<<endl;
    cout<<"Radius of unit cell           : "<<re<<endl;
    cout<<"Radius of PVDs                : "<<rw<<endl;
    cout<<"Horizontal Consolidation Coeff: "<<Cvr<<endl;
    cout<<"Initial pore pressure         : "<<p0<<endl;
    cout<<"------------------------------------"<<endl;
}

//-------------------------------------------------
//With Smear Zone
double PvdLib::PvdSmearBessel(double x, double kh, double ks, double s, double n)
{
    double fx;
    double U0as=besselj(0,x*s)*bessely(1,x*n)-besselj(1,x*n)*bessely(0,x*s);
    double U1as=besselj(1,x*s)*bessely(1,x*n)-besselj(1,x*n)*bessely(1,x*s);
    fx=ks*U0as/(kh*x*s*log(s))+U1as;
    return fx;
}

double PvdLib::rootFunctionSmearBiSection(double initial, double sizeStep, double kh, double ks, double s, double n)
{
    double x0,x1,x;
    double fx0,fx1,fx;
    double tol=1e-6;
    x0=initial;
    x1=x0+sizeStep;
    fx0=PvdSmearBessel(x0,kh,ks,s,n);
    fx1=PvdSmearBessel(x1,kh,ks,s,n);

    if(fx0>0)
    {
        int jj=0;
        while(fx1>0)
        {
            x1=x0+jj*sizeStep;
            fx1=PvdSmearBessel(x1,kh,ks,s,n);
            jj++;
        }
    }

    if(fx0<0)
    {
        int jj=0;
        while(fx1<0)
        {
            x1=x0+jj*sizeStep;
            fx1=PvdSmearBessel(x1,kh,ks,s,n);
            jj++;
        }
    }
    x=0.5*(x0+x1);
    fx=PvdSmearBessel(x,kh,ks,s,n);

    while(abs(fx)>tol)
    {
        double xOld=x;

        if((fx0*fx)>0)
        {
            x0=x;
        }
        else if ((fx1*fx)>0)
        {
            x1=x;
        }

        fx0=PvdSmearBessel(x0,kh,ks,s,n);
        fx1=PvdSmearBessel(x1,kh,ks,s,n);
        x=0.5*(x0+x1);
        fx=PvdSmearBessel(x,kh,ks,s,n);

        if(abs(x-xOld)<tol)
        {

            break;
            return x;
        }
    }
    return x;
}

void PvdLib::findRootSmear(vector<double> &alpha, double kh, double ks, double s, double n)
{
    double initial=rootFunctionSmearBiSection(0.001,0.001,kh,ks,s,n);
    int j=0;
    alpha.resize(10);
    alpha[j]=initial;
    while(j<(alpha.size()-1))
    {
        double val=rootFunctionSmearBiSection(alpha[j]+initial,initial,kh,ks,s,n);
        if(abs(val-alpha[j])>initial)
        {
            j=j+1;
            alpha[j]=val;
        }
    }
    rootSmearCheck=true;
    rootNoSmearCheck=false;
}

double PvdLib::analyticalSmear(double r,double time)
{
    if(rootSmearCheck==false)
    {
        findRootSmear(alpha,kh,ks,s,n);
    }
    if(time==0){return p0;}
    else
    {
        if(r<=rw){return 0;}
        else
        {
            double de=2*re;
            double Tr=Cvr*time/de/de;
            double pi=3.1415926535897932384626433832795029;
            double pore=0;
            double a,U0as,U1as,U0ar,m,ut,um,u;
            for (int i=0;i<alpha.size();i++)
            {
                a=alpha[i];
                U0as=besselj(0,a*s)*bessely(1,a*n)-besselj(1,a*n)*bessely(0,a*s);
                U1as=besselj(1,a*s)*bessely(1,a*n)-besselj(1,a*n)*bessely(1,a*s);
                U0ar=besselj(0,a*r/rw)*bessely(1,a*n)-besselj(1,a*n)*bessely(0,a*r/rw);
                m=-4*n*n*a*a*Tr;
                ut=-2*U1as*U0ar*exp(m)/s/a;
                um=4/(pi*pi*a*a*s*s)-U0as*U0as-U1as*U1as;
                u=ut/um;
                pore=pore+u;
            }
            return pore*p0;
        }
    }
}

void PvdLib::setParameterSmear(double re, double rw, double rs, double Cvr, double kh, double ks, double p0)
{
    this->re=re;
    this->rw=rw;
    this->rs=rs;
    this->Cvr=Cvr;
    this->kh=kh;
    this->ks=ks;
    this->p0=p0;
    n=re/rw;
    s=rs/rw;
    rootNoSmearCheck=false;
    rootSmearCheck=false;
    printInforSmear();
}

void PvdLib::printInforSmear()
{
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"INPUT PARAMTERS- WITH SMEAR ZONE"<<endl;
    cout<<"Radius of unit cell                        : "<<re<<endl;
    cout<<"Radius of PVDs                             : "<<rw<<endl;
    cout<<"Radius of Smear Zone                       : "<<rs<<endl;
    cout<<"Horizontal Consolidation Coeff             : "<<Cvr<<endl;
    cout<<"Soil horizontal hydraulic conductivity     : "<<kh<<endl;
    cout<<"SmearZone horizontal hydraulic conductivity: "<<ks<<endl;
    cout<<"Initial pore pressure                      : "<<p0<<endl;
    cout<<"--------------------------------------------------------"<<endl;
}

void PvdLib::pauseSystem()
{
    do {
        cout << '\n' << "Press the Enter key to continue.";
    } while (cin.get() != '\n');

}

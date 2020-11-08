#include "bslipne.h"
#include <cmath>

#include <iostream>


namespace
{

    class S
    {
    public:

        static double s( double x)
        {
            if (x>-2 && x<=-1){return pow(2+x,3)/4.;}
            else if (x>-1 && x<=0){return (pow(2+x,3)-4*pow(1+x,3))/4.;}
            else if (x>0 && x<=1){return (pow(2-x,3)-4*pow(1-x,3))/4.;}
            else if (x>1 && x<=2){return pow(2-x,3)/4.;}
            else return 0;
        }

        static double ds( double x)
        {
            if (x>-2 && x<=-1) return 3*pow(2+x,2)/4.;
            else if (x>-1 && x<=0) return 3*(pow(2+x,2)-4*pow(1+x,2))/4.;
            else if (x>0 && x<=1) return -3*(pow(2-x,2)-4*pow(1-x,2))/4.;
            else if (x>1 && x<=2) return -3*pow(2-x,2)/4.;
            else return 0;
        }
    };

}

Phi::Phi(double h)
{
    this->h = h;
}

double Phi::operator()(double)
{
    std::cout << "phi ";
    return 0;
}

double Phi::DPhi(double)
{
    return 0;
}

Phi_0::Phi_0(double h)
{
    this->h = h;
}

double Phi_0::operator()(double x)
{
    return S::s(x/h)-4*S::s((x+h)/h);
}

double Phi_0::DPhi(double x)
{
    return (1./h)*(S::ds(x/h)-4*S::ds((x+h)/h));
}



Phi_1::Phi_1(double h)
{
    this->h = h;
}

double Phi_1::operator()(double x)
{
    return S::s((x-h)/h)-S::s((x+h)/h);
}

double Phi_1::DPhi(double x)
{
    return (1./h)*(S::ds((x-h)/h)-S::ds((x+h)/h));
}

Phi_i::Phi_i(double h, int i)
{
    this->h = h;
    this->i = i;
}

double Phi_i::operator()(double x)
{
    return S::s((x-i*h)/h);
}

double Phi_i::DPhi(double x)
{
    return (1./h)*(S::ds((x-i*h)/h));
}

Phi_n::Phi_n(double h, int n)
{
    this->h = h;
    this->n = n;
}

double Phi_n::operator()(double x)
{
    return S::s((x-n*h)/h)-S::s((x-(n+2)*h)/h);
}

double Phi_n::DPhi(double x)
{
    return (1./h)*(S::ds((x-n*h)/h)-S::ds((x-(n+2)*h)/h));
}

Phi_n_1::Phi_n_1(double h, int n)
{
    this->h = h;
    this->n = n;
}

double Phi_n_1::operator()(double x)
{
    return S::s((x-(n+1)*h)/h)-4*S::s((x-(n+2)*h)/h);
}

double Phi_n_1::DPhi(double x )
{
    return (1./h)*(S::ds((x-(n+1)*h)/h)-4*S::ds((x-(n+2)*h)/h));
}





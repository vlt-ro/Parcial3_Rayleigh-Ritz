#include "csrayleighritz.h"
#include <cstdlib>
#include <cmath>
#include "../../utils/integracion.h"
#include "../../utils/sistemaslineales.h"

using namespace integracion;
using namespace sistemas_lineales;
using namespace std;

csRayleighRitz::csRayleighRitz(std::size_t n)
{
    this->n = n;
}
csRayleighRitz::~csRayleighRitz()
{
}

vector<double> csRayleighRitz::solve(double(*p)(double),
                     double(*q)(double),
                     double(*f)(double))
{
    vector<vector<double>> A(n+2,vector<double>(n+2)); // Matriz A de dimension n+2 x n+2
    vector<double> b(n+2); // Vector b de dimension n+2
    vector<double> c(n+2); // Vector c de dimension n+2


    /*******************************************/
    /*   Step 1: Cálculo del coeficientes h    */
    /*******************************************/
    h = 1/static_cast<double>(1+n);

    // Los pasos 2 y 3 se representan con funciones

    /**************************************************/
    /* Step 4: Definición de los elementos de la base */
    /**************************************************/
    for(unsigned i=0; i<= n+1; ++i)
        phi.push_back( csBasis(i,n,h));

    /***************************************/
    /* Step 5-8: Calculo de las componentes*/
    /* de la matriz A                      */
    /***************************************/
    for (int i = 0; i < n+2; ++i)
    {
      for (int j = i; j <= std::min(i+3,int(n+1)); ++j)
      {
        double L = std::max(x_i(j-2), 0.);
        double U = std::min(x_i(i+2), 1.);

        A[i][j] =  simpson_rule( L, U, 20, [i,j,p,q,this](float x)->float {return p(x)*this->phi[i].dPhi(x)*this->phi[j].dPhi(x)+q(x)*this->phi[i](x)*this->phi[j](x);});

        if (i!=j) A[j][i] = A[i][j]; // La matriz es simetrica
      }

      if (i>=4) { for (int j = 0; j <= i-4; ++j) A[i][j] = 0;}
      if (i<=n-3) { for (int j = i+4; j <= n+1; ++j) A[i][j] = 0;}
    }

    /*****************************************/
    /* Step 9: Calculo de las componentes del*/
    /* vector b                              */
    /*****************************************/
    for (int i = 0; i < n+2; ++i)
    {
      double L = std::max(x_i(i-2),(double)(0));
      double U = std::min(x_i(i+2),(double)(1));

      b[i] = simpson_rule(L,U,20,[i,f,this](float x)->float {return f(x)*this->phi[i](x);});
    }

    /****************************************/
    /* Step 10: Se encuentra c, resolviendo */
    /* el sistema de ecuaciones             */
    /****************************************/
    gauss_jordan(n+2, A, b, c);

    this->c = c;
    return c;
}

vector<csBasis> &csRayleighRitz::getBasis()
{
    return phi;
}

double csRayleighRitz::eval(double x)
{
    double rs = 0;

    for(size_t i=0; i<c.size(); ++i)
        rs += c[i] * phi[i](x);

    return rs;
}

double csRayleighRitz::x_i(int i)
{
    if( i < 0) return 0;
    else if(i < n+2) return h*i;
    else return 1;
}




csBasis::csBasis(int i,int n, double h): phi(nullptr)
{
    setMembers(i,n,h);
}

csBasis::csBasis(const csBasis &old): phi(nullptr)
{
    setMembers(old.i,old.n,old.h);
}

csBasis::~csBasis()
{
    // Liberar memória
    if(phi)
        delete phi;
}

double csBasis::operator()(double x)
{
    return (*phi)(x);
}

double csBasis::dPhi(double x)
{
    return phi->DPhi(x);
}


void csBasis::setMembers(int i, int n, double h)
{
    this->i = i;
    this->n = n;
    this->h = h;

    /* Seleccionar un elemento de la base específico */
    if(i==0)
        phi = new Phi_0(h);
    else if(i==1)
        phi = new Phi_1(h);
    else if(i<n)
        phi = new Phi_i(h,i);
    else if(i==n)
        phi = new Phi_n(h,n);
    else if(i==(n+1))
        phi = new Phi_n_1(h,n);
    else
        throw "Error: Parametros invalidos";
}

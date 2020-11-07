#include <stdlib.h>
#include<iomanip>
#include <iostream>
#include <vector>
#include <cmath>

#include "../utils/Cholesky/cholesky.h"
#include "../utils/Trapezoide/trapezoide.h" 
#include "../utils/BSpline/b-spline.h" 
#include "problem.h"
#include "rayleigh-ritz.h"

using namespace std;

Rayleigh_Ritz::Rayleigh_Ritz(int n)
{
  setParameters(n);
}

void Rayleigh_Ritz::setParameters(int n){ this->n = n; }

float Rayleigh_Ritz::x_j(int j)
{
  /***********************************/
  /* Step 2: Se establecen los nodos */
  /***********************************/
  if (j==-2 || j==-1) return 0;
  else if (j==n+2 || j==n+3) return 1;
  else return j*h;  
}

float Rayleigh_Ritz::y(float x, const vector<float> c)
{
  float y_value = 0;

  for (int i = 0; i < n+2; ++i)
  {
    Basis_j.setParameters(i, n);
    y_value += c[i]*Basis_j.phi(x);
  }

  return y_value; 
}


float Rayleigh_Ritz::eval_aij(float x){ return myProblem.p(x)*Basis_i.dphi(x)*Basis_j.dphi(x)+myProblem.q(x)*Basis_i.phi(x)*Basis_j.phi(x); }

float Rayleigh_Ritz::Aij(int i, int j, float lim_inf, float lim_sup)
{
  Basis_i.setParameters(i, n);
  Basis_j.setParameters(j, n);
  
  return Integral.solve(lim_inf,lim_sup,[this](float x)->float {return this->eval_aij(x);});
}

float Rayleigh_Ritz::eval_bi(float x) { return myProblem.f(x)*Basis_i.phi(x);}

float Rayleigh_Ritz::bi(int i, float lim_inf, float lim_sup)
{
  Basis_i.setParameters(i, n);

  return Integral.solve(lim_inf,lim_sup,[this](float x)->float {return this->eval_bi(x);});
}


void Rayleigh_Ritz::solve()
{
  vector<vector<float>> A(n+2,vector<float>(n+2)); // Matriz A de dimension n+2 x n+2
  vector<float> b(n+2); // Vector b de dimension n+2
  vector<float> c(n+2); // Vector c de dimension n+2

  /*************************************/
  /* Step 1: Calculo del coeficiente h */
  /*************************************/
  h = 1./(1+n);


  /*****************************************/
  /* Step 6,7,8: Calculo de las componentes*/
  /* de la matriz A                        */
  /*****************************************/
  for (int i = 0; i < n+2; ++i) 
  {
    for (int j = i; j <= min(i+3,n+1); ++j) 
    {
      L = max(x_j(j-2), (float)(0));
      U = min(x_j(i+2), (float)(1));

      A[i][j] = Aij(i,j, L, U);

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
    L = max(x_j(i-2),(float)(0));
    U = min(x_j(i+2),(float)(1));
  
    b[i] = bi(i,L,U);
  }


  /****************************************/
  /* Step 10: Se encuentra c, resolviendo */
  /* el sistema de ecuaciones             */
  /****************************************/
  Cholesky SisEcuaciones(n+2, A, b, c);
  SisEcuaciones.solve();

  /*************************************/
  /* Se imprimen los resultados        */
  /*************************************/
  cout<<"|  i|    c_i      |    x_i      |   y(x_i)    |y_exacta(x_i)||y(x_i)-y_exacta(x_i)||"<<endl;

  for (int i = 0; i < n+2; ++i)
  {
    //cout<<y_i<<" ,";
    float y_i = y(x_j(i), c);
    float y_exacta = sin(M_PI*x_j(i));

    cout<<"|"<<setw(3)<<i;
    cout<<"|"<<setw(13)<<c[i];
    cout<<"|"<<setw(13)<<x_j(i);
    cout<<"|"<<setw(13)<<y_i;
    cout<<"|"<<setw(13)<<y_exacta;
    cout<<"|"<<setw(21)<<abs(y_exacta-y_i)<<" |\n";
  }  

};
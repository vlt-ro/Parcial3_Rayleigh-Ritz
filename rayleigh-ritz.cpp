/**
 * @author Valentina Roquemen Echeverry
 * @brief  Solucion de sistemas de ecuaciones de la
 *         forma Ac=b usando el metodo de Cholesky
 *         A es una matriz definida positiva
 */

#include <stdlib.h>
#include <time.h>
#include<iomanip>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>    

#include "utils/Cholesky/cholesky.h"
#include "utils/Trapezoide/trapezoide.h" 
#include "utils/BSpline/b-spline.h" 
#include "problem.h"
#include "rayleigh-ritz.h"

using namespace std;

Rayleigh_Ritz::Rayleigh_Ritz(int ns)
{
	setParameters(ns);
}

void Rayleigh_Ritz::setParameters(int ns)
{
	n = ns;
	h = 1./(1+n);
}

float Rayleigh_Ritz::eval_aij(float x)
{
	return myProblem.p(x)*myBasis1.dphi(x)*myBasis2.dphi(x)+myProblem.q(x)*myBasis1.phi(x)*myBasis2.phi(x);
}

float Rayleigh_Ritz::Aij(int i, int j, float lim_inf, float lim_sup)
{
  myBasis1.setParameters(i, n);
  myBasis2.setParameters(j, n);
  
  return myIntegral.solve(lim_inf,lim_sup,15,[this](float x)->float {return this->eval_aij(x);});
}


float Rayleigh_Ritz::eval_bi(float x)
{
	return myProblem.f(x)*myBasis1.phi(x);
}

float Rayleigh_Ritz::bi(int i, float lim_inf, float lim_sup)
{
  myBasis1.setParameters(i, n);

  return myIntegral.solve(lim_inf,lim_sup,15,[this](float x)->float {return this->eval_bi(x);});

}


float Rayleigh_Ritz::x_j(int j)
{
	if (j==-2 || j==-1) return 0;
	else if (j==n+2 || j==n+3) return 1;
	else return j*h;	
}

float Rayleigh_Ritz::y(float x, const vector<float> c)
{
	float y_value = 0;

	for (int i = 0; i < n+2; ++i)
	{
        myBasis2.setParameters(i, n);

		y_value += c[i]*myBasis2.phi(x);
	}
	return y_value; 
}

void Rayleigh_Ritz::solve()
{
  vector<vector<float>> A(n+2,vector<float>(n+2));
  vector<float> b(n+2);
  vector<float> c(n+2);



  for (int i = 0; i < n+2; ++i) 
  {

  	  for (int j = i; j <= min(i+3,n+1); ++j) 
	  {

	  	L = max(x_j(j-2), (float)(0));
	  	U = min(x_j(i+2), (float)(1));

	  	A[i][j] = Aij(i,j, L, U);

	  	if (i!=j) {A[j][i] = A[i][j]; }

	  }

	  if (i>=4)
	  {
	  	for (int j = 0; j <= i-4; ++j)
	  	{
	  		A[i][j] = 0;
	  	}
	  }

	  if (i<=n-3)
	  {
	  	for (int j = i+4; j <= n+1; ++j)
	  	{
	  		A[i][j] = 0;
	  	}
	  }

  }


  for (int i = 0; i < n+2; ++i)
  {
	L = max(x_j(i-2),(float)(0));
	U = min(x_j(i+2),(float)(1));
	
	b[i] = bi(i,L,U);

  }

  Cholesky mySisEcua(n+2, A, b, c);

  mySisEcua.solve();


  cout<<"|  i|	c_i	  |	x_i	|   y(x_i)    |y_exacta(x_i)||y(x_i)-y_exacta(x_i)||"<<endl;

  float y_exacta[n+2]={0,0.30901699,0.58778525,0.80901699,0.95105652,1,0.95105652,0.8091699,0.58778525,0.30901699,0};
  for (int i = 0; i < n+2; ++i)
  {
    myBasis1.setParameters(i, n);
    float y_i = y(x_j(i), c);

  	cout<<"|"<<setw(3)<<i;

  	cout<<"|"<<setw(13)<<c[i];
  	cout<<"|"<<setw(13)<<x_j(i);
  	//cout<<y_i<<" ,";

  	cout<<"|"<<setw(13)<<y_i;

  	cout<<"|"<<setw(13)<<y_exacta[i];

  	cout<<"|"<<setw(20)<<abs(y_exacta[i]-y_i)<<" |\n";


  }  
};
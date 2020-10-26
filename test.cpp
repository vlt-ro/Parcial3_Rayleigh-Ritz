/**
 * Programa que prueba el metodo Cholesky
 * @author Valentina Roquemen
 */


#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

//#include "utils/Cholesky/cholesky.h"
//#include "utils/Trapezoide/trapezoide.h" 
//#include "utils/BSpline/b-spline.h" 
//#include "problem.h"
#include "rayleigh-ritz.h"

using namespace std;

/*class f_x: public Function
{
public:
   float eval(float x)
   {
       return 1/x;
   }
};
*/
  float eval(float x)
   {
       return 1/x;
   }
int main()
{ 
/*
  int n = 3;
  vector<vector<float>> A = {{4,-1,1},{-1,4.25,2.75},{1,2.75,3.5}};
  vector<float> b = {1,0,2};
  vector<float> c(n);

  Cholesky myCholesky(n,A,b,c);

  myCholesky.solve();
    
  for (int i = 0; i < n; ++i)
  {
    cout<<c[i]<<endl;
  }
  
  Problem myProblem;

  cout<<myProblem.f(5);
  
  
  //f_x fx;

  float lim_inf = 1, lim_sup = 2.;
  int k = 7;
  float integral;

  MetodoTrapezoide myIntegral;

  integral = myIntegral.solve(lim_inf,lim_sup,k,eval);
  cout<<integral<<endl;

  

    
  B_Spline base;
  base.setParameters(3,5);
  integral = myIntegral.solve(lim_inf,lim_sup,k,&base);

  cout<<integral<<endl;

  cout<<base.eval(0.6)<<endl;
  */

  Rayleigh_Ritz myProblem(20);

  myProblem.solve();
  
  return 0;
}
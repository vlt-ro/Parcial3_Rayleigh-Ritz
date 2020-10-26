/**
 * @author Valentina Roquemen Echeverry
 * @brief  Solucion de sistemas de ecuaciones de la
 *         forma Ac=b usando el metodo de Cholesky
 *         A es una matriz definida positiva
 */
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "utils/Cholesky/cholesky.h"
#include "utils/Trapezoide/trapezoide.h" 
#include "utils/BSpline/b-spline.h" 
#include "problem.h"

using namespace std;

class Rayleigh_Ritz 
{
public:

  Rayleigh_Ritz(int);
  void setParameters(int);
  void solve();
  float Aij(int, int, float, float);
  float bi(int,float, float);
  float x_j(int);
  float y(float, const vector<float>);
  float eval_bi(float);
  float eval_aij(float);

private:

	int n;
	float h;
	float L,U;
	B_Spline myBasis1;
	B_Spline myBasis2;
    Problem myProblem;
    MetodoTrapezoide myIntegral;

};
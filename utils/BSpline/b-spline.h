#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
//#include "../Trapezoide/trapezoide.h"
using namespace std;

#ifndef B_SPLINE_BASIS_H
#define B_SPLINE_BASIS_H

class B_Spline
{
public:

  B_Spline();
  void setParameters(int,int);
  float S(float);
  float dS(float);
  float phi(float);
  float dphi(float);
private:

	float x;
	float h; 
	int n;
	int i;

};

#endif // B_SPLINE_BASIS_H

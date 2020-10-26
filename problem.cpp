#include<cmath>
#include "problem.h"

using namespace std;

Problem::Problem()
{

}

float Problem::p(float x)
{
  return 1.;
}

float Problem::q(float x)
{
  return pow(M_PI,2);
  //return 0;

}

float Problem::f(float x)
{
  return 2*pow(M_PI,2)*sin(M_PI*x);
  //return x*pow(x-1,3);
  //return 3*x;
}
#include<cmath>
#include "problem.h"

Problem::Problem(){}

float Problem::p(float x){return 1.;}

float Problem::q(float x){return pow(M_PI,2);}

float Problem::f(float x){return 2*pow(M_PI,2)*sin(M_PI*x);}
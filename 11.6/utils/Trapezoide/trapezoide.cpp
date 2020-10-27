#include <functional>
#include <cmath>
#include "trapezoide.h"

MetodoTrapezoide::MetodoTrapezoide()
{
  N = 15;
}

void MetodoTrapezoide::setParameters(float lim_inf, float lim_sup, const function<float(float)> & f)
{
  a = lim_inf;
  b = lim_sup;
  function_eval = f;
} 

void MetodoTrapezoide::deltaX_k(int k)
{
  dx = (b-a)/pow(2,k);
}

void MetodoTrapezoide::T_k(int k)
{
  float sumatoria=0;
  deltaX_k(k);

  for (int i = 1; i < pow(2,k); i+=2) sumatoria += function_eval(a+i*dx);

  t_k = 0.5*t_k + dx*sumatoria;
}

void MetodoTrapezoide::integralTrapezoide()
{
  t_k = 0.5*(function_eval(b)+function_eval(b-a))*(b-a); // Calculo de T0

  for (int k = 1; k <= N; ++k) T_k(k); // Se actualiza el valor de la integral

  integral = t_k;
}

float MetodoTrapezoide::solve(float lim_inf, float lim_sup, const function<float(float)> & f)
{
  setParameters(lim_inf,lim_sup,f);
  integralTrapezoide();
  
  return integral;
}

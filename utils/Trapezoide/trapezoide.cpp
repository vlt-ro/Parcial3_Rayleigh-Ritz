/**
 * @author Valentina Roquemen Echeverry
 * @brief  Esta libreria implemeta el metodo trapezoide
 *			par la integracion de funciones definidas. 
 *
 */

  
#include <iostream> 
#include <functional>
#include <iomanip>
#include <cmath>
#include "trapezoide.h"
using namespace std;


MetodoTrapezoide::MetodoTrapezoide()
{
}

void MetodoTrapezoide::establecerParametros(float lim_inf, float lim_sup, int K,const function<float(float)> & f)
{
	a = lim_inf;
	b = lim_sup;
	N = K;
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

	for (int i = 1; i < pow(2,k); i+=2)
	{
		sumatoria += function_eval(a+i*dx);
		
	}
	t_k = 0.5*t_k + dx*sumatoria;
}

void MetodoTrapezoide::integralTrapezoide()
{
	t_k = 0.5*(function_eval(b)+function_eval(b-a))*(b-a); // Calculo de T0

	for (int k = 1; k <= N; ++k)
	{	
		T_k(k); // Se actualiza el valor de la integral
	}

	integral = t_k;
}

float MetodoTrapezoide::solve(float lim_inf, float lim_sup, int K,const function<float(float)> & f)
{
	establecerParametros(lim_inf,lim_sup,K,f);
	integralTrapezoide();
	
    return integral;

}

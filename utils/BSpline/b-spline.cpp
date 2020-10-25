/**
 * @author Valentina Roquemen Echeverry
 * @brief  Solucion de sistemas de ecuaciones de la
 *		     forma Ac=b usando el metodo de Cholesky
 *		     A es una matriz definida positiva
 */

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include<cmath>
#include "b-spline.h"

using namespace std;

B_Spline::B_Spline()
{

}

void B_Spline::setParameters(int i_p, int n_p)
{
	n = n_p;
	i = i_p;
	h =1./(n+1);
}

float B_Spline::S(float y)
{
	if (y<=-2){return 0;}
	else if (y>=-2 && y<=-1){return pow(2+y,3)/4;}
	else if (y>=-1 && y<=0){return (pow(2+y,3)-4*pow(1+y,3))/4;}
	else if (y>-1 && y<=1){return (pow(2-y,3)-4*pow(1-y,3))/4;}
	else if (y>1 && y<=2){return pow(2-y,3)/4;}
	else if (y>2){return 0;}
}

float B_Spline::dS(float y)
{
	if (y<=-2){return 0;}
	else if (y>=-2 && y<=-1){return 3*pow(2+y,2)/4;}
	else if (y>=-1 && y<=0){return 3*(pow(2+y,2)-4*pow(1+y,2))/4;}
	else if (y>-1 && y<=1){return -3*(pow(2-y,2)-4*pow(1-y,2))/4;}
	else if (y>1 && y<=2){return -3*pow(2-y,2)/4;}
	else if (y>2){return 0;}
}

float B_Spline::phi(float x_point)
{
	x = x_point;

	if (i==0){return S(x/h)-4*S((x+h)/h);}
	else if (i==1){return S((x-h)/h)-S((x+h)/h);}
	else if (i>=2 && i<=n-1){return S((x-i*h)/h);}
	else if (i==n){return S((x-n*h)/h)-S((x-(n+2)*h)/h);}
	else if (i==n+1){return S((x-(n+1)*h)/h)-4*S((x-(n+2)*h)/h);}
}

float B_Spline::dphi(float x_point)
{
	x = x_point;
	
	if (i==0){return dS(x/h)-4*dS((x+h)/h);}
	else if (i==1){return dS((x-h)/h)-dS((x+h)/h);}
	else if (i>=2 && i<=n-1){return dS((x-i*h)/h);}
	else if (i==n){return dS((x-n*h)/h)-dS((x-(n+2)*h)/h);}
	else if (i==n+1){return dS((x-(n+1)*h)/h)-4*dS((x-(n+2)*h)/h);}
}


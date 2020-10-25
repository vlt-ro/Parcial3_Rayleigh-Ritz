/**
 * @author Valentina Roquemen Echeverry
 * @brief  Esta libreria implemeta el metodo trapezoide
 *			par la integracion de funciones definidas. 
 *
 */

#include <stdlib.h>
#include <functional>
#include <iostream>

using namespace std;
#ifndef TRAPEZOIDALRULE_H
#define TRAPEZOIDALRULE_H

class MetodoTrapezoide
{
	public :
	/**
     * @param lim_inf:  Limite inferior de la integral
     * @param lim_sup:  Limite superior de la integral
	 * @param k : Potencia maxima de dos sobre la se 
	 *            partira la funcion para la integracion
    */
	MetodoTrapezoide();

	/*
	 * Fija los parametros de a,b,N ingresados por el usuario 
	*/
	void establecerParametros(float, float, int, const function<float(float)> & f);

	/*
	 * Calcula el ancho del trapezoide en la iteracion k
	*/
	void deltaX_k(int);

	/*
	 * Calcula el termino T_k 
	*/
	void T_k(int);

	/*
	 * Realiza la integracion llegando hasta la iteracion N
	*/
	void integralTrapezoide();

	/*
	* 
	*/ 
	float solve(float,float,int,const function<float(float)> & f);

	private :     
	
	float a; // Limite inferior de la integral
	float b; // Limite superior de la integral
	int N; // Numero de iteraciones a realizar
	float dx; // Tamanio de bins
	float t_k; // Area de 2^k trapezoides
	float integral = 0; // Resultado de la integral
	function<float(float)> function_eval;

};

#endif // TRAPEZOIDALRULE_H

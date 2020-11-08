#ifndef TOOLS_H
#define TOOLS_H

#include <vector>

using std::vector;


//template <typename func_type>

///**
// * @brief Integración con el método de simposon
// * @param a Límite inferior
// * @param b Límite superior
// * @param n Número de intervalos
// * @param f Función a integrar
// * @return Integral de la funcion f en el rango [a,b]
// * @ref https://stackoverflow.com/questions/60005533/composite-simpsons-rule-in-c
// */
//double simpson_rule(double a, double b, int n, func_type f);

/**
* @brief
* @param n Dimension de la matriz cuadrada
* @param A Matriz A de dimension nxn
* @param b Vector b de dimension 1xn
* @param c Vector c de dimension 1xn que contendra la
*        solucion.
* @ref http://cod-ayu.blogspot.com/2015/10/solucion-de-sistemas-de-ecuaciones.html
*/
void gauss_jordan(int n, const vector<vector<double>>& A, const vector<double>& b, vector<double>& c);

#endif // TOOLS_H

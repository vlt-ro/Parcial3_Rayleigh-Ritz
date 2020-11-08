/**
 * @author Valentina Roquemen Echeverry
 * @brief  Clase que implementa el metodo trapezoidal para integrar.
 *         Se usa el algoritmo descrito en la sección 14.6 del libro
 *         de Bronson "C++ para ingeniería y ciencias" segunda edición.
 */

#include <stdlib.h>
#include <functional>

using namespace std;

#ifndef TRAPEZOIDALRULE_H
#define TRAPEZOIDALRULE_H

class MetodoTrapezoide
{
public :

  /*
  * @param K : Potencia maxima de dos sobre la se 
  *            partira la funcion para la integracion
  */
  MetodoTrapezoide();

  /*
  * @brief Resuelve la integral
  * @param lim_inf:  Limite inferior de la integral
  * @param lim_sup:  Limite superior de la integral
  * @param f :  Funcion a integrar
  * @return Valor de la integral
  */ 
  float solve(float,float,const function<float(float)> & f);

private :

  float a; // Limite inferior de la integral
  float b; // Limite superior de la integral
  int N; // Numero de iteraciones a realizar
  float dx; // Tamanio de bins
  float t_k; // Area de 2^k trapezoides
  float integral = 0; // Resultado de la integral
  function<float(float)> function_eval; // Funcion a integrar


  /*
   * @brief Fija los parametros de a,b,N ingresados por el usuario
   */
  void setParameters(float, float, const function<float(float)> & f);

  /*
   * @brief Calcula el ancho del trapezoide en la iteracion k
   */
  void deltaX_k(int);

  /*
   * @brief Calcula el termino T_k
   */
  void T_k(int);

  /*
   * @brief Realiza la integracion llegando hasta la iteracion N
   */
  void integralTrapezoide();

};

#endif // TRAPEZOIDALRULE_H

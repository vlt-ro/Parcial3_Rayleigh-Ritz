/**
 * @author Valentina Roquemen Echeverry
 * @brief  Clase en la cual el usario define las funciones p(x),
 *         q(x) y f(x) de la ecuacion diferencial que se va a solucionar
 *         EDO: -d(p(x)dy/dx)/dx+q(x)y=f(x) 0<=x<=1 y(0)=y(1)=0
 */

#ifndef PROBLEM_H
#define PROBLEM_H

class Problem
{
public:

  Problem();

  /**
  * @brief Definicion de la funcion p(x)
  * @param x Punto en el cual se evaluara la funcion base
  * @note p(x)>0 para 0<=x<=1, p E C[0,1]
  */
  float p(float);
  
  /**
  * @brief Definicion de la funcion q(x)
  * @param x Punto en el cual se evaluara la funcion base
  * @note p(x)>=0 para 0<=x<=1, p E C[0,1]
  */
  float q(float);

  /**
  * @brief Definicion de la funcion f(x)
  * @param x Punto en el cual se evaluara la funcion base
  * @note f E C[0,1]
  */
  float f(float);

private:

  float x; // Variable a evaluar
};

#endif // PROBLEM_H
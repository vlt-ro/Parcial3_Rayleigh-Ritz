/**
 * @author Valentina Roquemen Echeverry
 * @brief Clase que implementa la definicion de la base
 *        B-Spline de la seccion 11.5 del libro "Numerical
 *        Analysis" de Richard L. Burden y J. Douglas Faires.
 */

#ifndef B_SPLINE_BASIS_H
#define B_SPLINE_BASIS_H

class B_Spline
{

public:

  B_Spline();

  /**
  * @brief Establece los parametros de la base
  * @param i Elemento i-esimo de la base que se va a tomar
  * @param n La dimension de la base es n+2
  */  
  void setParameters(int,int);

  /**
  * @brief Definicion del polinomio cubico S(x)
  * @param x Punto en el cual se evaluara la funcion base
  * @note Ecuación 11.31
  */
  float S(float);

  /**
  * @brief Definicion de los elementos de la base phi_i(x)
  * @param x Punto en el cual se evaluara la funcion base
  */
  float dS(float);


  float phi(float);
  
  /**
  * @brief Definicion de la primera derivada respecto a x de
           los elementos de la base phi_i(x)
  * @param x Punto en el cual se evaluara la funcion base
  */
  float dphi(float);

private:

  float x; // Variable a evaluar
  float h; // Tamanio de los nodos 
  int n; // La base es de dimension n+2
  int i; // Elemento i-esimo de la base

};

#endif // B_SPLINE_BASIS_H
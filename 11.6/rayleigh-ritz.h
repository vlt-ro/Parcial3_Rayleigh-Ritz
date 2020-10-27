/**
 * @author Valentina Roquemen Echeverry
 * @brief Clase que implementa el metodo de Rayleigh-Ritz usando
 *        la base B-Spline. Se usal el algoritmo 11.6  del libro
 *        "Numerical Analysis" de Richard L. Burden y J. Douglas Faires.
 */
#include <stdlib.h>
#include <vector>

#include "utils/Trapezoide/trapezoide.h" 
#include "utils/BSpline/b-spline.h" 
#include "problem.h"

using namespace std;

class Rayleigh_Ritz 
{
public:

  /*
  * @param n : Se expandira la base de dimension n+2
  */
  Rayleigh_Ritz(int);

  /**
  * @brief Establece el parametro n
  * @param n 
  */
  void setParameters(int);

  /**
  * @brief Define el nodo j-esimo
  * @param j 
  */
  float x_j(int);

  /**
  * @brief Expande la funcion y(x) en la base BSpline
  * @param x Punto en el cual se evaluara la funcion base
  * @param c : Vector con los coeficientes de la funcion y(x) 
  */
  float y(float, const vector<float>);

  /**
  * @brief Calcula la componente Aij de la matriz A
  * @param i : Fila
  * @param j : Columna
  * @param L : Limite inferior de integracion
  * @param U : Limite superior de integracion
  */
  float Aij(int, int, float, float);

  /**
  * @brief Evalua la funcion del integrando
  * @param x : Punto en el cual se evaluara la funcion 
  * @note Ecuacion 11.28
  */
  float eval_aij(float);

  /**
  * @brief Calcula la componente bi del vector b
  * @param i : Fila
  * @param L : Limite inferior de integracion
  * @param U : Limite superior de integracion
  */
  float bi(int,float, float);

  /**
  * @brief Evalua la funcion del integrando
  * @param x : Punto en el cual se evaluara la funcion 
  * @note Ecuacion 11.28
  */
  float eval_bi(float);
  
  /*
  * @brief Se resuelve la ecuacion diferencial
  */
  void solve();

private:

  int n; // La dimension de la base es de n+2
  float h;
  float L,U; // Limite inferior (L) y superior (U) de la integral
  B_Spline Basis_i, Basis_j; // Clase que define la base
  Problem myProblem; // Clase que define el problema
  MetodoTrapezoide Integral; // Clase que realiza la integracion

};
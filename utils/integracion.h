#ifndef INTEGRACION_H
#define INTEGRACION_H

#include <iostream>

using namespace std;

namespace integracion
{

  /**
   * @brief Integración con el método de simposon
   * @param a Límite inferior
   * @param b Límite superior
   * @param n Número de intervalos
   * @param f Función a integrar
   * @return Integral de la funcion f en el rango [a,b]
   * @ref https://stackoverflow.com/questions/60005533/composite-simpsons-rule-in-c
   */
  template <typename func_type>
  double simpson_rule(double a, double b, int n, func_type f)
  {
    double h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    double sum_odds = 0.0;
    for (int i = 1; i < n; i += 2) { sum_odds += f(a + i * h);}
  
    double sum_evens = 0.0;
    for (int i = 2; i < n; i += 2) { sum_evens += f(a + i * h);}

    return (f(a) + f(b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
  }
}

#endif // INTEGRACION_H
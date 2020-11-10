#ifndef RAYLEIGHRITZ_H
#define RAYLEIGHRITZ_H
#include <vector>
using std::vector;

class RayleighRitz
{
public:

    /**
     * @brief Resuelve el problema
     * @return Vector con los coeficientes de expansión 
     */
    virtual vector<double> solve(double(*p)(double),
                         double(*q)(double),
                         double(*f)(double)) = 0;

    /**
     * @brief Evalua la función en un punto x
     * @return Función evaluada en un púnto específico
     */
    virtual double eval(double) = 0;


};

#endif // RAYLEIGHRITZ_H


#ifndef CSRAYLEIGHRITZ_H
#define CSRAYLEIGHRITZ_H

#include <vector>
#include "bslipne.h"

using std::vector;

class Basis;

class csRayleighRitz
{
public:
    csRayleighRitz(std::size_t n);

    /**
     * @brief Resuelve de forma aproximada la EDO
     *        -D( p(x) D(y) ) + q(x)y = f(x)
     *        con condiciones de frontera y(0)=y(1)=0
     * @param n
     *
     * @return
     */
    vector<double> solve(double(*p)(double),
                         double(*q)(double),
                         double(*f)(double));

    vector<Basis>& getBasis();

private:
    vector<Basis> phi; //B-Spline basis
    std::size_t n;
    double h;

    /**
     * Retorna el valor del j-iesimo nodo
     */
    double x_i(int i);
};


/**
 * Clase con la definición de la base de las
 * funciones usadas en el método lineal por partes
 * de Ryleigh-Ritz.
 */
class Basis
{
public:

    /**
     * @brief Constructor para inicializar un elemento de la base
     */
    Basis(int i, int n, double h);

    ~Basis();

    double operator()(double x);

    double dPhi(double);

private:
    Phi *phi;
};

#endif // CSRAYLEIGHRITZ_H

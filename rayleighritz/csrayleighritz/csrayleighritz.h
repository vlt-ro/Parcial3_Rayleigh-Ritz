#ifndef CSRAYLEIGHRITZ_H
#define CSRAYLEIGHRITZ_H

#include <vector>
#include "bspline.h"
#include "../rayleighritz.h"
using std::vector;


class csBasis;

class csRayleighRitz : public RayleighRitz
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
    ~csRayleighRitz(); // destructor
    vector<double> solve(double(*p)(double),
                         double(*q)(double),
                         double(*f)(double));

    vector<csBasis>& getBasis();

    double eval(double);

private:
    vector<csBasis> phi; //B-Spline basis
    vector<double> c; //Coeficientes de expansión
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
class csBasis
{
public:

    /**
     * @brief Constructor para inicializar un elemento de la base
     */
    csBasis(int i, int n, double h);

    /**
     * Sobrecarga del constructor de copia
     */
    csBasis(const csBasis &);

    ~csBasis();

    double operator()(double x);

    double dPhi(double);

private:
    Phi *phi;
    int i;
    int n;
    double h;

    void setMembers(int i, int n, double h);
};

#endif // CSRAYLEIGHRITZ_H

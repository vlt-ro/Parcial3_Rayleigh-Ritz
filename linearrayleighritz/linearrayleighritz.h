/**
 * @author Santiago Duque
 * @brief Clase que implementa el método de Rayleigh-Ritz mostrado
 *        en el algoritmo 11.5, de la sección 11.5 del libro "Numerical
 *        Analysis" de Richard L. Burden y J. Douglas Faires.
 */

#ifndef LINEARRAYLEIGHRITZ_H
#define LINEARRAYLEIGHRITZ_H

#include <vector>

class Basis;

class LinearRayleighRitz
{
public:

    LinearRayleighRitz(std::vector<double>& x);

    /**
     * @brief Resuelve de forma aproximada la EDO
     *        -D( p(x) D(y) ) + q(x)y = f(x)
     *        con condiciones de frontera y(0)=y(1)=0
     * @param x Arreglo entre [0,1] con los valores en x
     *          donde se calculará la aproximación a la función
     * @return
     */
    std::vector<double> solve(
                             double(*p)(double),
                             double(*q)(double),
                             double(*f)(double));

    /**
     * Retorna la base que expande la aproximación de la
     * solución.
     */
    std::vector<Basis> &getBasis();

private:
    std::vector<double> h; // Coeficientes h_i = x_{i+1} - x_i
    std::vector<Basis> phi; // Base lineal a trozos
    std::vector<double>& x;

    /*
     * Integrales a resolver (página 700)
     */
    double Q1(int , std::vector<double>& , double(*)(double));
    double Q2(int , std::vector<double>& , double(*)(double));
    double Q3(int , std::vector<double>& , double(*)(double));
    double Q4(int , std::vector<double>& , double(*)(double));
    double Q5(int , std::vector<double>& , double(*)(double));
    double Q6(int , std::vector<double>& , double(*)(double));
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
     * @brief Constructor para inicializar los parametros
     *        de un elemento de la base específico.
     * @param x_im1 Punto x_{i-1}
     * @param x_i   Punto x_i
     * @param x_ip1 Punto x_{i+1}
     */
    Basis(double x_im1, double x_i,double x_ip1);

    /**
     * @brief Sobrecarga del operator () para usar los objetos de
     *        esta clase como "functors"
     * @param x Punto en el cual se evaluará la función base     *
     * @note Ecuación 11.30
     */
    double operator()(double x);

private:
    /* Puntos que definen la función por partes */
    double x_im1, x_i, x_ip1;
};

#endif // LINEARRAYLEIGHRITZ_H

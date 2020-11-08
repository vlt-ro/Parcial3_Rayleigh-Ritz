#include "linearrayleighritz.h"
#include <iostream>

using std::vector;


namespace
{
    template <typename func_type>
    /**
     * @brief simpson_rule
     * @param a  Límite inferior
     * @param b  Límite superior
     * @param n  Número de intervalos
     * @param f  función que se quiere integrar
     * @return resultado de integrar la función en el rango [a,b]
     * @ref https://stackoverflow.com/questions/60005533/composite-simpsons-rule-in-c
     */
    double simpson_rule(double a, double b,
                        int n, // Number of intervals
                        func_type f)
    {
        double h = (b - a) / n;

        // Internal sample points, there should be n - 1 of them
        double sum_odds = 0.0;
        for (int i = 1; i < n; i += 2)
        {
            sum_odds += f(a + i * h);
        }
        double sum_evens = 0.0;
        for (int i = 2; i < n; i += 2)
        {
            sum_evens += f(a + i * h);
        }

        return (f(a) + f(b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
    }
}


Basis::Basis(double x_im1, double x_i,double x_ip1)
{
    this->x_im1 = x_im1;
    this->x_i = x_i;
    this->x_ip1 = x_ip1;
}

double Basis::operator()(double x)
{
    double result = 0;

    if( x_im1 < x && x <= x_i )
    {
        double h = x_i - x_im1;
        result = (x - x_im1) / h;
    }
    else if(x_i < x && x < x_ip1)
    {
        double h = x_ip1 - x_i;
        result = (x_ip1 - x) / h;
    }

    return result;
}


std::vector<double> LinearRayleighRitz::solve(std::vector<double> &x,
                                             double(*p)(double),
                                             double(*q)(double),
                                             double(*f)(double))
{
    int n = x.size() - 2;

    /*******************************************/
    /* Step 1: Cálculo de los coeficientes h_i */
    /*******************************************/
    for(int i=0; i<=n; ++i)
        h.push_back(x[i+1]-x[i]);

    /***************************************************/
    /* Step 2: Definir la base lineal por tramos phi_i */
    /***************************************************/
    for(int i=1; i <= n; ++i)
        phi.push_back(Basis(x[i-1],x[i], x[i+1]));

    /*************************************/
    /* Step 3: Calculo de las integrales */
    /*************************************/

    /* Aunque es un desperdicio de recursos, se va a
     * definir la posición 0 de los vectores que almacenan
     * el resultado de las integrales, con el fin de seguir
     * el algoritmo y no tener problemas con los índices usados
     * en el libro. La idea es facilitar la lectura del código
     * y la comparación con el libro. */
    vector<double> q1(1),q2(1),q3(1),q4(1),q5(1),q6(1);
    for(int i=1; i<=n; ++i)
    {
        if(i<n)
            q1.push_back(Q1(i,x,q));
        q2.push_back(Q2(i,x,q));
        q3.push_back(Q3(i,x,q));
        q4.push_back(Q4(i,x,p));
        q5.push_back(Q5(i,x,f));
        q6.push_back(Q6(i,x,f));
    }
    q4.push_back(Q4(n+1,x,p));


    /***************************************************************/
    /* Step 4 and 5: Calculo de los factores alpha_i, beta_i y b_i */
    /***************************************************************/
    // Se define el primer elemento, por la misma razón que
    // el valor de las integrales
    vector<double> alpha(1), beta(1), b(1);
    for(int i=1; i<=n; ++i)
    {
        alpha.push_back(q4[i] + q4[i+1] + q2[i] + q3[i]);
        b.push_back(q5[i] + q6[i]);

        if(i < n)
            beta.push_back(q1[i] - q4[i+1]);
    }

    /**************************************************************/
    /* Step 6, 7 and 8: Calculo de los factores a_i, zeta_i y z_i */
    /**************************************************************/
    vector<double> a(2), zeta(2), z(2);
    a[1] = alpha[1];
    zeta[1] = beta[1] / alpha[1];
    z[1] = b[1] / a[1];

    for(int i=2; i<=n; ++i)
    {
        a.push_back(alpha[i] - beta[i-1]*zeta[i-1]);
        z.push_back( (b[i] - beta[i-1]*z[i-1])/a[i] );

        if(i<n)
            zeta.push_back(beta[i]/a[i]);
    }

    /**************************************************/
    /* Step 9 and 10: Calculo de los factores c_i */
    /**************************************************/
    /* Coeficientes de la expansión en series*/
    vector<double> c(n); //c_i con i=0,...,n-1
    c[n-1] = z[n];
    for(int i=n-1; i>0; --i)
        c[i-1] = z[i] - zeta[i]*c[i];

    return c;
}

std::vector<Basis>& LinearRayleighRitz::getBasis()
{
    return phi;
}

double LinearRayleighRitz::Q1(int i, std::vector<double>& x, double(*q)(double))
{
    auto f = [&x,q,i](float _x)->float{return q(_x)*(x[i+1] - _x)*(_x - x[i]);};
    auto r = simpson_rule(x[i],x[i+1],10,f);
    return r/(h[i]*h[i]);
}

double LinearRayleighRitz::Q2(int i, std::vector<double>& x, double(*q)(double))
{
    auto f = [&x,q,i](float _x)->float{return q(_x)*(_x - x[i-1])*(_x - x[i-1]);};
    auto r = simpson_rule(x[i-1],x[i],10,f);
    return r/(h[i-1]*h[i-1]);
}

double LinearRayleighRitz::Q3(int i, std::vector<double>& x, double(*q)(double))
{
    auto f = [&x,q,i](float _x)->float{return q(_x)*(x[i+1]-_x)*(x[i+1]-_x);};
    auto r = simpson_rule(x[i],x[i+1],10,f);
    return r/(h[i]*h[i]);
}

double LinearRayleighRitz::Q4(int i, std::vector<double>& x, double(*p)(double))
{
    auto f = [p](float _x)->float{return p(_x);};
    auto r = simpson_rule(x[i-1],x[i],10,f);
    return r/(h[i-1]*h[i-1]);
}

double LinearRayleighRitz::Q5(int i, std::vector<double>& x, double(*f)(double))
{
    auto fun = [&x,f,i](float _x)->float{return f(_x)*(_x-x[i-1]);};
    auto r = simpson_rule(x[i-1],x[i],10,fun);
    return r/h[i-1];
}

double LinearRayleighRitz::Q6(int i, std::vector<double>& x, double(*f)(double))
{
    auto fun = [&x,f,i](float _x)->float{return f(_x)*(x[i+1]-_x);};
    auto r = simpson_rule(x[i],x[i+1],10,fun);
    return r/h[i];
}

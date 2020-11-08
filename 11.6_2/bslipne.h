#ifndef BSLIPNE_H
#define BSLIPNE_H


class Phi
{
public:

    Phi(double h=0);

    /**
     * sobrecarga del operador ()
     */
    virtual double operator()(double) = 0;

    /**
     * @brief Derivada de la funcion
     * @return Derivada de la función evaluada en un púnto específico
     */
    virtual double DPhi(double) = 0;

protected:
    double h;
};


class Phi_0:public Phi
{
public:
    Phi_0(double h);
    double operator()(double);
    double DPhi(double);
};


class Phi_1:public Phi
{
public:
    Phi_1(double h);
    double operator()(double);
    double DPhi(double);
};

class Phi_i:public Phi
{
public:
    Phi_i(double h, int i);
    double operator()(double);
    double DPhi(double);
private:
    int i;
};

class Phi_n:public Phi
{
public:
    Phi_n(double h, int n);
    double operator()(double);
    double DPhi(double);
private:
    int n;
};

class Phi_n_1:public Phi
{
public:
    Phi_n_1(double h, int n);
    double operator()(double);
    double DPhi(double);
private:
    int n;
};

#endif // BSLIPNE_H

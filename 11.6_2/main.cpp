#include <iostream>
#include "csrayleighritz.h"
#include <vector>
#include <cmath>

#include <fstream>

using namespace std;

double p(double)
{
    return 1;
}

double q(double)
{
    return pow(M_PI,2);
}

double f(double x)
{
    return 2*pow(M_PI,2)*sin(M_PI*x);
}

double expand(double x, vector<double>& c, vector<Basis>&basis)
{
    double rs = 0;

    for(size_t i=0; i<c.size(); ++i)
        rs += c[i] * basis[i](x);

    return rs;
}

int main()
{
    vector<double>x;
    int n = 9;
    for(int i=0; i<=n+1; ++i)
        x.push_back((1.0/(float)(n+1))*i);

    csRayleighRitz rr(9);
    auto c = rr.solve(p, q, f);
    auto phi = rr.getBasis();

    cout << "i\tc\tc/2" << endl;
    for(size_t i=0; i<c.size(); ++i)
        cout << i << "\t" << c[i] << "\t" << c[i]/2 << endl;

    cout << "x\tphi" << endl;
    for(size_t i=0; i<x.size(); ++i)
        cout << x[i] << ",\t" << expand(x[i], c, phi) << "," << endl;



    return 0;
}

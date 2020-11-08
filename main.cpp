#include <iostream>
#include "csrayleighritz/csrayleighritz.h"
#include "linearrayleighritz/linearrayleighritz.h"
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

int main()
{
    vector<double>x;
    int n = 9;
    for(int i=0; i<=n+1; ++i)
        x.push_back((1.0/(float)(n+1))*i);

//    csRayleighRitz rr(9);
    LinearRayleighRitz rr(x);
    auto c = rr.solve(p, q, f);

    cout << "i\tc\t" << endl;
    for(size_t i=0; i<c.size(); ++i)
        cout << i << "\t" << c[i] <<  endl;

    cout << "x\tphi" << endl;
    for(size_t i=0; i<x.size(); ++i)
        cout << x[i] << ",\t" << rr.eval(x[i]) << "," << endl;

    return 0;
}

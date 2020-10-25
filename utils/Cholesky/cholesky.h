/**
 * @author Valentina Roquemen Echeverry
 * @brief  Solucion de sistemas de ecuaciones de la
 *         forma Ac=b usando el metodo de Cholesky
 *         A es una matriz definida positiva
 */
#include <stdlib.h>
#include <vector>

using namespace std;

#ifndef CHOLESKY_H
#define CHOLESKY_H

class Cholesky
{
public:

  Cholesky(int, const vector<vector<float>>&, const vector<float>&, vector<float>&);
  void factorL();
  void solve();

private:

  const vector<vector<float>>* A;
  const vector<float>* b;
  vector<float>* c;
  vector<vector<float>> LPtr;

  int dim ;
  float sum;
};

#endif // CHOLESKY_H

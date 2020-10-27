#include <vector>
#include<cmath>

#include "cholesky.h"

Cholesky::Cholesky(int n, const vector<vector<float>>& A, const vector<float>& b, vector<float>& c)
{
  this->A = &A;
  this->b = &b;
  this->c = &c;
  dim = n;
}

void Cholesky::factorL()
{
  vector<vector<float>> L(dim,vector<float>(dim));

  L[0][0] = sqrt((*A)[0][0]);

  for (int j = 1; j < dim; ++j) L[j][0] = (*A)[j][0]/L[0][0];

  for (int i = 1; i < dim-1; ++i)
  {
    sum = 0;

    for (int k = 0; k <= i-1; ++k) sum += pow(L[i][k],2);

    L[i][i] = sqrt((*A)[i][i]-sum);
  }

  for (int i = 1; i < dim-1; i++)
  {
    for (int j = i+1; j < dim; ++j)
    {
      sum = 0;

      for (int k = 0; k <= i-1; ++k) sum += L[j][k]*L[i][k];
      
      L[j][i] = ((*A)[j][i]-sum)/L[i][i];
    }
  }

  sum = 0;

  for (int k = 0; k < dim-1; ++k) sum += pow(L[dim-1][k],2);
    
  L[dim-1][dim-1] = sqrt((*A)[dim-1][dim-1]-sum);

  this->L = L;
}

void Cholesky::solve()
{
  factorL();
  
  vector<float> y(dim);

  y[0] = (*b)[0]/L[0][0];

  for (int i = 1; i < dim; ++i)
  {
    sum = 0;

    for (int j = 0; j <= i-1; ++j) sum += L[i][j]*y[j];
  
    y[i] = ((*b)[i]-sum)/L[i][i];  
  }

  (*c)[dim-1] = y[dim-1]/L[dim-1][dim-1];

  for (int i = dim-2; i>=0; --i)
  {
    sum = 0;

    for (int j = i+1; j < dim; ++j) sum += (*c)[j]*L[j][i];
    
    (*c)[i] = (y[i]-sum)/L[i][i];
  }
};
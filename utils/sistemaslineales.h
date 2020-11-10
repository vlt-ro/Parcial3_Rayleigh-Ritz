#ifndef SISTEMASLINEALES_H
#define SISTEMASLINEALES_H

#include <vector>

using std::vector;

namespace sistemas_lineales

{
  /**
   * @brief
   * @param n Dimension de la matriz cuadrada
   * @param A Matriz A de dimension nxn
   * @param b Vector b de dimension 1xn
   * @param c Vector c de dimension 1xn que contendra la
   *        solucion.
   * @ref http://cod-ayu.blogspot.com/2015/10/solucion-de-sistemas-de-ecuaciones.html
   */
  void gauss_jordan(int n,  vector<vector<double>>& A,  vector<double>& b, vector<double>& c)
  {
    vector<vector<double>> M(n,vector<double>(n+1)); // Matriz aumentada (A|b)
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n+1; ++j)
      {
        if (j==n) { M[i][j] = b[i];}
        else { M[i][j]=A[i][j];}
      }
    }
    
    double may;//variable para almacenar el mayor de la columna k
    int ind;//indice del mayor-->indice de may
    double aux;
    double pivote;
      
    for(int k=0;k<n;k++) //recorrer columnas de la matriz reducida
    {
      may=abs(M[k][k]);
      ind=k;
      
      for(int l=k+1;l<n;l++) //recorrer filas de la columna k para buscar el indice del mayor
      {
        if(may<abs(M[l][k]))
        {
          may=abs(M[l][k]);
          ind=l;
        }
      }
      
      //cambiar filas
      if(k!=ind) 
      {
        for(int i=0;i<n+1;i++)
        {
          aux=M[k][i];
          M[k][i]=M[ind][i];
          M[ind][i]=aux;
        }
      }
      
      if(M[k][k]==1e-14)
      {
        cout<<"No tiene solucion";
        break;
      }
      else
      {
        for(int i=0;i<n;i++) //recorrer fila
        {
          if(i!=k)
          {
            pivote=-M[i][k];
            for(int j=k;j<n+1;j++) { M[i][j]=M[i][j]+pivote*M[k][j]/M[k][k];}
          }
          else
          {
            pivote=M[k][k];
            for(int j=k;j<n+1;j++) { M[i][j]=M[i][j]/pivote;}
          }
        }
      }
    }

    for(int i=0;i<n;i++){ c[i] = M[i][n];}
  }
}

#endif // SISTEMASLINEALES_H
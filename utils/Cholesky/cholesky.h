/**
 * @author Valentina Roquemen Echeverry
 * @brief  Clase que soluciona un sistemas de ecuaciones de la
 *         forma Ac=b usando el metodo de Cholesky del algoritmo
 *         6.6 del libro "Numerical Analysis" de Richard L. Burden 
 *         y J. Douglas Faires.
 */
#include <stdlib.h>
#include <vector>

using namespace std;

#ifndef CHOLESKY_H
#define CHOLESKY_H

class Cholesky
{
public:

  /**
  * @brief Establece la matriz A y los vectores c y b
  * @param n Dimension de la matriz cuadrada
  * @param A Matriz A de dimension nxn
  * @param b Vector b de dimension 1xn
  * @param c Vector c de dimension 1xn que contendra la
  *        solucion.
  * @note A es definida positiva
  */ 
  Cholesky(int, const vector<vector<float>>&, const vector<float>&, vector<float>&);

  /**
  * @brief Encuentra la matriz de factorizacion L de A = L_t*L
  */ 
  void factorL();

  /**
  * @brief Usando L encuentra la solucion del sistema de ecuaciones
  *        y almacena por referencia la solucion en el vector c
  */ 
  void solve();

private:

  const vector<vector<float>>* A; //Puntero a la matriz A
  const vector<float>* b; // Puntero al vector b
  vector<float>* c; // Puntero al vector c
  vector<vector<float>> L; // Matriz L

  int dim ; // Dimension de la matriz
  float sum; //Variable temporal que sera usada en los ciclos
};

#endif // CHOLESKY_H
/**
 * @author Valentina Roquemen
 * @brief Programa que prueba el metodo de Rayleigh-Ritz
 */

#include "rayleigh-ritz.h"

int main()
{ 
  Rayleigh_Ritz myProblem(9);

  myProblem.solve();
  
  return 0;
}
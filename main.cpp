#include "rayleighritz/csrayleighritz/csrayleighritz.h"
#include "rayleighritz/linearrayleighritz/linearrayleighritz.h"
#include "rayleighritz//rayleighritz.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

double p(double) { return 1;}

double q(double) { return pow(M_PI,2);}

double f(double x) { return 2*pow(M_PI,2)*sin(M_PI*x);}

int main()
{
  char choose;
  int n=9;
  vector<double> x;
  vector<double> y;
  
  RayleighRitz *rr;

  cout<<"------------------------------------------------------------\n";
  cout<<"Se va a solucionar una ecuación diferencial de la forma:"<<endl;
  cout<<"-D( p(x) D(y) ) + q(x)y = f(x) con condiciones de frontera y(0)=y(1)=0"<<endl;
  cout<<"Se usa el método de Rayleigh-Ritz"<<endl;
  cout<<"-------------------------------------------------------------\n"; 
  cout<<"Elija la base con la que quiere que se solucione el problema:"<<endl;
  cout<<"1. Lineal a tramos.\n";
  cout<<"2. Cúbica spline.\n";
  cin>>choose;
  cout<<"-------------------------------------------------------------\n"; 
  cout<<"Ingrese el número de nodos que quiere que tenga su solución:";
  cin.clear(); // Se limpia lo que pueda haber en el buffer
  cin.ignore(10000, '\n');
  cin>>n;
  switch (choose)
  {
    case '1':

        for(int i=0; i<=n+1; ++i)
          x.push_back((1.0/(float)(n+1))*i);

        rr = new LinearRayleighRitz(x);
      
      break;


    case '2':
      
        rr = new csRayleighRitz(n);
      
      break;
      
    default:
      cout<<"Por favor ingresa una opcion valida\n";
      exit(1);
  }

  rr->solve(p, q, f);

  cout<<"-------------------------------------------------------------\n"; 
  cout<<"Elija cómo desea obtener los resultados:"<<endl;
  cout<<"1. Guardar en un archivo.\n";
  cout<<"2. Imprimir en pantalla.\n";
  cin>>choose;
  cout<<"-------------------------------------------------------------\n"; 
  cout<<"Ingrese el número de puntos que quiere que tenga su solución:";
  cin.clear(); // Se limpia lo que pueda haber en el buffer
  cin.ignore(10000, '\n');
  cin>>n;
      
  for(int i=0; i<n; ++i)
    y.push_back((1.0/(float)(n-1))*i);
      
  switch (choose)
  {
    case '1':
    {
      ofstream file;
      //file.open ("resultado.csv");
      file.open ("resultado_lineal5.csv");
      file <<"#x,y\n";
      for(size_t i=0; i<y.size(); ++i)
      {
        if (i<y.size()-1)
          file << y[i] << ",\t" << rr->eval(y[i]) << "\n";
        else
          file << y[i] << ",\t" << rr->eval(y[i]);
      }
      file.close();
      
      cout << "Se ha guardado en 'resultado.csv'" << endl;
    }
      break;

    case '2':
      
      cout<<  "---------------------"<<endl;
      cout << "|    x    |    y    |" << endl;
      cout<<  "---------------------"<<endl;

      for(size_t i=0; i<y.size(); ++i)
      {
        cout<<"|"<<setw(9)<<y[i]<<"|"; 
        cout<<setw(9)<<rr->eval(y[i])<<"|\n"; 
      }
      cout<<  "---------------------"<<endl;
      
      break;
      
    default:
      cout<<"Por favor ingresa una opcion valida\n";
      exit(1);
  }

  return 0;
}
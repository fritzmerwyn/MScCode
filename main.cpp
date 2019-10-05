/* This program *** calculates *** A(p^2) through DSE (iterative Eqs.) and then uses it
    to calculate the mass through DSE.
*/

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "definitions.h"
#include "help_numerics.h"
#include "DysonSchwinger.h"




int main(){
  std::cout.precision(17); // set precision to 17 values after comma

  // ##### Constants ###### //

  double m_g=132.0;
  double epsilon=10e-9;
  // double epsilon2=10e-2;
  double m_c = 1e-200;
  // double m_c_max=10.0;
  // double m_g_max= 250.0;
  // double m_g_min= 120.0;
  // double h=(m_c_max-m_c)/max_step;
  // double h=(m_g_max-m_g_min)/max_step;
  // double c_1;
  // c_1 = 1.0/(3.0*pow(m_g,2.0)*pow(pi,2.0));

  // double angint2;

  // Weights and abscissas are generated by the gauleg function in help_numerics.h. gauleg then returns a 2d array, which is split below.
  // The gauleg function is called twice, thought, because of angular integration *and* momentum integration.

  double** x_and_w = nullptr;
  double** x_and_w_ang = nullptr;
  double* absciss_x = nullptr;
  double* weights_w = nullptr;
  double* absciss_ang = nullptr;
  double* weights_ang = nullptr;

  double** vals = nullptr;

  // Splitting 2d array up (momentum)
  x_and_w = gauleg(log(LAMBDA_MIN*LAMBDA_MIN),log(LAMBDA*LAMBDA), absciss_points);
  absciss_x = x_and_w[0]; // x_and_w[0];
  weights_w = x_and_w[1]; // x_and_w[1];

  // Splitting 2d array up (angular)
  x_and_w_ang = gauleg(-1.0,1.0, absciss_points);
  absciss_ang = x_and_w_ang[0];
  weights_ang = x_and_w_ang[1];

  // Save A and B to array "vals"
  vals = iterate_dressing_functions(epsilon,m_c,m_g,absciss_x,weights_w,absciss_ang,weights_ang);
  // angint2 =  qgaus1(angkern2, absciss_ang, weights_ang);

  // ProgressBar pb(max_step, "Doing stuff");

  // Print A (= vals[0][i]) and B (= vals[1][i]) values saved in "vals"
  for(int i=1;i<=absciss_points;i++)
  {
    std::cout<< i << "\t"<< vals[0][i] << "\t" << vals[1][i]<< "\t" << exp(absciss_x[i])<< std::endl;
  }

  // Save A and B values from "vals" to File

  std::ofstream  fileout;
  fileout.open("DressingFunctions_A_and_B_mc_1e-200.dat");
  for(int j=1;j<=absciss_points;j++){
    fileout<< exp(absciss_x[j]) << " " << vals[0][j] << " " << vals[1][j] << std::endl;
  }
  fileout.close ();

  return 0;
}

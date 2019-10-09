
// This program *** calculates *** A(p^2) and B(p^2) from the DSE (iterative Eqs.) and then uses it
// to calculate the mass, which is B(p^2)/A(p^2).


#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "definitions.h"
#include "help_numerics.h"
#include "DysonSchwinger.h"
#include "progressbar.hpp"




int main(){
  std::cout.precision(17); // ##### set precision to 17 values after comma ######

  // ##### Constants ###### //

  double m_g=0.132;
  double epsilon=10e-9;
  // double epsilon2=10e-2;
  double m_c = 5e-3;
  double eta = 1.8;
  double g_squared = 1.0;
  // double m_c_max=10.0;
  // double m_g_max= 250.0;
  // double m_g_min= 120.0;
  // double h=(m_c_max-m_c)/max_step;
  // double h=(m_g_max-m_g_min)/max_step;
  // double c_1;
  // c_1 = 1.0/(3.0*pow(m_g,2.0)*pow(pi,2.0));

  // double angint2;

  // ##### Weights and abscissas are generated by the gauleg function in help_numerics.h. gauleg then returns a 2d array, which is split below. #####
  // ##### The gauleg function is called twice, thought, because of angular integration *and* momentum integration. #####

  double** x_and_w = nullptr;
  double** x_and_w_ang = nullptr;
  double* weights_w = nullptr;
  double* absciss_x = nullptr;
  double* absciss_ang = nullptr;
  double* weights_ang = nullptr;

  double** vals = nullptr;

  #ifdef MarisTandy
    std::cout<<std::endl<< "Now Using Maris Tandy Model!" << std::endl;
  #else
    std::cout<<std::endl<< "Now Using Contact Model!" << std::endl;
  #endif

  // ##### Creating abscissas and weights. Using a Log-grid. #####
  #ifdef loggrid

  std::cout<<std::endl<< "Log-Grid!" << std::endl << std::endl;
  x_and_w = gauleg(log(LAMBDA_MIN*LAMBDA_MIN),log(LAMBDA_squared), absciss_points);

  #else

  std::cout<<std::endl<< "Linear-Grid!" << std::endl << std::endl;
  x_and_w = gauleg(0.0,LAMBDA, absciss_points);

  #endif

  // ##### Splitting 2d array up (momentum) #####
  absciss_x = x_and_w[0]; // x_and_w[0];
  weights_w = x_and_w[1]; // x_and_w[1];

  x_and_w_ang = gauleg(-1.0,1.0, absciss_points);

  // ##### Splitting 2d array up (angular) #####

  absciss_ang = x_and_w_ang[0];
  weights_ang = x_and_w_ang[1];

  // for(int i=1;i<=absciss_points;i++)
  // {
  //   std::cout<< i << "\t"<< weights_w[i] << "\t" << weights_ang[i] << std::endl;
  // }

  // ##### Save A and B to array "vals" #####
  vals = iterate_dressing_functions(epsilon,m_c,m_g,absciss_x,weights_w,absciss_ang,weights_ang, g_squared, eta);

  // ##### Preparing for mass function M(p^2) = B(p^2)/A(p^2). #####
  double* a_vals = vals[0];
  double* b_vals = vals[1];
  double* m_vals = nullptr;
  m_vals = new double[absciss_points];

  for(int i = 1; i<= absciss_points; i++){
    if(a_vals[i] == 0.0){
      m_vals[i] = 0.0;
    }
    else{
      m_vals[i] = b_vals[i]/a_vals[i];
    }
  }

  // ProgressBar pb(max_step, "Doing stuff");

  // ##### Print A (= vals[0][i]) and B (= vals[1][i]) values saved in "vals" #####
  for(int i=1;i<=absciss_points;i++)
  {
    std::cout<< i << "\t"<< vals[0][i] << "\t" << vals[1][i]<< "\t" << m_vals[i] << "\t" << absciss_x[i] << std::endl;
  }

  // ##### Save A and B values from "vals" to File #####
  ProgressBar sd(absciss_points, "Saving Data to File");
  std::ofstream  fileout;
  fileout.open("Data/DressingFunctions_A_and_B_MarisTandy_new_gsq_1.0_mc78_log.dat");
  for(int j=1;j<=absciss_points;j++){
    ++sd;
    fileout<< exp(absciss_x[j]) << " " << vals[0][j] << " " << vals[1][j] << " " << m_vals[j]<< std::endl;
  }
  fileout.close ();

  return 0;
}

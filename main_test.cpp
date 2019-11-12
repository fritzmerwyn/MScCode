
// This program *** calculates *** A(p^2) and B(p^2) from the DSE (iterative Eqs.) and then uses it
// to calculate the mass, which is B(p^2)/A(p^2).


#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include "definitions.h"
#include "help_numerics.h"
#include "Dyson_test.h"
#include "progressbar.hpp"
#include <complex>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>




int main(){
  std::cout.precision(17); // ##### set precision to 17 values after comma ######

  // ##### Constants ###### //
 // std::ofstream gmorpionmass;
 // gmorpionmass.open("Data/gmorpionmass.dat");
 // gmorpionmass<<"#Stromquarkmass and Pionmass"<<std::endl;
 // for(int i=0; i<10;i++){
  double m_g=0.132;
  double epsilon=1e-10;
  double m_c = 0.0037;
  double eta = 1.8;
  double g_squared = 1.0;


  // ##### Weights and abscissas are generated by the gauleg function in help_numerics.h. gauleg then returns a 2d array, which is split below. #####
  // ##### The gauleg function is called twice, thought, because of angular integration *and* momentum integration. #####

  double** x_and_w = nullptr;
  double** x_and_w_ang = nullptr;
  double* weights_w = nullptr;
  double* absciss_x = nullptr;
  double* absciss_ang = nullptr;
  double* weights_ang = nullptr;

  double** vals = nullptr;

  std::complex<double> imag = {0.0,1.0};

  std::cout<< "imaginary component is: " << imag.imag() << "\t" << "real component is: "<< imag.real()<<std::endl;

  #ifdef MarisTandy
    std::cout<<std::endl<< "Now Using Maris Tandy Model!" << std::endl;
  #else
    std::cout<<std::endl<< "Now Using Contact Model!" << std::endl;
  #endif

  // ##### Creating abscissas and weights. Using a Log-grid. #####
  #ifdef loggrid

  std::cout<<std::endl<< "Log-Grid!" << std::endl << std::endl;
  x_and_w = gauleg(log(LAMBDA_MIN*LAMBDA_MIN),log(LAMBDA), absciss_points);

  #else

  std::cout<<std::endl<< "Linear-Grid!" << std::endl << std::endl;
  x_and_w = gauleg(0.0,LAMBDA, absciss_points);

  #endif

  // ##### Splitting 2d array up (momentum) #####
  absciss_x = x_and_w[0] + 1; // x_and_w[0];
  weights_w = x_and_w[1] + 1; // x_and_w[1];

  x_and_w_ang = gauleg(0.0,M_PI, ang_absciss_points);

  // ##### Splitting 2d array up (angular) #####

  absciss_ang = x_and_w_ang[0] + 1;
  weights_ang = x_and_w_ang[1] + 1;

  std::cout<<"Weights and abscissae calculated"<<std::endl;

  // ##### Save A and B to array "vals" #####
  vals = iterate_dressing_functions(epsilon,m_c,m_g,absciss_x,weights_w,absciss_ang,weights_ang, g_squared, eta, mu_renorm);

  // ##### Preparing for mass function M(p^2) = B(p^2)/A(p^2). #####
  double* a_vals = vals[0];
  double* b_vals = vals[1];
  double* renorm_constants = vals[2];
  double* m_vals = nullptr;
  m_vals = new double[absciss_points];
  // std::cout<<" here ok " << std::endl;

  for(int i = 0; i < absciss_points; i++){
    if(a_vals[i] == 0.0){
      m_vals[i] = 0.0;
    }
    else{
      m_vals[i] = b_vals[i]/a_vals[i];
    }
  }
  // ProgressBar pb(max_step, "Doing stuff");

  // ##### Print A (= vals[0][i]) and B (= vals[1][i]) values saved in "vals" #####
  // for(int i=0;i<absciss_points;i++)
  // {
  //   std::cout<< i << "\t"<< vals[0][i] << "\t" << vals[1][i] << std::endl;
  // }

  std::complex<double>* renormpointabvals = interpolation_cmplx(mu_renorm, m_c, renorm_constants, a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta);

  std::cout<<"A(mu) = " << renormpointabvals[0] << " B(mu) = " << renormpointabvals[1] << " M(mu) = "<< renormpointabvals[2] <<std::endl;

  std::cout<<"Z2 is: "<<  renorm_constants[0] << " Zm is: " << renorm_constants[1] << std::endl;

  std::cout<<"PionMass-Input is: "<<  PionMass << " Alpha-Angle is : " << alpha_angle << std::endl;
//
// std::cout << "# Z2 = " << renorm[0] << std::endl << "# Zm = " << renorm[1] << std::endl;
// for (int i = 0; i < INT_STEPS; ++i) {
// #ifdef LOGGRID
//     std::cout << exp(x_absci[i]) << "\t" << avals[i] << "\t" << bvals[i] << "\t" << mvals[i] << std::endl;
// #else // LOGGRID
//     std::cout << x_absci[i] << "\t" << avals[i] << "\t" << bvals[i] << "\t" << mvals[i] << std::endl;
// #endif
// }
// return 0;
// }

  // std::ofstream fileout2;
  // fileout2.open("Data/Zm(m_c)_4GeV.dat");
  // fileout2 << "# Parameters used: " <<" LAMBDA(GeV) in UV-Cuttoff in log(LAMBDA*LAMBDA): "<< LAMBDA << "LAMBDA_MIN(GeV) in IR-Cuttoff in log(LAMBDA_MIN*LAMBDA_MIN): "<< LAMBDA_MIN<< " gamma_m: "<< gamma_fun(N_C, N_F) <<std::endl;
  // fileout2 << "# mu(GeV): "<< mu_renorm << " Lambda_QCD(GeV): "<< Lambda_QCD << " Lambda_t(GeV): "<< Lambda_t << " Lambda_0 " <<Lambda_0<<std::endl;
  // fileout2 << "# q-abscissae used: "<< absciss_points <<" ang_abscissae used: "<< ang_absciss_points <<std::endl;
  // fileout2 << "# m_c " << "\t" << " zm"<<std::endl;
  // for(int i=1; i<200;i++){
  //   // ##### Save A and B to array "vals" #####
  //   vals = iterate_dressing_functions(epsilon,log(1.0+0.0001*i),m_g,absciss_x,weights_w,absciss_ang,weights_ang, g_squared, eta, mu_renorm);
  //
  //   // ##### Preparing for mass function M(p^2) = B(p^2)/A(p^2). #####
  //   double* renorm_constants = vals[2];
  //   fileout2<< log(1.0+0.0001*i) << " " << renorm_constants[1] << std::endl;
  //
  // }
  // fileout2.close();

  // ##### Save A and B values from "vals" to File #####
  ProgressBar sd(absciss_points, "Saving Data to File");
  //Total geiles neues feature
  std::ofstream  fileout;
  fileout.open("Data/DressingFunctions_A_and_B_and_M_log_600_128ang_test3.dat");
  fileout << "# Parameters used: " << "mc(GeV): "<< m_c<<" LAMBDA(GeV) in UV-Cuttoff in log(LAMBDA*LAMBDA): "<< LAMBDA << "LAMBDA_MIN(GeV) in IR-Cuttoff in log(LAMBDA_MIN*LAMBDA_MIN): "<< LAMBDA_MIN<< " gamma_m: "<< gamma_fun(N_C, N_F) <<std::endl;
  fileout << "# mu(GeV): "<< mu_renorm << " Lambda_QCD(GeV): "<< Lambda_QCD << " Lambda_t(GeV): "<< Lambda_t << " Lambda_0 " <<Lambda_0<<std::endl;
  fileout << "# q-abscissae used: "<< absciss_points <<" ang_abscissae used: "<< ang_absciss_points <<std::endl;
  fileout << "# z2 is: " << renorm_constants[0] << " zm is: " << renorm_constants[1]<<std::endl;
  fileout << "# p^2"<< " "<< "A(p^2)"<< " "<< "B(p^2)"<< " "<< "M(p^2)" << std::endl;
  for(int j=0;j<absciss_points;j++){
    ++sd;
    fileout<< exp(absciss_x[j]) << " " << vals[0][j] << " " << vals[1][j] << " " << m_vals[j] << std::endl;
  }
  fileout.close ();


  double eigenvaluebse=0.0;
  double pimass = 1.0;


  // std::complex<double>** mother3 = initialize_mother_matrix(pimass, m_c, renorm_constants,a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta, alpha_angle);

  // for(int i=1; i<10; i++){
  //   pimass = pimass + 0.0002*i;
  //   eigenvaluebse = bse_root(pimass, m_c, renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta);
  //   std::cout<< "Pionmass Input = "<< pimass << "eigenvalue = " << eigenvaluebse <<std::endl;
  // }
  double status = 0.0;
  double status2 = 0.0;
  // status = bse_root_eigenlib(0.12094, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta);

  std::complex<double>*** theta_matrix = initialize_theta_matrix(renorm_constants, absciss_x, absciss_ang,weights_w,  weights_ang, eta, alpha_angle);
  status = regulaFalsi(0.110,0.140,1e-8,m_c, renorm_constants, a_vals, b_vals, absciss_x,absciss_ang, weights_w, weights_ang, eta, theta_matrix);

  // bse_root_eigenlib(0.11982047108638232, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta);
  bse_root_eigenlib(status, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta, theta_matrix);
//   std::ofstream gmorpionmass;
//   gmorpionmass.open("Data/gmorpionmass.dat");
//   gmorpionmass<<"#Stromquarkmass and Pionmass"<<std::endl;
//   for(int i=0; i<10;i++){
//     double stromquark=0.0001+i*0.001;
//   double status=0.0;
//   status = regulaFalsi(0.001,0.5,1e-5,stromquark, renorm_constants, a_vals, b_vals, absciss_x,
//     absciss_ang, weights_w, weights_ang, eta);
//  // status = regulaFalsitest(5.0,6.0,0.0000001);
//   std::cout<< std::endl<<"Pion mass again is: " << status << std::endl<<std::endl;
  // gmorpionmass<<m_c<<" "<<status<<std::endl;
// }
// gmorpionmass.close();
 // double q, z, psi, theta, routing_plus, routing_minus;
 // std::complex<double> matrix_entry,q_plus_q_minus, q_plus_squared, q_minus_squared, p, k_squared;
 // std::complex<double> Imag = {0.0,1.0};
 //
 // routing_plus = 0.5;
 // routing_minus = routing_plus - 1.0;

 // std::complex<double>** mother_temp = nullptr;
 // mother_temp = new std::complex<double>*[absciss_points];
 // std::ofstream  fileouta;
 // fileouta.open("Data/DressingFunctions_A_and_B_PLUS_complex.dat");
 // fileouta<<"qp"<<"\t"<<"A+"<<"\t"<<"B+"<<"qm"<<"\t"<<"A-"<<"\t"<<"B-"<<std::endl;
 // for(int q_idx = 0; q_idx < absciss_points; q_idx++){
 //
 //   q = exp(0.5*absciss_x[q_idx]);
 //
 //     for(int psi_idx=0; psi_idx < ang_absciss_points; psi_idx++){
 //
 //         psi = absciss_ang[psi_idx];
 //         z = std::cos(psi);
 //
 //         q_plus_q_minus = std::pow(q,2.0) + routing_minus*q*Imag*m_pion*z + routing_plus*q*Imag*m_pion*z - routing_plus*routing_minus*m_pion*m_pion;
 //         q_plus_squared = std::pow(q,2.0) + 2.0*routing_plus*q*Imag*m_pion*z - std::pow(routing_plus,2.0)*m_pion*m_pion;
 //         q_minus_squared = std::pow(q,2.0) + 2.0*routing_minus*q*Imag*m_pion*z - std::pow(routing_minus,2.0)*m_pion*m_pion;
 //
 //         std::complex<double> abs_qp = std::sqrt(q_plus_squared);
 //         std::complex<double> abs_qm = std::sqrt(q_minus_squared);
 //
 //         // std::cout<< "sqrt(qp**2) = " << abs_qp << " sqrt(qm**2) = " << abs_qm <<std::endl;
 //
 //         std::complex<double>* plus_template = interpolation_cmplx(abs_qp, m_c, renorm_constants, a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta);
 //         std::complex<double>* minus_template = interpolation_cmplx(abs_qm, m_c, renorm_constants, a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta);
 //
 //         std::complex<double> a_plus = plus_template[0];
 //         std::complex<double> b_plus = plus_template[1];
 //         std::complex<double> a_minus = minus_template[0];
 //         std::complex<double> b_minus = minus_template[1];
 //
 //         fileouta<<q_plus_squared<<" "<< a_plus <<" "<<b_plus<<" "<<q_minus_squared<<" "<< a_minus <<" "<<b_minus<<std::endl;
 //       }
 //     }
 //     fileouta.close();


     // std::complex<double>* abtest = interpolation_cmplx(0.121, m_c, renorm_constants, a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta);
     // std::cout<< "interpolation test at x=2.133 = "<< abtest[0] <<std::endl;

  return 0;
}




// Read about BSE in the pdfs.
// Implement the BSE with Maris Tandy in the code.

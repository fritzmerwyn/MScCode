#include <iostream>
#include <math.h>
#include <omp.h>
#include "help_numerics.h"
#include "DysonSchwinger.h"
#include "progressbar.hpp"

double gamma_fun(double color = N_C, double flavor = N_F){
  return 12.0/(11.0 * color - 2.0 * flavor);
}

double running_coupling_MarisTandy(double k_squared, double eta){
  double infrared, ultraviolet;
  infrared = M_PI * ((pow(eta,7.0))/(pow(Lambda_0,4.0))) * (k_squared * k_squared) * exp( (-(eta*eta)) * (k_squared/(Lambda_0*Lambda_0)) );
  ultraviolet = (2.0 * M_PI * gamma_fun() * (1.0 - exp(-(k_squared)/(Lambda_t*Lambda_t)) ) ) / (log( M_E*M_E - 1.0 + pow(1.0 + (k_squared/(Lambda_QCD*Lambda_QCD) ),2.0) ));
  return  infrared + ultraviolet;
}

double*** initialize_matrix(double epsilon, double m_c, double* absciss_x, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu){

  double q, c2_a, c2_b, p, z, yota, k_squared, s0_a, s0_b;

  c2_a = 4.0/(2.0*3.0*pow(M_PI,2.0));
  c2_b = 4.0/(2.0*pow(M_PI,2.0));


  double*** temp_matrix = nullptr;
  temp_matrix = new double**[2];
  temp_matrix[0] = nullptr;
  temp_matrix[1] = nullptr;
  temp_matrix[0] = new double*[absciss_points];
  temp_matrix[1] = new double*[absciss_points];

#pragma omp parallel for default(none) shared(temp_matrix)
  for(int i=0;i<absciss_points;i++){
    temp_matrix[0][i] = nullptr;
    temp_matrix[1][i] = nullptr;
    temp_matrix[0][i] = new double[absciss_points + 1];
    temp_matrix[1][i] = new double[absciss_points + 1];
  }

  for(int q_idx = 0; q_idx < absciss_points; q_idx++){

    q = exp(0.5*absciss_x[q_idx]);

    for(int p_idx = 0; p_idx < absciss_points + 1; p_idx++){

      if(p_idx == absciss_points){
        p = mu;
      }

      else{
      p = exp(0.5*absciss_x[p_idx]);
      }
      s0_a=0.0;
      s0_b=0.0;

#pragma omp parallel for private(z, yota, k_squared) default(none) shared(p, q, weights_ang, absciss_ang, absciss_x, c2_a, c2_b, eta, mu) reduction(+:s0_a, s0_b)
            for(int ang_idx=0;ang_idx<ang_absciss_points;ang_idx++){ //START AT J=1 because 0th abscissa is 0. And 0 no good.(look at bottom of main)
              yota = absciss_ang[ang_idx];
              z = cos(yota);
              k_squared = p*p + q*q - 2.0*p*q*z;
              if(p==0.0){

              }
              else{
              s0_a += (c2_a * q*q*q*q)* (1.0/(p*p)) *
                      weights_ang[ang_idx] * sin(yota)*sin(yota) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z))   *
                      (running_coupling_MarisTandy(k_squared, eta) / (k_squared));

              s0_b += (c2_b * q*q*q*q) *
                      weights_ang[ang_idx] * (sin(yota)*sin(yota) *running_coupling_MarisTandy(k_squared,eta)/(k_squared));
              }
            }

      temp_matrix[0][q_idx][p_idx] = s0_a;
      temp_matrix[1][q_idx][p_idx] = s0_b;

    }
  }

  std::cout<<std::endl;
  std::cout<<"Angular Matrix initialized"<< std::endl;
  return temp_matrix;

}

double** initialize_dressing_functionAB(double a0, double b0){
  double** dress2d = nullptr;
  dress2d = new double*[2];
  dress2d[0] = nullptr;
  dress2d[1] = nullptr;
  dress2d[0] = new double[absciss_points];
  dress2d[1] = new double[absciss_points];

  ProgressBar idf(absciss_points, "Initializing Dressing Functions");

  for(int i=0;i<absciss_points;i++,++idf){
    dress2d[0][i] = a0;
    dress2d[1][i] = b0;
  }
  std::cout<<std::endl;
  std::cout<<"Dressing Functions initialized as A = "<< a0 <<" and B = "<< b0 <<std::endl;
  return dress2d;
}

double int_coupled_a(double*** angular_matrix, double* absciss_x, double* weights_w, double* a_vals, double* b_vals, int p_idx, double m_g){
  double s0_a = 0.0;
  double c1,c2; //Prefactor of Integral
  c1 = 2.0/(3.0*pow(m_g,2.0)*pow(M_PI,3.0));
  double q, z, yota, k_squared, angularpta;

#pragma omp parallel for private(q) default(none) shared(weights_w, absciss_x, a_vals, b_vals, p_idx, angular_matrix) reduction(+:s0_a)
  for(int q_idx=0;q_idx<absciss_points;q_idx++){

#ifdef MarisTandy

      #ifdef loggrid

      // c2 = g_squared/(2.0*48.0*pow(M_PI,3.0));

      q = exp(0.5*absciss_x[q_idx]);

      s0_a += (weights_w[q_idx] * angular_matrix[0][q_idx][p_idx] * a_vals[q_idx] / (pow(q*a_vals[q_idx],2.0) + pow(b_vals[q_idx],2.0)));

      #else

      c2 = g_squared/(3.0*pow(M_PI,3.0));
      q = absciss_x[q_idx];
      for(int j=0;j<absciss_points;j++){
        z = cos(absciss_ang[j]);
        yota = absciss_points[j];
        k_squared = p*p + q*q - 2.0*p*q*z;
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[q_idx] *(1.0/(p*p))*( (c2 * q*q*q * a_vals[q_idx]) / ((q*q * pow(a_vals[q_idx],2.0) + pow(b_vals[q_idx],2.0))) ) *
                weights_ang[j] * sin(yota)*sin(yota) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z)) *
                (running_coupling_MarisTandy(k_squared,eta) / (k_squared));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
        }
      }

      #endif


#else

      for(int j=0;j<absciss_points;j++){
        z = cos(absciss_ang[j]);
        yota = absciss_ang[j];
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[q_idx]*((c1 * q*q*q*q*q * a_vals[q_idx])/(p*(q*q*pow(a_vals[q_idx],2)+pow(b_vals[q_idx],2))))*
                (weights_ang[j]*sin(yota)*sin(yota)*z);
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
        }
      }
#endif
  }
  return s0_a;
}

double int_coupled_b(double*** angular_matrix, double* absciss_x, double* weights_w, double* a_vals, double* b_vals, int p_idx, double m_g){
  double s0_b = 0.0;
  double c1,c2;
  c1 = 2.0/(3.0*pow(m_g,2.0)*pow(M_PI,3.0));
  double q, z, yota, k_squared;

#pragma omp parallel for private(q) default(none) shared(weights_w, absciss_x, a_vals, b_vals, p_idx, angular_matrix) reduction(+:s0_b)
  for(int q_idx=0;q_idx<absciss_points;q_idx++){


#ifdef MarisTandy

      #ifdef loggrid

        q = exp(0.5*absciss_x[q_idx]);

        s0_b += (weights_w[q_idx] * angular_matrix[1][q_idx][p_idx] * b_vals[q_idx] / (pow(q * a_vals[q_idx],2.0) + pow(b_vals[q_idx],2.0)));

      #else

          c2 = (g_squared/pow(M_PI,3.0));

          q = absciss_x[q_idx];

          for(int j=0;j<absciss_points;j++){
            z = cos(absciss_ang[j]);
            yota = absciss_ang[j];
            k_squared = p*p + q*q - 2.0*p*q*z;
            if(p==0.0){
              return 0.0;
            }
            else{
            s0_b += (weights_w[q_idx] * (c2 * q*q*q * b_vals[q_idx] / (q*q * pow(a_vals[q_idx],2) + pow(b_vals[q_idx],2)))) *
                    (weights_ang[j] * (sin(yota)*sin(yota)*(running_coupling_MarisTandy(k_squared, eta)/(k_squared)))) ;
            // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
            // m0B = m_c + c_1*s0_b;
            // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
            // return 0;
            }
          }

      #endif


#else

    for(int j=0;j<absciss_points;j++){

      z = cos(absciss_ang[j]);
      yota = absciss_ang[j];
        // mend=m0B;
        s0_b += weights_w[q_idx]*((c1 * q*q*q*q * b_vals[q_idx])/(q*q*pow(a_vals[q_idx],2)+pow(b_vals[q_idx],2)))*
                (weights_ang[j]*sin(yota)*sin(yota));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;

      }

#endif

  }
  return s0_b;
}



double** iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu){
  // double temp_a=1e-20;
  double*** angular_matrix = initialize_matrix(epsilon,m_c,absciss_x, absciss_ang, weights_ang,g_squared,eta,mu);

  double** init= initialize_dressing_functionAB(1.0,m_c);
  double* a_vals= init[0];
  double* b_vals= init[1];
  double p;

  double new_b[absciss_points];
  double new_a[absciss_points];
  double z2, z2_old, zm, new_siga, new_sigb;
  double renorm_a[absciss_points];
  double renorm_b[absciss_points];

  // double a_end=0.1;
  double a_start=0.0;
  double b_start=0.0;
  double a_end=1.0;
  double b_end=1.0;

  double a_renorm_s=0.0;
  double b_renorm_s=0.0;
  double a_renorm_end=1.0;
  double b_renorm_end=1.0;

  zm = 1.0;
  z2 = 1.0;


  ProgressBar pb(max_iter, "Iterating Dressing Functions A and B");

  for(int j=0;j<=max_iter;j++){
      ++pb;
      if((abs((b_end-b_start)/(b_end+b_start))<epsilon) && (abs((a_end-a_start)/(a_end+a_start))<epsilon)) { //(abs((b_end-b_start)/(b_end+b_start))<epsilon)
        std::cout<<std::endl<<std::endl<< j << " iterations used. Maximum Iterations were set to "<< max_iter <<std::endl<<std::endl;
        // std::cout<<std::endl<<std::endl<< b_vals[2] << " is B(p^2)"<<std::endl<<std::endl;
        // return temp_b;
        // break;
        double** dress2d = nullptr;
        dress2d = new double*[2];
        dress2d[0] = a_vals;
        dress2d[1] = b_vals;
        return dress2d;

      }
      else{

        a_start = new_a[absciss_points-1];
        b_start = new_b[absciss_points-1];
        // ProgressBar ABupdate(absciss_points, "Calculating New A and B");

        for(int p_idx=0;p_idx<absciss_points;p_idx++){ //Integration over q
          // ++ABupdate;

          new_a[p_idx] = z2*1.0;
          new_b[p_idx] = z2*zm*m_c;

          new_a[p_idx] += z2*z2*int_coupled_a(angular_matrix, absciss_x, weights_w, a_vals, b_vals, p_idx, m_g);
          new_b[p_idx] += z2*z2*int_coupled_b(angular_matrix, absciss_x, weights_w, a_vals, b_vals, p_idx, m_g);

      }

          // new_siga = int_coupled_a(mu, m_c, m_g, absciss_x, weights_w, absciss_ang, weights_ang, a_vals, b_vals, g_squared, eta);
          new_siga = int_coupled_a(angular_matrix, absciss_x, weights_w, a_vals, b_vals, absciss_points, m_g);
          new_sigb = int_coupled_b(angular_matrix, absciss_x, weights_w, a_vals, b_vals, absciss_points, m_g);

          z2 = 1.0/(1.0 + z2*new_siga);
          zm = 1.0/z2 - z2*new_sigb/m_c;
          std::cout << std::endl <<"z2 is "<< z2 << "\t"<< "zm is " << zm <<std::endl;



        std::cout<<std::endl;
        if(j%1==0 && j !=0){
        std::cout<<std::endl<<"new_a[40] = "<< new_a[40] << "\t" << "new_b[40] = "<< new_b[40] << "\t" << "a_vals[40] = "<< a_vals[40] << "\t" << "b_vals[40] = "<< b_vals[40] <<std::endl;
        // std::cout<<"a_vals[40] = " << a_vals[40] << std::endl;
        }

        for(int k=0; k<absciss_points; k++){

          a_vals[k] = new_a[k];
          b_vals[k] = new_b[k];

        }
        // std::cout<<"a_vals[40] = " << a_vals[40] << std::endl;
      a_end = new_a[absciss_points-1];
      b_end = new_b[absciss_points-1];

    }
  }
  return 0;
}

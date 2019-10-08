#include <iostream>
#include <math.h>
#include "help_numerics.h"
#include "DysonSchwinger.h"
#include "progressbar.hpp"

double gamma_fun(double color = N_C, double flavor = N_F){
  return 12.0/(11.0 * color - 2.0 * flavor);
}

double running_coupling_MarisTandy(double k_squared, double c, double w){
  double term1, term2;
  term1 = ((c*M_PI) / pow(w,7)) * pow(k_squared/(Lambda_0*Lambda_0),2) * exp(- k_squared / (w*w*Lambda_0*Lambda_0));
  term2 = (M_PI*gamma_fun()*(1.0 - exp(-k_squared/(Lambda_0*Lambda_0)))) / (log(sqrt(M_E*M_E - 1.0 + pow(1.0 + pow(((k_squared) / (Lambda_QCD)),2),2))));
  return  term1 + term2;
}

double** initialize_dressing_functionAB(double a0, double b0){
  double** dress2d = nullptr;
  dress2d = new double*[2];
  dress2d[0] = nullptr;
  dress2d[1] = nullptr;
  dress2d[0] = new double[absciss_points];
  dress2d[1] = new double[absciss_points];

  ProgressBar idf(absciss_points, "Initializing Dressing Functions");

  for(int i=0;i<=absciss_points;i++,++idf){
    dress2d[0][i] = a0;
    dress2d[1][i] = b0;
  }
  std::cout<<std::endl;
  std::cout<<"Dressing Functions initialized as A = "<< a0 <<" and B = "<< b0 <<std::endl;
  return dress2d;
}

double int_coupled_a(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals, double g_squared, double c, double w){
  double s0_a = 0.0;
  double c1,c2; //Prefactor of Integral
  c1 = 2.0/(3.0*pow(m_g,2.0)*pow(M_PI,3.0));
  double q,z, k_squared;
  for(int i=0;i<=absciss_points;i++){




#ifdef MarisTandy

      #ifdef loggrid

      // c2 = g_squared/(2.0*48.0*pow(M_PI,3.0));
      c2 = g_squared/(2.0*3.0*pow(M_PI,3.0));

      q = exp(0.5*absciss_x[i]);
      for(int j=0;j<=absciss_points;j++){
        z = absciss_ang[j];
        k_squared = p*p + q*q - 2.0*p*q*z;
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[i] *((c2 * q*q*q*q * a_vals[i]) / (p*p*(q*q * pow(a_vals[i],2) + pow(b_vals[i],2)))) *
                weights_ang[j] * sqrt(1.0 - z*z) * (p*q*z + (2.0/(k_squared) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z))) *
                (running_coupling_MarisTandy(k_squared, c, w) / (k_squared));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
        }
      }

      #else

      c2 = g_squared/(48.0*pow(M_PI,3.0));
      // c2 = g_squared/(3.0*pow(M_PI,3.0));
      q = absciss_x[i];
      for(int j=0;j<=absciss_points;j++){
        z = absciss_ang[j];
        k_squared = p*p + q*q - 2.0*p*q*z;
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[i] *(1.0/(p*p))*( (c2 * q*q*q * a_vals[i]) / ((q*q * pow(a_vals[i],2) + pow(b_vals[i],2))) ) *
                weights_ang[j] * sqrt(1.0 - z*z) * (p*q*z + (2.0/(k_squared)) * (p*p*p*q*z - p*p*q*q - p*p*q*q*z*z + p*q*q*q*z)) *
                (running_coupling_MarisTandy(k_squared, c, w) / (k_squared));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
        }
      }

      #endif


#else

      for(int j=0;j<=absciss_points;j++){
        z = absciss_ang[j];
        if(p==0.0){
          return 0.0;
        }
        else{
        s0_a += weights_w[i]*((c1 * q*q*q*q*q * a_vals[i])/(p*(q*q*pow(a_vals[i],2)+pow(b_vals[i],2))))*
                (weights_ang[j]*sqrt(1.0-z*z)*z);
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

double int_coupled_b(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals, double g_squared, double c, double w){
  double s0_b = 0.0;
  double c1,c2;
  c1 = 2.0/(3.0*pow(m_g,2.0)*pow(M_PI,3.0));
  double q,z, k_squared;
  for(int i=0;i<=absciss_points;i++){


#ifdef MarisTandy

          #ifdef loggrid

          c2 = (g_squared/(2*pow(M_PI,3)));

          q = exp(0.5*absciss_x[i]);

          for(int j=0;j<=absciss_points;j++){
            z = absciss_ang[j];
            k_squared = p*p + q*q - 2.0*p*q*z;
            if(p==0.0){
              return 0.0;
            }
            else{
            s0_b += (weights_w[i] * (c2 * q*q*q*q * b_vals[i] / (q*q * pow(a_vals[i],2) + pow(b_vals[i],2)))) *
                    (weights_ang[j] * (sqrt(1.0 - z*z)*running_coupling_MarisTandy(k_squared, c, w)/(k_squared))) ;
            // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
            // m0B = m_c + c_1*s0_b;
            // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
            // return 0;
            }
          }

          #else

          c2 = (g_squared/pow(M_PI,3));

          q = absciss_x[i];

          for(int j=0;j<=absciss_points;j++){
            z = absciss_ang[j];
            k_squared = p*p + q*q - 2.0*p*q*z;
            if(p==0.0){
              return 0.0;
            }
            else{
            s0_b += (weights_w[i] * (c2 * q*q*q * b_vals[i] / (q*q * pow(a_vals[i],2) + pow(b_vals[i],2)))) *
                    (weights_ang[j] * (sqrt(1.0 - z*z)*(running_coupling_MarisTandy(k_squared, c, w)/(k_squared)))) ;
            // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
            // m0B = m_c + c_1*s0_b;
            // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
            // return 0;
            }
          }

          #endif


#else

    for(int j=0;j<=absciss_points;j++){

      z = absciss_ang[j];
        // mend=m0B;
        s0_b += weights_w[i]*((c1 * q*q*q*q * b_vals[i])/(q*q*pow(a_vals[i],2)+pow(b_vals[i],2)))*
                (weights_ang[j]*sqrt(1.0 - z*z));
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;

      }

#endif

  }
  return s0_b;
}

double** iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double g_squared, double c, double w){
  // double temp_a=1e-20;
  double** init= initialize_dressing_functionAB(1.0,m_c);
  double* a_vals= init[0];
  double* b_vals= init[1];
  double p;

  double new_b[absciss_points];
  double new_a[absciss_points];
  // double a_end=0.1;
  double a_start=0.0;
  double b_start=0.0;
  double a_end=1.0;
  double b_end=1.0;
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
        ProgressBar ABupdate(absciss_points, "Calculating New A and B");

        for(int i=0;i<=absciss_points;i++){
          ++ABupdate;

            #ifdef loggrid

            p = exp(0.5*absciss_x[i]);

            #else

            p = absciss_x[i];

            #endif



          new_b[i] = m_c;
          new_a[i] = 1.0;
          // a_end = temp_a;
          // temp_a = effectivemassA(epsilon, m_c,m_g,x,w,angint2,temp_b);
          // temp_b = effectivemassB(epsilon, m_c,m_g,x,w,angint,temp_a);
          new_b[i] += int_coupled_b(p, m_c, m_g, absciss_x, weights_w, absciss_ang, weights_ang, a_vals, b_vals, g_squared, c, w);
          new_a[i] += int_coupled_a(p, m_c, m_g, absciss_x, weights_w, absciss_ang, weights_ang, a_vals, b_vals, g_squared, c, w);
          // return 0;
      }
        std::cout<<std::endl;
        if(j%2==0 && j !=0){
        std::cout<<std::endl<<"new_a[40] = "<< new_a[40] << "\t" << "new_b[40] = "<< new_b[40] << "\t" << "a_vals[40] = "<< a_vals[40] << "\t" << "b_vals[40] = "<< b_vals[40] <<std::endl;
        // std::cout<<"a_vals[40] = " << a_vals[40] << std::endl;

       }
        for(int k=0;k<=absciss_points;k++){
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

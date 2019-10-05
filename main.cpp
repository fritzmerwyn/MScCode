/* This program *** calculates *** A(p^2) through DSE (iterative Eqs.) and then uses it
    to calculate the mass through DSE.
*/

#include <iostream>
#include <math.h>
#include <fstream>
#include "help_header.h"
#include <string>

#define N_rad 150
#define q_min 1e-100
#define q_max 873.0
#define max_step 150
#define max_iter 10000000.0
#define pi 3.14159265358979323846

double effectivemassA(double epsilon, double m_c, double m_g, double* x, double* w, double angint2);
double effectivemassB(double epsilon, double m_c, double m_g, double* x, double* w[], double angint);
double m0A=1e-200;
double m0B=1e-200;


double qgaus(double (*func)(double, double, double), double* x, double* w, double carry, double m0X ){
  int j;
  double s;
  s=0.0;
  for(j=1;j<=N_rad;j++){
    s += w[j]*((*func)(x[j], carry, m0X));
  }
  return s;
}

double qgaus1(double (*func)(double), double* x, double* w){
  int j;
  double s;
  s=0.0;
  for(j=1;j<=N_rad;j++){
    s += w[j]*((*func)(x[j]));
  }
  return s;
}

double massA(double x, double a, double b){
  return (exp(2.5*x)*a)/(exp(x)*a*a+b*b);
}

double massB(double x, double b, double a){
  return (exp(2*x)*b)/(exp(x)*a+b*b);
}

double angkern(double x){
  return sqrt(1-x*x);
}

double angkern2(double x){
  return sqrt(1-x*x)*x;
}

double effectivemassA(double epsilon, double m_c, double m_g, double* x, double* w, double angint2, double m0B){
  double s0_a = 0.0;
  // double m0A=1e-200;
  double mend=1.0;  //Set to 1.0, just to make the first if condition go through.
  double c_1;
  c_1 = 2.0/(3.0*pow(m_g,2.0)*pow(pi,3.0));

  for(int i=0;i<=max_iter;i++){
    // double q = absciss_x[i];
      // for(int j=0;j<=max_iter;j++){
      if(abs((mend-m0A)/(mend+m0A))<epsilon) {
        // std::cout<<std::endl<<std::endl<< i << " iterations used"<<std::endl<<std::endl;

        return m0A;
        break;
      }
      else{
        mend=m0A;
        s0_a = qgaus(massA,x,w,m0A,m0B);
        m0A = 1.0 + c_1*s0_a*angint2;
        // return 0;
      }
    }
  // }
}

double effectivemassB(double epsilon, double m_c, double m_g, double x[], double w[], double angint, double m0A){
  double s0=0.0;
  // double m0B=1e-200;
  double mend=1.0;  //Set to 1.0, just to make the first if condition go through.
  double c_1;
  c_1 = 2.0/(3.0*pow(m_g,2.0)*pow(pi,3.0));
  for(int i=0;i<=max_iter;i++){
    if(abs((mend-m0B)/(mend+m0B))<epsilon) {
      // std::cout<<std::endl<<std::endl<< i << " iterations used"<<std::endl<<std::endl;
      // std::cout<<std::endl<<std::endl<< m0A << " is A"<<std::endl<<std::endl;
      // std::cout<<std::endl<<std::endl<< m0B << " is B"<<std::endl<<std::endl;
      return m0B;
      break;
    }
    else{
      mend=m0B;
      s0 = qgaus(massB,x,w,m0B,m0A);
      m0B = m_c + c_1*s0*angint;
      // return 0;
    }
  }
}

double** initialize_dressing_functionAB(double a0, double b0){
  double** dress2d = nullptr;
  dress2d = new double*[2];
  dress2d[0] = nullptr;
  dress2d[1] = nullptr;
  dress2d[0] = new double[max_step];
  dress2d[1] = new double[max_step];

  for(int i=0;i<max_step;i++){
    dress2d[0][i] = a0;
    dress2d[1][i] = b0;
  }
  std::cout<<"Dressing Functions initialized as A = "<<a0<<" and B = "<<b0<<std::endl;
  return dress2d;
}

double int_coupled_b(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* b_vals){
  double s0_b = 0.0;
  double mend=1.0;
  double c_1;
  c_1 = 2.0/(3.0*pow(m_g,2.0)*pow(pi,3.0));
  double q,z;
  for(int i=1;i<=max_step;i++){
      q = absciss_x[i];
        mend=m0B;
        s0_b = qgaus(massB,absciss_x,weights_w,b_vals[i],1.0);
        m0B = m_c + c_1*s0_b*qgaus1(angkern,absciss_ang,weights_ang);
        // s0_b = weights_w[i]*((exp(2*q)*b_vals[i])/(exp(q)+pow(b_vals[i],2)))*qgaus1(angkern,absciss_ang,weights_ang);
        // m0B = m_c + c_1*s0_b;
        // std::cout<<std::endl<<std::endl<< j << " iterations used"<<std::endl<<std::endl;
        // return 0;
  }
}

// double dynamicmass(double epsilon, double m_c, double m_g, double x[], double w[], double angint, double angint2){
double* iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang){
  // double temp_a=1e-20;
  double** init= initialize_dressing_functionAB(1.0,m_c);
  double* b_vals= init[1];

  double new_b[max_step];
  // double a_end=0.1;
  double b_end=1.0;


  for(int i=0;i<=max_step;i++){

    new_b[i] = m_c;

    if((abs((b_end-new_b[i])/(b_end+new_b[i]))<epsilon)) {
    // std::cout<<std::endl<<std::endl<< temp_a << " is A(p^2)"<<std::endl<<std::endl;
    // std::cout<<std::endl<<std::endl<< temp_b << " is B(p^2)"<<std::endl<<std::endl;
    // return temp_b;
    // break;
      for(int k=1;k<=max_step;k++){
        b_vals[k] = new_b[k];
      }
    }
    else{
      // a_end = temp_a;
      b_end = new_b[i];
      // temp_a = effectivemassA(epsilon, m_c,m_g,x,w,angint2,temp_b);
      // temp_b = effectivemassB(epsilon, m_c,m_g,x,w,angint,temp_a);
      new_b[i] += int_coupled_b(epsilon, m_c, m_g, absciss_x, weights_w, absciss_ang, weights_ang, b_vals);

      // return 0;
    }
  }
}


int main(){
  std::cout.precision(17);
  // ##### Constants ###### //
  double m_g=132.0;
  double epsilon=10e-9;
  double epsilon2=10e-2;
  double m_c = 0.0;
  double m_c_max=10.0;
  double m_g_max= 250.0;
  double m_g_min= 120.0;
  // double h=(m_c_max-m_c)/max_step;
  double h=(m_g_max-m_g_min)/max_step;
  // double c_1;
  // c_1 = 1.0/(3.0*pow(m_g,2.0)*pow(pi,2.0));

  // double x[N_rad];
  // double w[N_rad];
  // double x1[N_rad];
  // double w1[N_rad];

  double Mend;
  double angint;
  double angint2;

  double** x_and_w = nullptr;
  double** x_and_w_ang = nullptr;
  double* absciss_x = nullptr;
  double* weights_w = nullptr;
  double* absciss_ang = nullptr;
  double* weights_ang = nullptr;

  double* vals = nullptr;


  x_and_w = gauleg(log(q_min*q_min),log(q_max*q_max), N_rad);
  absciss_x = x_and_w[0]; // x_and_w[0];
  weights_w = x_and_w[1]; // x_and_w[1];

  x_and_w_ang = gauleg(-1.0,1.0, N_rad);
  absciss_ang = x_and_w_ang[0];
  weights_ang = x_and_w_ang[1];

  // angint = qgaus1(angkern,absciss_ang,weights_ang); //Angular Integral for B(p^2)
  // angint2 = qgaus1(angkern2,x1,w1); //Angular Integral for A(p^2)

  vals = iterate_dressing_functions(epsilon, m_c, m_g, absciss_x, weights_w, absciss_ang, weights_ang);

  ProgressBar pb(max_step, "Doing stuff");

  for(int i=0;i<=max_step;i++){
    std::cout<< vals[i] << "\t" << std::endl;
  }
  std::cout<<std::endl<<"angular integral is " << angint <<" "<<std::endl;


  // Mend = effectivemassB(epsilon, m_c, m_g, x, w, angint);


  // std::ofstream  fileout;
  // fileout.open("gluonmass_dressingB78.dat");
  // for(int j=0;j<=max_step;j++,++pb){
  //   fileout<< m_g_min + j*h << " " << dynamicmass(epsilon, m_c , m_g_min + j*h , x, w, angint, angint2)<<std::endl;
  // }
  // fileout.close ();

  // std::cout<<std::endl<<"effective mass is " << temp_b <<" MeV"<<std::endl;
}

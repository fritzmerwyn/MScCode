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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


double regulaFalsi(double low_mass, double high_mass, double epsilon, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x,
  double* absciss_ang, double* weights_w, double* weights_ang, double eta)
{
      double f_low = bse_root(low_mass, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta);
      double f_high = bse_root(high_mass, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta);

      double x_low, x_high, delta_x, delta_mass, c, possibleroot;

    if ( f_low * f_high >= 0.0)
    {
        std::cout << "Variables low_mass and high_mass not assumed correctly. Try again\n";
        return 0;
    }

    if(f_low < 0.0){
      x_low = low_mass;
      x_high = high_mass;
    }
    else{
      x_low = high_mass;
      x_high = low_mass;
    }

    delta_x = high_mass - low_mass;

    for (int i=0; i < max_iter; i++)
    {
        // Find the point that touches x axis
        c = ((x_low*f_high - x_high*f_low )/ (f_high - f_low));
        possibleroot = bse_root(c, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x, absciss_ang,  weights_w,  weights_ang,  eta);

        // Check if the above found point is root
        if(possibleroot*f_high < 0.0){
          delta_mass = x_low - c;
          x_low = c;
          f_low = possibleroot;
        }
        else{
          delta_mass = x_high - c;
          x_high = c;
          f_high = possibleroot;
        }
        delta_x = x_high - x_low;

        std::cout<<std::endl<<possibleroot<<std::endl;
        // Decide the side to repeat the steps
        if(abs(possibleroot) < epsilon){
          std::cout<< "Difference: "<< abs(possibleroot) - epsilon << std::endl;

          std::cout << "The value of root is : " << c << std::endl;
          return c;
          break;
        }
    }
    std::cout<< "Reached maximum Iterations for regulaFalsitest! " <<std::endl;
    return 0;
}


double regulaFalsitest(double low_mass, double high_mass, double epsilon)
{
      double f_low = low_mass*low_mass-4.0*low_mass-10.0;
      double f_high = high_mass*high_mass-4.0*high_mass-10.0;
      double x_low, x_high, delta_x, delta_mass, c, possibleroot;

    if ( f_low * f_high >= 0.0)
    {
        std::cout << "Variables low_mass and high_mass not assumed correctly. Try again\n";
        return 0;
    }

    if(f_low < 0.0){
      x_low = low_mass;
      x_high = high_mass;
    }
    else{
      x_low = high_mass;
      x_high = low_mass;
    }

    delta_x = high_mass - low_mass;

    for (int i=0; i < max_iter; i++)
    {
        // Find the point that touches x axis
        c = ((x_low*f_high - x_high*f_low )/ (f_high - f_low));
        possibleroot = c*c-4.0*c-10.0;

        // Check if the above found point is root
        if(possibleroot*f_high < 0.0){
          delta_mass = x_low - c;
          x_low = c;
          f_low = possibleroot;
        }
        else{
          delta_mass = x_high - c;
          x_high = c;
          f_high = possibleroot;
        }
        delta_x = x_high - x_low;

        std::cout<<std::endl<<possibleroot<<std::endl;
        // Decide the side to repeat the steps
        if(abs(possibleroot) < epsilon){

          std::cout << "The value of root is : " << c << std::endl;
          return c;
          break;
        }
    }
    std::cout<< "Reached maximum Iterations for regulaFalsitest! " <<std::endl;
    return 0;
}


struct quadratic_params
  {
    double m_c,eta;
    double* renorm_constants,a_vals,b_vals,absciss_x,absciss_ang,weights_w,weights_ang;
  };

// double bsefunction(double mass, void* params){
//   struct quadratic_params *p = (struct quadratic_params *) params;
//
//   double m_c = p-> m_c;
//   double eta = p-> eta;
//   double* renorm_constants = p-> renorm_constants;
//   double* a_vals = p-> a_vals;
//   double* b_vals = p-> b_vals;
//   double* absciss_x = p-> absciss_x;
//   double* absciss_ang = p-> absciss_ang;
//   double* weights_w = p-> weights_w;
//   double* weights_ang = p-> weights_ang;
//
//   return bse_root(mass, m_c,  renorm_constants,  a_vals,  b_vals,  absciss_x,
//      absciss_ang,  weights_w,  weights_ang,  eta);
// }
//
// double findRoot(double low_mass, double high_mass, double m_c, double eta, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x,
//   double* absciss_ang, double* weights_w, double* weights_ang){
//
// int status;
//   int iter = 0;
//   const gsl_root_fsolver_type *T;
//   gsl_root_fsolver *s;
//   double r = 0, r_expected = sqrt (5.0);
//   double x_lo = low_mass, x_hi = high_mass;
//   gsl_function F;
//   struct quadratic_params params = {m_c,eta,renorm_constants,a_vals,b_vals,absciss_x,absciss_ang,weights_w, weights_ang};
//
//   F.function = &bsefunction;
//   F.params = &params;
//
//   T = gsl_root_fsolver_brent;
//   s = gsl_root_fsolver_alloc (T);
//   gsl_root_fsolver_set (s, &F, x_lo, x_hi);
//
//   printf ("using %s method\n",
//           gsl_root_fsolver_name (s));
//
//   printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//           "iter", "lower", "upper", "root",
//           "err", "err(est)");
//
//   do
//     {
//       iter++;
//       status = gsl_root_fsolver_iterate (s);
//       r = gsl_root_fsolver_root (s);
//       x_lo = gsl_root_fsolver_x_lower (s);
//       x_hi = gsl_root_fsolver_x_upper (s);
//       status = gsl_root_test_interval (x_lo, x_hi,
//                                        0, 0.001);
//
//       if (status == GSL_SUCCESS)
//         printf ("Converged:\n");
//
//       printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//               iter, x_lo, x_hi,
//               r, r - r_expected,
//               x_hi - x_lo);
//     }
//   while (status == GSL_CONTINUE && iter < max_iter);
//
//   gsl_root_fsolver_free (s);
//
//   return status;
// }

double bse_root(double pionmass, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x,
  double* absciss_ang, double* weights_w, double* weights_ang, double eta){

std::cout<<std::endl<< "Generating Mother Matrix" << std::endl;
std::complex<double>** mother = initialize_mother_matrix(pionmass, m_c, renorm_constants,a_vals, b_vals, absciss_x, absciss_ang, weights_w, weights_ang, eta, alpha_angle);

// double* gsl_mother_temp = nullptr;
// gsl_mother_temp = new double[absciss_points*absciss_points];

// for(int i = 0; i< absciss_points; i ++){
//   for(int j = 0; j< absciss_points; j++){
//     std::cout<< mother[i][j].real() << " " << mother[i][j].imag() << std::endl;
//   }
// }



// ##### WRITE FILE ######
// std::ofstream fileout3;
// fileout3.open("Data/MotherMatrix_temp.dat");
// for(int i = 0; i< absciss_points; i ++){
//   for(int j = 0; j< absciss_points; j++){
//     fileout3 << mother[i][j].real() << std::endl;
//   }
// }
// fileout3.close();
//
// // ##### READ FILE #####
// std::ifstream input("Data/MotherMatrix_temp.dat");
//  if (!input) {
//    std::cout << "Cannot open file.\n";
//    return 0;
//  }
//  for (int i = 0; i < absciss_points*absciss_points; i++) {
//    input >> gsl_mother_temp[i];
//  }
//  input.close();
// #####  #####

// ##### SOLVE EIGENVALUE PROBLEM WITH GSL LIBRARY #####
// gsl_matrix_view gsl_mother
//   = gsl_matrix_view_array (gsl_mother_temp, absciss_points, absciss_points);

gsl_vector_complex *eval = gsl_vector_complex_alloc (absciss_points);
gsl_matrix_complex *evec = gsl_matrix_complex_alloc (absciss_points, absciss_points);

gsl_eigen_nonsymmv_workspace * mother_workspace =
  gsl_eigen_nonsymmv_alloc (absciss_points);

gsl_matrix* gsl_mother_temp = gsl_matrix_alloc(absciss_points,absciss_points);


for(int i=0; i<absciss_points; i++){
  for(int j=0; j<absciss_points; j++){
    gsl_matrix_set(gsl_mother_temp,i,j,mother[i][j].real());
  }
}

// gsl_eigen_nonsymmv_params(0, 1, mother_workspace);
gsl_eigen_nonsymmv (gsl_mother_temp, eval, evec, mother_workspace);

gsl_eigen_nonsymmv_free (mother_workspace);

gsl_eigen_nonsymmv_sort (eval,evec,
                         GSL_EIGEN_SORT_ABS_DESC);

std::cout<< "first eigenvalue of mother-matrix is: "<< GSL_REAL(gsl_vector_complex_get(eval,0))<<std::endl;

double root = GSL_REAL(gsl_vector_complex_get(eval,0));

gsl_vector_complex_view evec_0
           = gsl_matrix_complex_column (evec, 0);

for (int j = 0; j < absciss_points; ++j)
          {
            gsl_complex z =
              gsl_vector_complex_get(&evec_0.vector, j);
            std::cout<<GSL_REAL(z)<<" + i* "<< GSL_IMAG(z)<<std::endl;
          }

return root - 1.0;

// std::ofstream fileout3;
// fileout3.open("Data/VectorMatrixPion_temp_1st.dat");

// fileout3.close();

gsl_vector_complex_free(eval);
gsl_matrix_complex_free(evec);
gsl_matrix_free(gsl_mother_temp);

// return 0;
}

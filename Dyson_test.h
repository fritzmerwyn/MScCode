#include "definitions.h"
#include <complex>


double** initialize_dressing_functionAB(double a0, double b0);

double*** initialize_matrix(double epsilon, double m_c, double* absciss_x, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu);

double** initialize_dressing_functionAB(double a0, double b0);

std::complex <double>** initialize_mother_matrix(double m_pion, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x, double* absciss_ang, double* weights_w, double* weights_ang, double eta, double alpha);

double int_coupled_a(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals, double g_squared, double eta);

double int_coupled_b(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals, double g_squared, double eta);

double** iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu);

double gamma_fun(double NC, double NF);

std::complex<double>* interpolation_cmplx(std::complex<double> p, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x, double* absciss_ang, double* weights_w, double* weights_ang, double eta);

double regulaFalsi(double low_mass, double high_mass, double epsilon, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x,
  double* absciss_ang, double* weights_w, double* weights_ang, double eta);

double* bse_root(double pionmass, double m_c, double* renorm_constants, double* a_vals, double* b_vals, double* absciss_x,
  double* absciss_ang, double* weights_w, double* weights_ang, double eta);

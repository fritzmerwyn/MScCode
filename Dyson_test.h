#include "definitions.h"



double** initialize_dressing_functionAB(double a0, double b0);

double int_coupled_a(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals, double g_squared, double eta);

double int_coupled_b(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals, double g_squared, double eta);

double** iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double g_squared, double eta, double mu);

double gamma_fun(double NC, double NF);

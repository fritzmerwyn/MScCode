#include "definitions.h"



double** initialize_dressing_functionAB(double a0, double b0);

double int_coupled_a(double p, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals);

double int_coupled_b(double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang, double* a_vals, double* b_vals);

double** iterate_dressing_functions(double epsilon, double m_c, double m_g, double* absciss_x, double* weights_w, double* absciss_ang, double* weights_ang);

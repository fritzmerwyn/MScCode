#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include "help_numerics.h"
#include "DysonSchwinger.h"
// #include "progressbar.hpp"

std::complex<double> funccomplex(std::complex<double> z);

int precalculation(double*absciss_x, double*absciss_ang, double*weights_w, double*weights_ang);

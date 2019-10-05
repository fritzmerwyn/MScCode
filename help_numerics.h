#ifndef help_H    // So you dont declare header/function twice
#define help_H

// Gauss Legendre Code from Numerical Recipes in C
// and more numerical functions for the integration.
#include "definitions.h"
#include <iostream>
#include <math.h>
#define EPS 3.0e-11 // EPS is the relative precision.

double** gauleg(double a,double b,int e);
double angkern(double x);
double angkern2(double x);
double qgaus1(double (*func)(double), double* x, double* w);

#endif

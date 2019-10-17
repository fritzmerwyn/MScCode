#include <math.h>
#include <iostream>
#include "help_numerics.h"
#define EPS 3.0e-11 // EPS is the relative precision.

double qgaus1(double (*func)(double), double* x, double* w){
  int j;
  double s;
  s=0.0;
  for(j=0;j<=absciss_points;j++){
    s += w[j]*((*func)(x[j]));
  }
  return s;
}

double angkern(double x){
  return sqrt(1-x*x);
}

double angkern2(double x){
  return sqrt(1-x*x)*x;
}

double** gauleg(double x1, double x2, int n)
// Given the lower and upper limits of integration x1 and x2, and given n,
// this routine returns arrays x[1..n] and w[1..n] of length n, containing
// the abscissas and weights of the Gauss- Legendre n-point quadrature formula.
{
int m,j,i;
double z1,z,xm,xl,pp,p3,p2,p1;
double* x = nullptr;
double* w = nullptr;

x = new double[n];
w = new double[n];

m=(n+1)/2;
xm=0.5*(x2+x1);
xl=0.5*(x2-x1);
for(i=1;i<=m;i++) {
  z=cos(M_PI*(i-0.25)/(n+0.5));
  // Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.
  do {
    p1=1.0;
    p2=0.0;
    for (j=1;j<=n;j++) {
      p3=p2;
      p2=p1;
      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
    };
    // p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2, the polynomial of one lower order.
    pp=n*(z*p1-p2)/(z*z-1.0);
    z1=z;
    z=z1-p1/pp;
  } while (fabs(z-z1) > EPS);
  x[i]=xm-xl*z;
  x[n+1-i]=xm+xl*z;
  w[i]=2.0*xl/((1.0-z*z)*pp*pp);
  w[n+1-i]=w[i];
// Scale the root to the desired interval, and put in its symmetric counterpart. Compute the weight
// and its symmetric counterpart.
  }
  // Initialises a pointer to a 2d array (abcissas_gl2d) which contains abscissas x[j] and weights
  // w[j] in [0] and [1] respectively.
  // This array is then returned.
  double** abcissas_gl2d = nullptr;
  abcissas_gl2d = new double*[2];
  abcissas_gl2d[0] = x;
  abcissas_gl2d[1] = w;

  return abcissas_gl2d;
}

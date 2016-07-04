/**

   \file sphericalbessel.c

   gives Spherical Bessel Functions j_l(x) and i_l(x) = i^(-l) j_l(i x)
   uses lowest order Taylor expansion for small arguments
   use recursion relations for list of bessel functions

   simple and fast, only work for reasonably small l

   alternative implementation using coulcc: higher precision but 
   almost 10 times slower
   
   (c) 2012 Thomas Neff
*/

#include <math.h>
#include <complex.h>

#include "sphericalbessel.h"
#include "coulcc.h"


#define SQR(x)	((x)*(x))


double besselj0(double x)
{
  const double EPSILON = 1E-8;

  if (fabs(x) < EPSILON)
    return (1.0 - SQR(x)/6.0);
  else
    return sin(x)/x;
}


complex double cbesselj0(complex double x)
{
  const double EPSILON = 1E-8;

  if (cabs(x) < EPSILON)
    return (1.0 - SQR(x)/6.0);
  else
    return csin(x)/x;
}


// recursion upwards becomes numerically instable for |x| < l
// for |x| > 0.5 they are ok for not too large l
// for smaller x use taylor expansion 

void besseljs(int lmax, double x, double* jl)
{
  const double EPSILON = 0.5;
  double x2 = SQR(x);

  // low order Taylor for small x
  if (fabs(x) < EPSILON) {
    jl[0] = 1.0 - x2/6.0* (1.0+x2/20.0);
    for (int l=1; l<=lmax; l++)
      jl[l] = jl[l-1]*x/(2*l+1) *(1.0+x2/((2l+1)*(2l+3))* (1.0+x2/((2l+1)*(2l+5))));
  } else {
    jl[0] = sin(x)/x;
    if (lmax==0) return;
    jl[1] = (sin(x)-x*cos(x))/SQR(x);
    for (int l=2; l<=lmax; l++)
      jl[l] = (2*l-1)*jl[l-1]/x-jl[l-2];
  }
}


void cbesseljs(int lmax, complex double x, complex double* jl)
{
  const double EPSILON = 0.5;
  complex double x2 = SQR(x);

  // low order Taylor for small x
  if (cabs(x) < EPSILON) {
    jl[0] = 1.0 - x2/6.0* (1.0+x2/20.0);
    for (int l=1; l<=lmax; l++)
      jl[l] = jl[l-1]*x/(2*l+1) *(1.0+x2/((2l+1)*(2l+3))* (1.0+x2/((2l+1)*(2l+5))));
  } else {
    jl[0] = csin(x)/x;
    if (lmax==0) return;
    jl[1] = (csin(x)-x*ccos(x))/SQR(x);
    for (int l=2; l<=lmax; l++)
      jl[l] = (2*l-1)*jl[l-1]/x-jl[l-2];
  }
}


/*
// using coulcc
void cbesseljs(int lmax, complex double x, complex double* jl)
{
  complex double XX = x;
  complex double ETA1 = 0.0;
  complex double ZLMIN = 0;
  int NL = lmax+1;
  complex double FC[NL];
  complex double GC[1];
  complex double FCP[1];
  complex double GCP[1];
  complex double SIG[1];
  int MODE1 = 4; // possible scaling with -4
  int KFN = 1; 
  int IFAIL = 1; // print error messages

  FORTRAN(coulcc)(&XX, &ETA1, &ZLMIN, &NL, FC, GC, FCP, GCP, SIG, 
		  &MODE1, &KFN, &IFAIL); 
  
  int l;
  for (l=0; l<NL; l++)
    jl[l] = FC[l];
}
*/
    

double besseli0(double x)
{
  const double EPSILON = 1E-8;

  if (fabs(x) < EPSILON)
    return(1.0 + SQR(x)/6.0);
  else
    return sinh(x)/x;
}


complex double cbesseli0(complex double x)
{
  const double EPSILON = 1E-8;

  if (cabs(x) < EPSILON)
    return(1.0 + SQR(x)/6.0);
  else
    return csinh(x)/x;
}



// recursion upwards becomes numerically instable for |x| < l
// for |x| > 0.5 they are ok for not too large l
// for smaller x use taylor expansion 

void besselis(int lmax, double x, double* il)
{
  const double EPSILON = 0.5;
  double x2 = SQR(x);

  // Taylor expansion for small x
  if (fabs(x) < EPSILON) {
    il[0] = 1.0 + x2/6.0* (1.0+x2/20.0);
    for (int l=1; l<=lmax; l++)
      il[l] = il[l-1]*x/(2*l+1)*(1.0-x2/((2l+1)*(2l+3))*(1.0-x2/((2l+1)*(2l+5))));
  } else {
    il[0] = sinh(x)/x;
    if (lmax==0) return;
    il[1] = -(sinh(x)-x*cosh(x))/x2;
    for (int l=2; l<=lmax; l++)
      il[l] = -(2*l-1)*il[l-1]/x+il[l-2];
  }
}


void cbesselis(int lmax, complex double x, complex double* il)
{
  const double EPSILON = 0.5;
  complex double x2 = SQR(x);

  // Taylor expansion for small x
  if (cabs(x) < EPSILON) {
    il[0] = 1.0 + x2/6.0* (1.0+x2/20.0);
    for (int l=1; l<=lmax; l++)
      il[l] = il[l-1]*x/(2*l+1)*(1.0-x2/((2l+1)*(2l+3))*(1.0-x2/((2l+1)*(2l+5))));
  } else {
    il[0] = csinh(x)/x;
    if (lmax==0) return;
    il[1] = -(csinh(x)-x*ccosh(x))/x2;
    for (int l=2; l<=lmax; l++)
      il[l] = -(2*l-1)*il[l-1]/x+il[l-2];
  }
}

/*
// using coulcc
void cbesselis(int lmax, complex double x, complex double* il)
{
  // coulcc does not work for x == 0
  // use Taylor expansion for very small x

  double EPSILON=1E-8;

  if (cabs(x) < EPSILON) {

    il[0] = 1.0 + SQR(x)/6.0;
    for (int l=1; l<=lmax; l++)
      il[l] = il[l-1]*x/(2*l+1)*(1.0-SQR(x)/((2l+1)*(2l+3)));

  } else {

    complex double XX = I*x;
    complex double ETA1 = 0.0;
    complex double ZLMIN = 0;
    int NL = lmax+1;
    complex double FC[NL];
    complex double GC[1];
    complex double FCP[1];
    complex double GCP[1];
    complex double SIG[1];
    int MODE1 = 4; // -4: perfrom scaling with exponential factor
    int KFN = 1; 
    int IFAIL = 1; // 1: print error messages

    FORTRAN(coulcc)(&XX, &ETA1, &ZLMIN, &NL, FC, GC, FCP, GCP, SIG, 
		    &MODE1, &KFN, &IFAIL); 

    double complex itoml[] = {1, -I, -1, I}; 
    int l;

    for (l=0; l<NL; l++)
      il[l] = itoml[l%4]*FC[l];
  }
}
*/

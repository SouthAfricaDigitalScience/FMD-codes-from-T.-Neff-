/**

  \file wignerd.h

  calculates Wigner D-function

  (c) 2003 Thomas Neff

*/


#ifndef _WIGNERD_H
#define _WIGNERD_H

#include <complex.h>

#include "fortranc.h"


double FORTRAN(djmnb)(const int* j, const int* m, const int* n, 
		      const double* b);


inline static double djmk(int j, int m, int k, double beta)
{
  return (FORTRAN(djmnb)(&j, &m, &k, &beta));
}


inline static complex double Djmk(int j, int m, int k, 
				  double alpha, double beta, double gamma)
{
  return (cexp(-0.5*I*(m*alpha+k*gamma))*FORTRAN(djmnb)(&j, &m, &k, &beta));
}


/// j, m, k are twice the physical quantities: spin 1/2 -> j=1
inline static complex double Djmkstar(int j, int m, int k, 
				      double alpha, double beta, double gamma)
{
  return (cexp(0.5*I*(m*alpha+k*gamma))*FORTRAN(djmnb)(&j, &m, &k, &beta));
}

#endif

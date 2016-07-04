/**

   \file interpol.h

   interpolate data by cubic spline
   uses routine from toms 547

   2006 Thomas Neff
*/


#ifndef _INTERPOL_H
#define _INTERPOL_H


typedef struct {
  int n;
  double h;
  double* T;
  double* G;
  double* B;
  double* C;
  double* D;
} interpolpara;


void initInterpolzeroderivative(interpolpara* ipp, int n, double h, 
				double* x, double* f);

void freeInterpol(interpolpara* ipp);


double calcInterpol(const interpolpara* ipp, double x);


#endif

/**

   \file interpol.c

   interpolate data by cubic spline
   uses routine from toms 547

   2006 Thomas Neff
*/

#include <stdio.h>
#include <stdlib.h>

#include "fortranc.h"
#include "interpol.h"


void FORTRAN(dcsint)(int* IENT, double* H, int* N, double* TNODE, double* G,
		     double* END1, double* ENDN,
		     double* B, double* C, double* D);



void initInterpolzeroderivative(interpolpara* ipp, int n, double h, 
				double* x, double* f)
{
  ipp->n = n;
  ipp->h = h;
  
  ipp->T = malloc(n*sizeof(double));
  ipp->G = malloc(n*sizeof(double));
  ipp->B = malloc(n*sizeof(double));
  ipp->C = malloc(n*sizeof(double));
  ipp->D = malloc(n*sizeof(double));

  int i;
  for (i=0; i<n; i++)
    ipp->T[i] = x[i];
  for (i=0; i<n; i++)
    ipp->G[i] = f[i];

  int IENT=1;
  double END1=0.0; double ENDN=0.0;
  FORTRAN(dcsint)(&IENT, &h, &n, ipp->T, ipp->G, &END1, &ENDN, 
		  ipp->B, ipp->C, ipp->D);
}


void freeInterpol(interpolpara* ipp)
{
  free(ipp->T);
  free(ipp->G);
  free(ipp->B);
  free(ipp->C);
  free(ipp->D);
}


inline double sqr(double x) { return x*x; }
inline double cbe(double x) { return x*x*x; }

double calcInterpol(const interpolpara* p, double x)
{
  int i=0;

  if (x < p->T[0] || x > p->T[p->n-1]) {
    fprintf(stderr, "interpol: argument %f outside interpolation range\n", x);
    return 0.0;
  }

  while (i < p->n - 1) {
    if (x <= p->T[i+1])
      break;
    i++;
  }

  return (p->G[i] + 
	  p->B[i]*(x - p->T[i]) +
	  p->C[i]*sqr(x - p->T[i]) +
	  p->D[i]*cbe(x - p->T[i]));
}

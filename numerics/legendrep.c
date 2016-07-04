/**

   \file legendrep.c

   calculates Legendre polynomials
   
   (c) 2012 Thomas Neff
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "legendrep.h"


#define SQR(x)	((x)*(x))


double legendrep(int l, double x)
{
  switch (l) {
  case 0:
    return 1.0;
  case 1:
    return x;
  case 2:
    return (1.5*SQR(x)-0.5);
  case 3:
    return (2.5*SQR(x)-1.5)*x;
  case 4:
    return ((4.375*SQR(x)-3.75)*SQR(x)+0.375);
  default:
    fprintf(stderr, "legendrep.c: legendrep(l,x) not implemented for l>4\n");
    exit(-1);
  }
}


// use recursion relation
void legendreps(int lmax, double x, double* lp)
{
  lp[0] = 1.0;
  if (lmax==0) return;
  
  lp[1] = x;
  for (int l=2; l<=lmax; l++) 
    lp[l] = ((2*l-1)*x*lp[l-1]-(l-1)*lp[l-2])/l;
}


complex double clegendrep(int l, complex double x)
{
  switch (l) {
  case 0:
    return 1.0;
  case 1:
    return x;
  case 2:
    return (1.5*SQR(x)-0.5);
  case 3:
    return (2.5*SQR(x)-1.5)*x;
  case 4:
    return ((4.375*SQR(x)-3.75)*SQR(x)+0.375);
  default:
    fprintf(stderr, "legendrep.c: legendrep(l,x) not implemented for l>4\n");
    exit(-1);
  }
}


// use recursion relation
void clegendreps(int lmax, complex double x, complex double* lp)
{
  lp[0] = 1.0;
  if (lmax==0) return;
  
  lp[1] = x;
  for (int l=2; l<=lmax; l++) 
    lp[l] = ((2*l-1)*x*lp[l-1]-(l-1)*lp[l-2])/l;
}

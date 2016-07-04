/**

   \file legendrep.h

   calculates Legendre polynomials

   (c) 2012 Thomas Neff
*/

#ifndef _LEGENDREP_H
#define _LEGENDREP_H

#include <complex.h>

// for a given l
double legendrep(int l, double x);

complex double clegendrep(int l, complex double x);


// for all l up to lmax, result in lp[]
void legendreps(int lmax, double x, double* lp);

void clegendreps(int lmax, complex double x, complex double* lp);


#endif

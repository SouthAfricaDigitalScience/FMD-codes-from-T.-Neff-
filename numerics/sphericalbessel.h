/**

   \file sphericalbessel.h

   Spherical Bessel Functions j_l(x) and i_l(x) = i^(-l) j_l(i x)

   (c) 2012 Thomas Neff
*/

#ifndef _SPHERICALBESSEL_H
#define _SPHERICALBESSEL_H

#include <complex.h>

// l=0
double besselj0(double x);
complex double cbesselj0(complex double x);

// l=0 up to l=lmax, results in jl[]
void besseljs(int lmax, double x, double* jl);
void cbesseljs(int lmax, complex double x, complex double* jl);

// l=0
double besseli0(double x);
complex double cbesseli0(complex double x);

// l=0 up to l=lmax, results in il[]
void besselis(int lmax, double x, double* il);
void cbesselis(int lmax, complex double x, complex double* il);


#endif

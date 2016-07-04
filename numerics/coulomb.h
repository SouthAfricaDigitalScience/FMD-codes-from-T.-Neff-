/**

  \file coulomb.h

  calculates 1/sqrt(z) erf(sqrt(z)) for complex z
  using a pade approximation


  (c) 2003 Thomas Neff

*/


#ifndef _COULOMB_H
#define _COULOMB_H


#include <complex.h>


complex double zcoulomb(complex double z);

complex double dzcoulomb(complex double z);


#endif


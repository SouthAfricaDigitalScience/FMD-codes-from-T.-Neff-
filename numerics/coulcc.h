/**

   \file coulcc.h

   Coulomb and Bessel functions
   by A. R. Barnett and I. J. Thompson

   (c) 2012 Thomas Neff

*/

#ifndef _COULCC_H
#define _COULCC_H


#include "fortranc.h"

void FORTRAN(coulcc)(complex double* xx, complex double* eta1,
                     complex double* zlmin, int* nl,
                     complex double* fc, complex double* gc,
                     complex double* fcp, complex double* gcp,
                     complex double* sig,
                     int* mode1, int* kfn, int* fail);

#endif

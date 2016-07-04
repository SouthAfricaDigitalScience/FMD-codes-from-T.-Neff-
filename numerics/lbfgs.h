/**

  \file lbfgs.h

  Limited memory BFGS method for large scale optimization
  by Jorge Nocedal

  (c) 2003 Thomas Neff

*/


#ifndef _LBFGS_H
#define _LBFGS_H


#include "fortranc.h"


void FORTRAN(lbfgs)(const int* n, const int* m, double* x, const double* f,
		    const double* gradf, unsigned int* diagco, double* diag,
		    const int* iprint, const double* eps, const double* xtol,
		    double* w, int* iflag);

#endif

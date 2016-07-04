/**

   \file DensityMatrixHO.h

   calculate density matrix in harmonic oscillator basis


   (c) 2006 Thomas Neff

*/


#ifndef _DENSITYMATRIXHO_H
#define _DENSITYMATRIXHO_H

#include <complex.h>

#include "SlaterDet.h"
#include "HOBasis.h"


typedef struct {
  int nmax;
  int dim;
  double omega;
} DensityMatrixHOPar;


void calcDensityMatrixHO(const DensityMatrixHOPar* par, 
			 const SlaterDet* Q, const SlaterDetAux* X,
			 complex double* rho);

void calcDensityMatrixHOod(const DensityMatrixHOPar* par,
			   const SlaterDet* Q, const SlaterDet* Qp,
			   const SlaterDetAux* X,
			   complex double* rho);

#endif

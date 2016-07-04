/**

   \file ProjectedDensityMatrixHO.h

   calculate the projected density matrix in
   harmonic oscillator basis


   (c) 2006 Thomas Neff

*/


#ifndef _PROJECTEDDENSITYMATRIXHO_H
#define _PROJECTEDDENSITYMATRIXHO_H

#include <complex.h>

#include "SlaterDet.h"
#include "Projection.h"
#include "DensityMatrixHO.h"



extern ManyBodyOperator OpDiagonalDensityMatrixHO;


void initOpDiagonalDensityMatrixHO(DensityMatrixHOPar* par);

void calcDiagonalDensityMatrixHOod(DensityMatrixHOPar* par,
				   const SlaterDet* Q, const SlaterDet* Qp,
				   const SlaterDetAux* X,
				   complex double* rhodiag);

void writeprojectedDiagonalDensityMatrixHO(FILE* fp,
					   const Projection* P,
					   const DensityMatrixHOPar* dmpar,
					   int j, int pi, int a,
					   void* dmatrix,
					   const Eigenstates* E);

void writeprojectedDiagonalDensityMatricesHO(FILE* fp,
					     const Projection* P,
					     const DensityMatrixHOPar* dmpar,
					     void* dmatrix,
					     const Eigenstates* E);

void writeprojectedDiagonalTransitionDensityMatrixHO(FILE* fp,
						     const Projection* P,
						     const DensityMatrixHOPar* dmpar,
						     int j, int pi,
						     int beta, int alpha,
						     void* dmatrix,
						     const Eigenstates* E);
 
#endif

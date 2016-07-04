/**

  \file gradSlaterDet.h

  gradients for SlaterDet matrixelements


  (c) 2003 Thomas Neff

*/


#ifndef _GRADSLATERDET_H
#define _GRADSLATERDET_H

#include <stdio.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"


/// Auxiliaries for gradients of SlaterDeterminants
typedef struct {
  int ngauss;
  gradGaussianAux* dGaux;	///< Gaussian gradients dGaux(ngauss,ngauss)
  gradGaussian* dno;		///< overlap gradients dno(ngauss,A) 
} gradSlaterDetAux;


/// gradient of functionial f from Slater determinant
typedef struct {
  int ngauss;
  complex double val;
  gradGaussian* gradval;
} gradSlaterDet;


void allocategradSlaterDet(gradSlaterDet* G, int A);

/// init Slater determinant gradient
void initgradSlaterDet(const SlaterDet* Q, gradSlaterDet* G);

void freegradSlaterDet(gradSlaterDet* G);

void zerogradSlaterDet(gradSlaterDet* G);

/// 
void copygradSlaterDet(const gradSlaterDet* dG, gradSlaterDet* dGp);

void multgradSlaterDet(gradSlaterDet* G, complex double x);

void addtogradSlaterDet(gradSlaterDet* G, const gradSlaterDet* dG);

void addmulttogradSlaterDet(gradSlaterDet* G, const gradSlaterDet* dG, 
			    complex double x);

/// rotate gradient by Euler angles
void rotategradSlaterDet(gradSlaterDet* G, double alpha, double beta, double gamma);

/// allocate space for gradSlaterDetAux
void allocategradSlaterDetAux(gradSlaterDetAux* dX, int A);

/// init gradSlaterDetAux dX
void initgradSlaterDetAux(const SlaterDet* Q, gradSlaterDetAux* dX);

/// free memory used by gradSlaterDetAux dX
void freegradSlaterDetAux(gradSlaterDetAux* dX);

/// calculate SlaterDetAux for SlaterDet Q.
/// hermiticity not exploited yet
void calcgradSlaterDetAux(const SlaterDet* Q, const SlaterDetAux* X,
			  gradSlaterDetAux* dX);

/// calculate  SlaterDetAux for SlaterDet's Q and Qp
void calcgradSlaterDetAuxod(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X,
			    gradSlaterDetAux* dX);

/// calculate gradient of Ovlap <Q|Qp>
void calcgradOvlapod(const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X, const gradSlaterDetAux* dX,
		     gradSlaterDet* G);

/// calculate gradient of matrix element with SaterDet Q for 
/// one-body operator op defined for Gaussians
/// val and dval(ngauss) will be added up
void calcgradSlaterDetOBME(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   const gradOneBodyOperator* op, 
			   gradSlaterDet* G);

void calcgradSlaterDetOBMEod(const SlaterDet* Q, const SlaterDet* Qp, 
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const gradOneBodyOperator* op, 
			     gradSlaterDet* G);

/// calculate gradient of matrix element with SaterDet Q for 
/// two-body operator op defined for Gaussians
/// and dval(ngauss) will be added up
void calcgradSlaterDetTBME(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   const gradTwoBodyOperator* op, 
			   gradSlaterDet* G);

void calcgradSlaterDetTBMEod(const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const gradTwoBodyOperator* op, 
			     gradSlaterDet* G);

void calcgradSlaterDetTBMErowcol(const SlaterDet* Q, const SlaterDetAux* X,
				 const gradSlaterDetAux* dX,
				 const gradTwoBodyOperator* op, 
				 gradSlaterDet* G,
				 int k, int l);

void calcgradSlaterDetTBMEodrowcol(const SlaterDet* Q, const SlaterDet* Qp,
				   const SlaterDetAux* X,
				   const gradSlaterDetAux* dX,
				   const gradTwoBodyOperator* op, 
				   gradSlaterDet* G,
				   int k, int l);

#endif

/**

  \file SlaterDet.h

  Slater determinant build of Gaussians


  (c) 2003 Thomas Neff

*/


#ifndef _SLATERDET_H
#define _SLATERDET_H

#include <stdio.h>
#include <complex.h>

#include "Gaussian.h"


// in SlaterDet and SlaterDetAux keep space for up 
// to MAXNG Gaussians per nucleon
#define MAXNG 3


/// Slater Determinant

typedef struct {
  int A;		///< number of nucleons
  int Z;		///< number of protons
  int N;		///< number of neutrons
  int ngauss;		///< total number of gaussians
  int* idx;	        ///< index of nucleons
  int* ng;	        ///< number of gaussians per nucleon
  Gaussian* G;		///< Gaussians
} SlaterDet;


/// Auxiliaries for SlaterDeterminants
/// build from |Q> and |Q'>

typedef struct {
  GaussianAux* Gaux;	///< Gaussian Auxiliaries Gaux(ngauss, ngauss)
  complex double* n;	///< single particle ovlap matrix n(A,A)
  complex double* o;	///< inverse ovlap matrix o(A,A)
  complex double ovlap;	///< <Q|Q'>, only defined for non-diagonal
} SlaterDetAux;


/// allocate memory for SlaterDet
void allocateSlaterDet(SlaterDet* Q, int A);

/// initialize Slater determinant Qp
/// with Q
void initSlaterDet(const SlaterDet* Q, SlaterDet* Qp);

/// free memory used by Slater determinant Q
void freeSlaterDet(SlaterDet* Q);

/// write parameters of Slater determinant Q in file
void writeSlaterDet(FILE* fp, const SlaterDet* Q);

/// read parameters of Slater determinant Q from file
/// Q will be initialized
int readSlaterDet(FILE* fp, SlaterDet* Q);

/// read Slater determinant Q from file
/// Q will be initialized
int readSlaterDetfromFile(SlaterDet* Q, const char* fname);

/// copy Slater determinant Q into Qp
/// Qp must have been initialized already
void copySlaterDet(const SlaterDet* Q, SlaterDet* Qp);

/// copy Slater determinant Q into Qp
/// Qp will be initialized
void cloneSlaterDet(const SlaterDet* Q, SlaterDet* Qp);

/// join Slater determinants Qa and Qb into Q
/// Q will be initialized
void joinSlaterDets(const SlaterDet* Qa, const SlaterDet* Qb, SlaterDet* Q);

/// normalize SlaterDet
void normalizeSlaterDet(SlaterDet* Q, SlaterDetAux* X);

/// move Slater determinant by vector d
void moveSlaterDet(SlaterDet* Q, double d[3]);

/// boost Slater determinant with speed v
void boostSlaterDet(SlaterDet* Q, double v[3]);

/// rotate Slater determinant by euler angles
void rotateSlaterDet(SlaterDet* Q, double alpha, double beta, double gamma);

/// move and boost center of mass to origin 
void moveboostSlaterDet(SlaterDet* Q, SlaterDetAux* X);

/// move and boost center of mass to origin, orient main axes 
void moveboostorientSlaterDet(SlaterDet* Q, SlaterDetAux* X);

/// orient main axes
void orientSlaterDet(SlaterDet* Q, SlaterDetAux* X);

/// invert Slater determinant in origin
void invertSlaterDet(SlaterDet* Q);

/// invert Slater determinant with respect to cm position and velocity
void invertcmSlaterDet(SlaterDet* Q, SlaterDetAux* X);

/// time revert Slater determinant
void timerevertSlaterDet(SlaterDet* Q);

/// scale Slater determinant
void scaleSlaterDet(SlaterDet* Q, double kappa);

/// allocate memory for SlaterDetAux
void allocateSlaterDetAux(SlaterDetAux* X, int A);

/// init SlaterDetAux X
void initSlaterDetAux(const SlaterDet* Q, SlaterDetAux* X);

/// free memory used by SlaterDetAux X
void freeSlaterDetAux(SlaterDetAux* X);

/// calculate SlaterDetAux for SlaterDet Q.
/// hermiticity not exploited yet
void calcSlaterDetAux(const SlaterDet* Q, SlaterDetAux* X);


/// calculate SlaterDetAux for SlaterDet's Q and Qp.
/// inversion by LU decomposition 
/// does not work for singular overlap matrix
void calcSlaterDetAuxod(const SlaterDet* Q, const SlaterDet* Qp,
			SlaterDetAux* X);

/// calculate SlaterDetAux for SlaterDet's Q and Qp.
/// inversion by SVD decomposition
/// rank A-1 of overlap matrix assumed, ovlap set to 1.0
/// can only be used with One-Body Operators (GT)
void calcSlaterDetAuxodsingular(const SlaterDet* Q, const SlaterDet* Qp,
				SlaterDetAux* X);

/// calculate matrix element with SaterDet Q for 
/// one-body operator op defined for Gaussians
/// val will be overwritten 
void calcSlaterDetOBME(const SlaterDet* Q, const SlaterDetAux* X,
		       const OneBodyOperator* op, double val[]);


/// calculate matrix element with SaterDet Q for 
/// two-body operator op defined for Gaussians
/// val will be overwritten 
void calcSlaterDetTBME(const SlaterDet* Q, const SlaterDetAux* X,
		       const TwoBodyOperator* op, double val[]);

void calcSlaterDetTBMErowcol(const SlaterDet* Q, const SlaterDetAux* X,
			     const TwoBodyOperator* op, double val[], 
			     int k, int l);


/// calculate off-diagonal matrix element with SaterDets Q  and Qp for 
/// one-body operator op defined for Gaussians
/// val will be overwritten 
void calcSlaterDetOBMEod(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X,
			 const OneBodyOperator* op, complex double val[]);


/// calculate off-diagonal matrix element with SaterDets Q  and Qp for 
/// two-body operator op defined for Gaussians
/// val will be overwritten 
void calcSlaterDetTBMEod(const SlaterDet* Q, const SlaterDet* Qp,
			 const SlaterDetAux* X,
			 const TwoBodyOperator* op, complex double val[]);


void calcSlaterDetTBMEodrowcol(const SlaterDet* Q, const SlaterDet* Qp,
			       const SlaterDetAux* X,
			       const TwoBodyOperator* op, 
			       complex double val[], 
			       int k, int l);


/// calculate Hartree-Fock matrix elements for one-body operator
void calcSlaterDetOBHFMEs(const SlaterDet* Q, const SlaterDetAux* X,
			  const OneBodyOperator* op, void* val);


/// calculate Hartree-Fock matrix elements for two-body operator
void calcSlaterDetTBHFMEs(const SlaterDet* Q, const SlaterDetAux* X,
			  const TwoBodyOperator* op, void* val);


#endif

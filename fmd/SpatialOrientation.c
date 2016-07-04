/**

  \file SpatialOrientation.c

  determine euler angles of spatial orientation
  by diagonalizing the inertia tensor

  determine euler angles of intrinsic coordinate system
  by diagonalizing the < Ji Jj > tensor

  (c) 2003 Thomas Neff

*/


#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "CenterofMass.h"
#include "AngularMomenta.h"
#include "SpatialOrientation.h"

#include "numerics/cmath.h"
#include "numerics/lapack.h"


// inertia tensor 

static void ob_inertia(const double* Xcm,
		       const Gaussian* G1, const Gaussian* G2, 
		       const GaussianAux* X, complex double T[9])
{
  complex double rho2cm = 0.0;
  int i,j;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  for (j=0; j<3; j++)
    for (i=0; i<3; i++) {
      T[i+j*3] -= (X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*X->Q;
      if (i==j) T[i+j*3] += (2*X->alpha + rho2cm)*X->Q;
    } 
}


void calcInertiaTensor(const SlaterDet* Q, const SlaterDetAux* X,
		       double T[9])
{
  double Xcm[3];
  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_inertia = { dim: 9, opt: 1, par: Xcm, me: ob_inertia};
  calcSlaterDetOBME(Q, X, &op_ob_inertia, T);
}	


// octupole tensor

static void ob_octupole(const double* Xcm,
			const Gaussian*G1, const Gaussian* G2, 
			const GaussianAux* X, 
			complex double o[27])
{
  complex double rho2cm = 0.0;
  int i,j,k;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i] - Xcm[i]);

  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++) {
	o[i+j*3+k*9] += 
	  5*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*(X->rho[k]-Xcm[k])*X->Q;

	if (i==j) o[i+j*3+k*9] += - rho2cm* (X->rho[k]-Xcm[k])* X->Q; 
	if (j==k) o[i+j*3+k*9] += - rho2cm* (X->rho[i]-Xcm[i])* X->Q; 
	if (k==i) o[i+j*3+k*9] += - rho2cm* (X->rho[j]-Xcm[j])* X->Q; 
      }
}


void calcOctupoleTensor(const SlaterDet* Q, const SlaterDetAux* X, 
		      double O[27])
{
  double Xcm[3];
  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_octupole = {dim: 27, opt: 1, par: Xcm, me: ob_octupole};

  calcSlaterDetOBME(Q, X, &op_ob_octupole, O);
}


static void RotationMatrixToEulerAngles(double R[3][3], 
					double* a, double* b, double* c)
{
  const double eps=10.e-10;
  double arg;

  if (fabs(R[2][0]) < eps && fabs(R[2][1]) < eps) {
    *a = atan2(R[0][1], R[0][0]);
    arg = R[2][2];
    if (arg > 1.0) arg=1.0;
    if (arg < -1.0) arg=-1.0;
    *b = acos(arg);
    *c = 0.0;
  } else {
    *a = atan2(R[2][1], R[2][0]);
    arg = R[2][2];
    if (arg > 1.0) arg=1.0;
    if (arg < -1.0) arg=-1.0;
    *b = acos(arg);
    *c = atan2(R[1][2], -R[0][2]);
    
    if (fabs(R[1][2]) > fabs(R[0][2])) {
      if ((0.0 > sin(*b)*sin(*c)) != (0.0 > R[1][2])) {
        *b = -(*b);
      }
    } else {
      if ((0.0 > sin(*b)*cos(*c)) != (0.0 > -R[0][2])) {
        *b = -(*b);
      }
    }
  }

  // return alpha, gamma [0, 2Pi], beta [0,Pi]
  if (*a < 0.0) *a += 2*M_PI;
  if (*c < 0.0) *c += 2*M_PI;
}


static double det(double R[3][3])
{
  return (R[0][0]*R[1][1]*R[2][2] + R[0][1]*R[1][2]*R[2][0] + 
	  R[0][2]*R[1][0]*R[2][1] - R[0][0]*R[1][2]*R[2][1] - 
	  R[0][1]*R[1][0]*R[2][2] - R[0][2]*R[1][1]*R[2][0]);
}



void calcSpatialOrientation(const SlaterDet* Q, const SlaterDetAux* X, 
			    double* alpha, double* beta, double* gamma)
{
  double T[3][3];

  calcInertiaTensor(Q, X, T);

  char JOBZ='V';
  char UPLO='U';
  int N=3;
  double W[3];
  double WORK[15];
  int LWORK=15;
  int INFO;

  FORTRAN(dsyev)(&JOBZ, &UPLO, &N, (double*)T, &N, W, WORK, &LWORK, &INFO);

  // align supposed symmetry axis with z axis

  int i;
  double maxe=0.0;

  for (i=0; i<3; i++)
    maxe = fmax(maxe, fabs(W[i]));

  double dif[3], mindif = 1.0;
  int imin=-1;

  for (i=0; i<3; i++) {
    dif[i] = fabs(W[(i+1)%3]-W[(i+2)%3])/maxe;
    if (dif[i] < mindif) {
      mindif = dif[i];
      imin = i;
    }
  }

  // align z-axis with symmetry axis
  double dummy;
  for (i=0; i<3; i++) {
    dummy = T[2][i];
    T[2][i] = T[imin][i];
    T[imin][i] = dummy;
  }

  if (det(T) < 0.0)
    for (i=0; i<3; i++)
      T[0][i] *= -1;


  RotationMatrixToEulerAngles(T, alpha, beta, gamma);
}


void calcOrientedOrientation(const SlaterDet* Q, const SlaterDetAux* X, 
			     double* alpha, double* beta, double* gamma)
{
  double T[3][3];
  double O[3][3][3];
  int i,j,k;

  calcInertiaTensor(Q, X, T);
  calcOctupoleTensor(Q, X, O);

  char JOBZ='V';
  char UPLO='U';
  int N=3;
  double W[3];
  double WORK[15];
  int LWORK=15;
  int INFO;

  FORTRAN(dsyev)(&JOBZ, &UPLO, &N, (double*)T, &N, W, WORK, &LWORK, &INFO);

  // axis with smallest eigenvalue to become z axis

  int imin=0;
  double minw=fabs(W[0]);

  for (i=1; i<3; i++)
    if (fabs(W[i]) < minw) {
      imin = i;
      minw = fabs(W[i]);
    }

  // align z-axis with symmetry axis
  double dummy;
  for (i=0; i<3; i++) {
    dummy = T[2][i];
    T[2][i] = T[imin][i];
    T[imin][i] = dummy;
  }

  // we want a pure rotation
  if (det(T) < 0.0)
    for (i=0; i<3; i++)
      T[0][i] *= -1;

  // inertia tensor is invariant under Ry(180deg)
  // use sign of octupole moment to select orientation 
  double Ozzz = 0.0;
  
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      for (k=0; k<3; k++)
	Ozzz += T[2][i]*T[2][j]*T[2][k]*O[i][j][k];

  fprintf(stderr, "... Ozzz: %8.5f\n", Ozzz);

  // Ry(180deg)
  if (Ozzz < 0.0) {
    for (i=0; i<3; i++) {
      T[0][i] *= -1;
      T[2][i] *= -1;
    }
  }

  RotationMatrixToEulerAngles(T, alpha, beta, gamma);
}


/// spin orientation with respect to z-axis
void calcSpinOrientation(const SlaterDet* Q, const SlaterDetAux* X, 
			 double* alpha, double* beta, double* gamma)
{
  double J[3];

  calcJ(Q, X, J);

  double Jxy=sqrt(J[0]*J[0]+J[1]*J[1]);
  double J2=sqrt(J[0]*J[0]+J[1]*J[1]+J[2]*J[2]);

  *alpha = (J[1] > 0.0 ? +1 : -1)* acos(J[0]/Jxy);
  *beta = acos(J[2]/J2);
  *gamma = 0.0;
}

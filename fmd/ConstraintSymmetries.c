/**

  \file ConstraintSymmetries.c

  constrain intrinsic state regarding discrete symmetries


  (c) 2009 Thomas Neff

*/

#include <assert.h>
#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Parity.h"

#include "ConstraintSymmetries.h"


// positive parity
Constraint ConstraintParityP = {
  val : 0.0,
  label : "ParP",
  me : calcConstraintParP,
  gradme : calcgradConstraintParP,
  output : outputConstraintZero
};

// negative parity
Constraint ConstraintParityN = {
  val : 0.0,
  label : "ParN",
  me : calcConstraintParN,
  gradme : calcgradConstraintParN,
  output : outputConstraintZero
};

Constraint ConstraintParityCharge = {
  val : 0.0,
  label: "ParCh",
  me : calcConstraintParCh,
  gradme : calcgradConstraintParCh,
  output : outputConstraintZero
};


// rotation by 120deg around z-axis
Constraint ConstraintD3H = {
  val : 0.0,
  label : "D3H",
  me : calcConstraintD3H,
  gradme : calcgradConstraintD3H,
  output : outputConstraintZero
};

// reflection on xz-plane
Constraint ConstraintSxz = {
  val : 0.0,
  label : "Sxz",
  me : calcConstraintSxz,
  gradme : calcgradConstraintSxz,
  output : outputConstraintZero
};

// reflection on xz-plane
Constraint ConstraintSyz = {
  val : 0.0,
  label : "Syz",
  me : calcConstraintSyz,
  gradme : calcgradConstraintSyz,
  output : outputConstraintZero
};


// workspace
static SlaterDet SQ;
static SlaterDet SdagQ;
static gradSlaterDet dN;
static gradSlaterDet dS;
static gradSlaterDet dSdag;

static void initworkspace(const SlaterDet* Q)
{
  initSlaterDet(Q, &SQ);
  initSlaterDet(Q, &SdagQ);

  initgradSlaterDet(Q, &dN);
  initgradSlaterDet(Q, &dS);
  initgradSlaterDet(Q, &dSdag);  
}                    


// parity operation is hermitian

// constraint on positive parity
void calcConstraintParP(const SlaterDet* Q, const SlaterDetAux* X,
                        double* parpconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);
    
  copySlaterDet(Q, &SQ);
  // parity operation
  invertSlaterDet(&SQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;

  calcSlaterDetAuxod(Q, &SQ, X);
  double S = X->ovlap;

  *parpconstraint = S/N-1.0;
}

void calcgradConstraintParP(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dparpconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // parity operation
  invertSlaterDet(&SQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;
  calcgradSlaterDetAuxod(Q, Q, X, dX);
  calcgradOvlapod(Q, Q, X, dX, &dN);

  calcSlaterDetAuxod(Q, &SQ, X);
  double S = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SQ, X, dX);
  calcgradOvlapod(Q, &SQ, X, dX, &dS);

  addmulttogradSlaterDet(dparpconstraint, &dS, 1.0/N);
  addmulttogradSlaterDet(dparpconstraint, &dN, -S/(N*N));
  dparpconstraint->val = S/N - 1.0;

}


// constraint on negative parity
void calcConstraintParN(const SlaterDet* Q, const SlaterDetAux* X,
                        double* parnconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);
    
  copySlaterDet(Q, &SQ);
  // parity operation
  invertSlaterDet(&SQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;

  calcSlaterDetAuxod(Q, &SQ, X);
  double S = X->ovlap;

  *parnconstraint = S/N + 1.0;
}

void calcgradConstraintParN(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dparnconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // parity operation
  invertSlaterDet(&SQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;
  calcgradSlaterDetAuxod(Q, Q, X, dX);
  calcgradOvlapod(Q, Q, X, dX, &dN);

  calcSlaterDetAuxod(Q, &SQ, X);
  double S = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SQ, X, dX);
  calcgradOvlapod(Q, &SQ, X, dX, &dS);

  addmulttogradSlaterDet(dparnconstraint, &dS, 1.0/N);
  addmulttogradSlaterDet(dparnconstraint, &dN, -S/(N*N));
  dparnconstraint->val = S/N + 1.0;

}


// rotate by pi in isospin space
static void RisoppiSlaterDet(SlaterDet* Q)
{
  assert(Q->Z == Q->N);
  
  int i;
  for (i=0; i<Q->ngauss; i++) {
    if (Q->G[i].xi == +1) {
      Q->G[i].xi = -1;
    } else {
      Q->G[i].xi = +1;
      Q->G[i].chi[0] *= -1.0;
      Q->G[i].chi[1] *= -1.0;
    }
  }
}

// rotate by -pi in isospin space
static void RisompiSlaterDet(SlaterDet* Q)
{
  assert(Q->Z == Q->N);
  
  int i;
  for (i=0; i<Q->ngauss; i++) {
    if (Q->G[i].xi == -1) {
      Q->G[i].xi = +1;
    } else {
      Q->G[i].xi = -1;
      Q->G[i].chi[0] *= -1.0;
      Q->G[i].chi[1] *= -1.0;
    }
  }
}


// constraint on symmetry under parity and charge
void calcConstraintParCh(const SlaterDet* Q, const SlaterDetAux* X,
                         double* parchconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // parity and charge
  invertSlaterDet(&SQ);
  RisoppiSlaterDet(&SQ);

  copySlaterDet(Q, &SdagQ);
  // parity and charge
  invertSlaterDet(&SdagQ);
  RisompiSlaterDet(&SdagQ);
  
  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;

  *parchconstraint = (Sdag*S)/(N*N) - 1.0;
}


void calcgradConstraintParCh(const SlaterDet* Q, const SlaterDetAux* X,
                             const gradSlaterDetAux* dX,
                             gradSlaterDet* dparchconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // parity and charge
  invertSlaterDet(&SQ);
  RisoppiSlaterDet(&SQ);

  copySlaterDet(Q, &SdagQ);
  // parity and charge
  invertSlaterDet(&SdagQ);
  RisompiSlaterDet(&SdagQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;
  calcgradSlaterDetAuxod(Q, Q, X, dX);
  calcgradOvlapod(Q, Q, X, dX, &dN);

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SQ, X, dX);
  calcgradOvlapod(Q, &SQ, X, dX, &dS);

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SdagQ, X, dX);
  calcgradOvlapod(Q, &SdagQ, X, dX, &dSdag);

  addmulttogradSlaterDet(dparchconstraint, &dSdag, S/(N*N));
  addmulttogradSlaterDet(dparchconstraint, &dS, Sdag/(N*N));
  addmulttogradSlaterDet(dparchconstraint, &dN, -2*Sdag*S/(N*N*N));
  dparchconstraint->val = Sdag*S/(N*N) - 1.0;
}


// constraint on rotation by 120deg around z-axis
void calcConstraintD3H(const SlaterDet* Q, const SlaterDetAux* X,
                       double* pard3hconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);
    
  copySlaterDet(Q, &SQ);
  // rotation by 120deg
  rotateSlaterDet(&SQ, 0.0, 0.0, 2.0/3.0*M_PI);

  copySlaterDet(Q, &SdagQ);
  // rotation by -120deg
  rotateSlaterDet(&SdagQ, 0.0, 0.0, -2.0/3.0*M_PI);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;

  // do not ignore phase
  // *pard3hconstraint = (Sdag/N - 1.0)*(S/N - 1.0);

  // ignore phase
  *pard3hconstraint = Sdag*S/(N*N) - 1.0;
}

void calcgradConstraintD3H(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dd3hconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // rotation by 120deg
  rotateSlaterDet(&SQ, 0.0, 0.0, 2.0/3.0*M_PI);

  copySlaterDet(Q, &SdagQ);
  // rotation by -120deg
  rotateSlaterDet(&SdagQ, 0.0, 0.0, -2.0/3.0*M_PI);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;
  calcgradSlaterDetAuxod(Q, Q, X, dX);
  calcgradOvlapod(Q, Q, X, dX, &dN);

  calcSlaterDetAuxod(Q, &SQ, X);
  double S = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SQ, X, dX);
  calcgradOvlapod(Q, &SQ, X, dX, &dS);

  calcSlaterDetAuxod(Q, &SdagQ, X);
  double Sdag = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SdagQ, X, dX);
  calcgradOvlapod(Q, &SdagQ, X, dX, &dSdag);

  // ignore phase
  addmulttogradSlaterDet(dd3hconstraint, &dSdag, S/(N*N));
  addmulttogradSlaterDet(dd3hconstraint, &dS, Sdag/(N*N));
  addmulttogradSlaterDet(dd3hconstraint, &dN, -2*Sdag*S/(N*N*N));
  dd3hconstraint->val = Sdag*S/(N*N) - 1.0;
}


// constraint on reflection symmetry on xz-plane
void calcConstraintSxz(const SlaterDet* Q, const SlaterDetAux* X,
                       double* sxzconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);
    
  copySlaterDet(Q, &SQ);
  // invert and rotate by 180deg around y-axis
  invertSlaterDet(&SQ);
  rotateSlaterDet(&SQ, 0.0, M_PI, 0.0);

  copySlaterDet(Q, &SdagQ);
  // invert and rotate by -180deg around y-axis
  rotateSlaterDet(&SdagQ, 0.0, -M_PI, 0.0);
  invertSlaterDet(&SdagQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;

  *sxzconstraint = (Sdag*S)/(N*N) - 1.0;
}

void calcgradConstraintSxz(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dsxzconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // invert and rotate by 180deg around y-axis
  invertSlaterDet(&SQ);
  rotateSlaterDet(&SQ, 0.0, M_PI, 0.0);

  copySlaterDet(Q, &SdagQ);
  // invert and rotate by -180deg around y-axis
  rotateSlaterDet(&SdagQ, 0.0, -M_PI, 0.0);
  invertSlaterDet(&SdagQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;
  calcgradSlaterDetAuxod(Q, Q, X, dX);
  calcgradOvlapod(Q, Q, X, dX, &dN);

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SQ, X, dX);
  calcgradOvlapod(Q, &SQ, X, dX, &dS);

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SdagQ, X, dX);
  calcgradOvlapod(Q, &SdagQ, X, dX, &dSdag);

  addmulttogradSlaterDet(dsxzconstraint, &dSdag, S/(N*N));
  addmulttogradSlaterDet(dsxzconstraint, &dS, Sdag/(N*N));
  addmulttogradSlaterDet(dsxzconstraint, &dN, -(2*Sdag*S)/(N*N*N));
  dsxzconstraint->val = (Sdag*S)/(N*N) - 1.0;

}


// constraint on reflection symmetry on xz-plane
void calcConstraintSyz(const SlaterDet* Q, const SlaterDetAux* X,
                       double* syzconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);
    
  copySlaterDet(Q, &SQ);
  // invert and rotate by 180deg around x-axis
  invertSlaterDet(&SQ);
  rotateSlaterDet(&SQ, 0.5*M_PI, M_PI, -0.5*M_PI);

  copySlaterDet(Q, &SdagQ);
  // invert and rotate by -180deg around x-axis
  rotateSlaterDet(&SdagQ, 0.5*M_PI, -M_PI, -0.5*M_PI);
  invertSlaterDet(&SdagQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;

  *syzconstraint = (Sdag*S)/(N*N) - 1.0;
}

void calcgradConstraintSyz(const SlaterDet* Q, const SlaterDetAux* X,
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dsyzconstraint)
{
  if (SQ.A == 0)
    initworkspace(Q);

  copySlaterDet(Q, &SQ);
  // invert and rotate by 180deg around x-axis
  invertSlaterDet(&SQ);
  rotateSlaterDet(&SQ, 0.5*M_PI, M_PI, -0.5*M_PI);

  copySlaterDet(Q, &SdagQ);
  // invert and rotate by -180deg around x-axis
  rotateSlaterDet(&SdagQ, 0.5*M_PI, -M_PI, -0.5*M_PI);
  invertSlaterDet(&SdagQ);

  calcSlaterDetAuxod(Q, Q, X);
  double N = X->ovlap;
  calcgradSlaterDetAuxod(Q, Q, X, dX);
  calcgradOvlapod(Q, Q, X, dX, &dN);

  calcSlaterDetAuxod(Q, &SQ, X);
  complex double S = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SQ, X, dX);
  calcgradOvlapod(Q, &SQ, X, dX, &dS);

  calcSlaterDetAuxod(Q, &SdagQ, X);
  complex double Sdag = X->ovlap;
  calcgradSlaterDetAuxod(Q, &SdagQ, X, dX);
  calcgradOvlapod(Q, &SdagQ, X, dX, &dSdag);

  addmulttogradSlaterDet(dsyzconstraint, &dSdag, S/(N*N));
  addmulttogradSlaterDet(dsyzconstraint, &dS, Sdag/(N*N));
  addmulttogradSlaterDet(dsyzconstraint, &dN, -(2*Sdag*S)/(N*N*N));
  dsyzconstraint->val = (Sdag*S)/(N*N) - 1.0;

}



double outputConstraintZero(double val)
{
  return (val);
}

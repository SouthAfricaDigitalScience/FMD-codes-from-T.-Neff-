/**

  \file ConstraintNOsci.c

  Constrain number of oscillator quanta
  for an isotropic oscillator


  (c) 2006 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "gradCMSlaterDet.h"
#include "CenterofMass.h"

#include "Constraint.h"
#include "ConstraintNOsci.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


Constraint ConstraintNOsci = {
  val : 0.0,
  label : "NOsci",
  me : calcConstraintNOsci,
  gradme : calcgradConstraintNOsci,
  output : outputConstraintNOsci
};

Constraint ConstraintPNOsci = {
  val : 0.0,
  label : "PNOsci",
  me : calcConstraintPNOsci,
  gradme : calcgradConstraintPNOsci,
  output : outputConstraintNOsci
};

Constraint ConstraintNNOsci = {
  val : 0.0,
  label : "NNOsci",
  me : calcConstraintNNOsci,
  gradme : calcgradConstraintNNOsci,
  output : outputConstraintNOsci
};


typedef struct {
  double Xcm[3];
  double Vcm[3];
  double x2;
  double k2;
} noscipara;


// 0hbw oscillator quanta

int nq[] = {  0,                                 // no nucleons
              0, 0,                              // s-shell
              1, 2, 3, 4, 5, 6,                  // p-shell
              8,10,12,14,16,18,20,22,24,26,      // sd-shell
              29,32,35,38,41,44,47,50,53,56,     // pf-shell
              59,62,65,68,71,74,77,80,83,86 };

static inline double nquanta0(const SlaterDet* Q)
{
  return (1.5*Q->A + nq[Q->Z] + nq[Q->N]);
}


static void ob_x2(const double* Xcm,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, 
		  complex double *x2cm)
{
  int i;
  complex double rho2cm;

  rho2cm = 3*X->alpha;
  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  *x2cm += rho2cm* X->Q;
}

static void ob_px2(const double* Xcm,
                   const Gaussian*G1, const Gaussian* G2, 
                   const GaussianAux* X, 
                   complex double *x2cm)
{
  int i;
  complex double rho2cm;

  if (G1->xi == 1) {

    rho2cm = 3*X->alpha;
    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    *x2cm += rho2cm* X->Q;
  }
}

static void ob_nx2(const double* Xcm,
                   const Gaussian*G1, const Gaussian* G2, 
                   const GaussianAux* X, 
                   complex double *x2cm)
{
  int i;
  complex double rho2cm;

  if (G1->xi == -1) {

    rho2cm = 3*X->alpha;
    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    *x2cm += rho2cm* X->Q;
  }
}


static void ob_k2(const double* Vcm,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, 
		  complex double *k2cm)
{
  int i;
  complex double pi2cm;

  pi2cm = 3*X->lambda;
  for (i=0; i<3; i++)
    pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);

  *k2cm += pi2cm* X->Q;
}

static void ob_pk2(const double* Vcm,
                   const Gaussian*G1, const Gaussian* G2, 
                   const GaussianAux* X, 
                   complex double *k2cm)
{
  int i;
  complex double pi2cm;

  if (G1->xi == 1) {

    pi2cm = 3*X->lambda;
    for (i=0; i<3; i++)
      pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);

    *k2cm += pi2cm* X->Q;
  }
}

static void ob_nk2(const double* Vcm,
                   const Gaussian*G1, const Gaussian* G2, 
                   const GaussianAux* X, 
                   complex double *k2cm)
{
  int i;
  complex double pi2cm;

  if (G1->xi == -1) {

    pi2cm = 3*X->lambda;
    for (i=0; i<3; i++)
      pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);

    *k2cm += pi2cm* X->Q;
  }
}


static void gob_nosci(const noscipara* oscipar,
		      const Gaussian* G1, const Gaussian* G2,
		      const GaussianAux* X, const gradGaussianAux* dX,
		      complex double* nosci, gradGaussian* dnosci)
{
  double *Xcm = oscipar->Xcm;
  double *Vcm = oscipar->Vcm;
  double x2cm = oscipar->x2;
  double k2cm = oscipar->k2;

  int i;
  complex double rho2cm, pi2cm;

  rho2cm = 3*X->alpha;
  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  pi2cm = 3*X->lambda;
  for (i=0; i<3; i++)
    pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);

  gradGaussian drho2cm, dpi2cm;
 
  drho2cm.a = 3*dX->dalpha;
  for (i=0; i<3; i++)
    drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
  for (i=0; i<3; i++)
    drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

  dpi2cm.a = 3*dX->dlambda;
  for (i=0; i<3; i++)
    dpi2cm.a += 2*(X->pi[i]-mass(G1->xi)*Vcm[i])*dX->dpi.a[i];
  for (i=0; i<3; i++)
    dpi2cm.b[i] = 2*(X->pi[i]-mass(G1->xi)*Vcm[i])*dX->dpi.b;

  *nosci += (rho2cm*k2cm + x2cm*pi2cm)*X->Q;

  for (i=0; i<2; i++)
    dnosci->chi[i] += 
      (rho2cm*k2cm + x2cm* pi2cm)* dX->dQ.chi[i];
  dnosci->a += 
    ((drho2cm.a*k2cm + x2cm*dpi2cm.a)* X->Q +
     (rho2cm*k2cm + x2cm*pi2cm)* dX->dQ.a);
  for (i=0; i<3; i++)
    dnosci->b[i] += 
      ((drho2cm.b[i]*k2cm + x2cm*dpi2cm.b[i])* X->Q +
       (rho2cm*k2cm + x2cm*pi2cm)* dX->dQ.b[i]);
}

static void gob_pnosci(const noscipara* oscipar,
                       const Gaussian* G1, const Gaussian* G2,
                       const GaussianAux* X, const gradGaussianAux* dX,
                       complex double* nosci, gradGaussian* dnosci)
{
  double *Xcm = oscipar->Xcm;
  double *Vcm = oscipar->Vcm;
  double x2cm = oscipar->x2;
  double k2cm = oscipar->k2;

  int i;
  complex double rho2cm, pi2cm;

  if (G1->xi == 1) {

    rho2cm = 3*X->alpha;
    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    pi2cm = 3*X->lambda;
    for (i=0; i<3; i++)
      pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);

    gradGaussian drho2cm, dpi2cm;
 
    drho2cm.a = 3*dX->dalpha;
    for (i=0; i<3; i++)
      drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
    for (i=0; i<3; i++)
      drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

    dpi2cm.a = 3*dX->dlambda;
    for (i=0; i<3; i++)
      dpi2cm.a += 2*(X->pi[i]-mass(G1->xi)*Vcm[i])*dX->dpi.a[i];
    for (i=0; i<3; i++)
      dpi2cm.b[i] = 2*(X->pi[i]-mass(G1->xi)*Vcm[i])*dX->dpi.b;

    *nosci += (rho2cm*k2cm + x2cm*pi2cm)*X->Q;

    for (i=0; i<2; i++)
      dnosci->chi[i] += 
        (rho2cm*k2cm + x2cm* pi2cm)* dX->dQ.chi[i];
    dnosci->a += 
      ((drho2cm.a*k2cm + x2cm*dpi2cm.a)* X->Q +
       (rho2cm*k2cm + x2cm*pi2cm)* dX->dQ.a);
    for (i=0; i<3; i++)
      dnosci->b[i] += 
        ((drho2cm.b[i]*k2cm + x2cm*dpi2cm.b[i])* X->Q +
         (rho2cm*k2cm + x2cm*pi2cm)* dX->dQ.b[i]);
  }
}

static void gob_nnosci(const noscipara* oscipar,
                       const Gaussian* G1, const Gaussian* G2,
                       const GaussianAux* X, const gradGaussianAux* dX,
                       complex double* nosci, gradGaussian* dnosci)
{
  double *Xcm = oscipar->Xcm;
  double *Vcm = oscipar->Vcm;
  double x2cm = oscipar->x2;
  double k2cm = oscipar->k2;

  int i;
  complex double rho2cm, pi2cm;

  if (G1->xi == -1) {

    rho2cm = 3*X->alpha;
    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    pi2cm = 3*X->lambda;
    for (i=0; i<3; i++)
      pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);

    gradGaussian drho2cm, dpi2cm;
 
    drho2cm.a = 3*dX->dalpha;
    for (i=0; i<3; i++)
      drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
    for (i=0; i<3; i++)
      drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

    dpi2cm.a = 3*dX->dlambda;
    for (i=0; i<3; i++)
      dpi2cm.a += 2*(X->pi[i]-mass(G1->xi)*Vcm[i])*dX->dpi.a[i];
    for (i=0; i<3; i++)
      dpi2cm.b[i] = 2*(X->pi[i]-mass(G1->xi)*Vcm[i])*dX->dpi.b;

    *nosci += (rho2cm*k2cm + x2cm*pi2cm)*X->Q;

    for (i=0; i<2; i++)
      dnosci->chi[i] += 
        (rho2cm*k2cm + x2cm* pi2cm)* dX->dQ.chi[i];
    dnosci->a += 
      ((drho2cm.a*k2cm + x2cm*dpi2cm.a)* X->Q +
       (rho2cm*k2cm + x2cm*pi2cm)* dX->dQ.a);
    for (i=0; i<3; i++)
      dnosci->b[i] += 
        ((drho2cm.b[i]*k2cm + x2cm*dpi2cm.b[i])* X->Q +
         (rho2cm*k2cm + x2cm*pi2cm)* dX->dQ.b[i]);
  }
}


static void gcmob_nosci(const noscipara* oscipar,
			const Gaussian* G1, const Gaussian* G2,
			const GaussianAux* X,
			gradCM* dCMnosci)
{
  double *Xcm = oscipar->Xcm;
  double *Vcm = oscipar->Vcm;
  double x2cm = oscipar->x2;
  double k2cm = oscipar->k2;

  int i;

  for (i=0; i<3; i++)
    dCMnosci->X[i] -=
      (X->rho[i]-Xcm[i])*X->Q* k2cm;
  for (i=0; i<3; i++)
    dCMnosci->V[i] -=
      x2cm* mass(G1->xi)*(X->pi[i]-mass(G1->xi)*Vcm[i])*X->Q;
}

static void gcmob_pnosci(const noscipara* oscipar,
                         const Gaussian* G1, const Gaussian* G2,
                         const GaussianAux* X,
                         gradCM* dCMnosci)
{
  double *Xcm = oscipar->Xcm;
  double *Vcm = oscipar->Vcm;
  double x2cm = oscipar->x2;
  double k2cm = oscipar->k2;

  int i;

  if (G1->xi == 1) {

    for (i=0; i<3; i++)
      dCMnosci->X[i] -=
        (X->rho[i]-Xcm[i])*X->Q* k2cm;
    for (i=0; i<3; i++)
      dCMnosci->V[i] -=
        x2cm* mass(G1->xi)*(X->pi[i]-mass(G1->xi)*Vcm[i])*X->Q;
  }
}

static void gcmob_nnosci(const noscipara* oscipar,
                         const Gaussian* G1, const Gaussian* G2,
                         const GaussianAux* X,
                         gradCM* dCMnosci)
{
  double *Xcm = oscipar->Xcm;
  double *Vcm = oscipar->Vcm;
  double x2cm = oscipar->x2;
  double k2cm = oscipar->k2;

  int i;

  if (G1->xi == -1) {

    for (i=0; i<3; i++)
      dCMnosci->X[i] -=
        (X->rho[i]-Xcm[i])*X->Q* k2cm;
    for (i=0; i<3; i++)
      dCMnosci->V[i] -=
        x2cm* mass(G1->xi)*(X->pi[i]-mass(G1->xi)*Vcm[i])*X->Q;
  }
}


#define SQR(x) (x)*(x)

void calcConstraintNOsci(const SlaterDet* Q, const SlaterDetAux* X, 
			 double* nosci)
{
  double Xcm[3], Vcm[3];
  double x2, k2;

  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  OneBodyOperator op_ob_x2 = {dim: 1, opt: 1, par: Xcm, me: ob_x2};
  OneBodyOperator op_ob_k2 = {dim: 1, opt: 1, par: Vcm, me: ob_k2};

  calcSlaterDetOBME(Q, X, &op_ob_x2, &x2);
  calcSlaterDetOBME(Q, X, &op_ob_k2, &k2);
  
  *nosci = x2*k2;
}

void calcConstraintPNOsci(const SlaterDet* Q, const SlaterDetAux* X, 
                          double* nosci)
{
  double Xcm[3], Vcm[3];
  double x2, k2;

  // use proton cm instead ?

  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  OneBodyOperator op_ob_px2 = {dim: 1, opt: 1, par: Xcm, me: ob_px2};
  OneBodyOperator op_ob_pk2 = {dim: 1, opt: 1, par: Vcm, me: ob_pk2};

  calcSlaterDetOBME(Q, X, &op_ob_px2, &x2);
  calcSlaterDetOBME(Q, X, &op_ob_pk2, &k2);
  
  *nosci = x2*k2;
}

void calcConstraintNNOsci(const SlaterDet* Q, const SlaterDetAux* X, 
                          double* nosci)
{
  double Xcm[3], Vcm[3];
  double x2, k2;

  // use neutron cm instead ?

  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  OneBodyOperator op_ob_nx2 = {dim: 1, opt: 1, par: Xcm, me: ob_nx2};
  OneBodyOperator op_ob_nk2 = {dim: 1, opt: 1, par: Vcm, me: ob_nk2};

  calcSlaterDetOBME(Q, X, &op_ob_nx2, &x2);
  calcSlaterDetOBME(Q, X, &op_ob_nk2, &k2);
  
  *nosci = x2*k2;
}


void calcgradConstraintNOsci(const SlaterDet* Q, const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     gradSlaterDet* dnosci)
{
  double Xcm[3], Vcm[3];
  double x2, k2;
  noscipara opar;
  int i;

  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  OneBodyOperator op_ob_x2 = {dim: 1, opt: 1, par: Xcm, me: ob_x2};
  OneBodyOperator op_ob_k2 = {dim: 1, opt: 1, par: Vcm, me: ob_k2};

  calcSlaterDetOBME(Q, X, &op_ob_x2, &x2);
  calcSlaterDetOBME(Q, X, &op_ob_k2, &k2);

  for (i=0; i<3; i++)
    opar.Xcm[i] = Xcm[i];
  for (i=0; i<3; i++)
    opar.Vcm[i] = Vcm[i];
  opar.x2 = x2;
  opar.k2 = k2;

  gradOneBodyOperator gop_ob_nosci = { opt: 1, par: &opar, me: gob_nosci };

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_nosci, dnosci);

  OneBodyOperator gcmop_ob_nosci = { dim: 6, opt: 1, par: &opar, 
				     me: gcmob_nosci };

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_nosci, dnosci);

  dnosci->val = x2*k2;
}

void calcgradConstraintPNOsci(const SlaterDet* Q, const SlaterDetAux* X,
                              const gradSlaterDetAux* dX,
                              gradSlaterDet* dnosci)
{
  double Xcm[3], Vcm[3];
  double x2, k2;
  noscipara opar;
  int i;

  // proton cm instead ?

  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  OneBodyOperator op_ob_px2 = {dim: 1, opt: 1, par: Xcm, me: ob_px2};
  OneBodyOperator op_ob_pk2 = {dim: 1, opt: 1, par: Vcm, me: ob_pk2};

  calcSlaterDetOBME(Q, X, &op_ob_px2, &x2);
  calcSlaterDetOBME(Q, X, &op_ob_pk2, &k2);

  for (i=0; i<3; i++)
    opar.Xcm[i] = Xcm[i];
  for (i=0; i<3; i++)
    opar.Vcm[i] = Vcm[i];
  opar.x2 = x2;
  opar.k2 = k2;

  gradOneBodyOperator gop_ob_pnosci = { opt: 1, par: &opar, me: gob_pnosci };

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_pnosci, dnosci);

  OneBodyOperator gcmop_ob_pnosci = { dim: 6, opt: 1, par: &opar, 
                                      me: gcmob_pnosci };

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_pnosci, dnosci);

  dnosci->val = x2*k2;
}

void calcgradConstraintNNOsci(const SlaterDet* Q, const SlaterDetAux* X,
                              const gradSlaterDetAux* dX,
                              gradSlaterDet* dnosci)
{
  double Xcm[3], Vcm[3];
  double x2, k2;
  noscipara opar;
  int i;

  // neutron cm instead ?

  calcCMPosition(Q, X, Xcm);
  calcCMVelocity(Q, X, Vcm);

  OneBodyOperator op_ob_nx2 = {dim: 1, opt: 1, par: Xcm, me: ob_nx2};
  OneBodyOperator op_ob_nk2 = {dim: 1, opt: 1, par: Vcm, me: ob_nk2};

  calcSlaterDetOBME(Q, X, &op_ob_nx2, &x2);
  calcSlaterDetOBME(Q, X, &op_ob_nk2, &k2);

  for (i=0; i<3; i++)
    opar.Xcm[i] = Xcm[i];
  for (i=0; i<3; i++)
    opar.Vcm[i] = Vcm[i];
  opar.x2 = x2;
  opar.k2 = k2;

  gradOneBodyOperator gop_ob_nnosci = { opt: 1, par: &opar, me: gob_nnosci };

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_nnosci, dnosci);

  OneBodyOperator gcmop_ob_nnosci = { dim: 6, opt: 1, par: &opar, 
                                      me: gcmob_nnosci };

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_nnosci, dnosci);

  dnosci->val = x2*k2;
}


double outputConstraintNOsci(double val)
{
  return sqrt(val);
}




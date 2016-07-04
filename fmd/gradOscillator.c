/**

  \file gradOscillator.c

  external quadratic and quartic oscillator


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "gradOscillator.h"

#include "numerics/cmath.h"
#include "misc/physics.h"



static void gob_osci2(double* kosci,
		      const Gaussian* G1, const Gaussian* G2, 
		      const GaussianAux* X, const gradGaussianAux* dX, 
		      complex double* h, gradGaussian* dh)
{	
  int i;

  *h += 0.5* *kosci* (3*X->alpha + X->rho2)* X->Q;

  for (i=0; i<2; i++)
    dh->chi[i] += 0.5* *kosci* X->T* dX->dS.chi[i]* 
      (3*X->alpha + X->rho2)* X->R;

  dh->a += 0.5* *kosci* X->T* X->S* 
    ((3*dX->dalpha+dX->drho2.a)* X->R + (3*X->alpha + X->rho2)*dX->dR.a);
  for (i=0; i<3; i++)
    dh->b[i] += 0.5* *kosci* X->T* X->S* (dX->drho2.b[i]* X->R +
					  (3*X->alpha + X->rho2)* dX->dR.b[i]);
}


static void gob_osci4(double* kosci,
		      const Gaussian*G1, const Gaussian* G2, 
		      const GaussianAux* X, const gradGaussianAux* dX, 
		      complex double* h, gradGaussian* dh)
{	
  int i;

  *h += 0.5* *kosci* (15*csqr(X->alpha)+10*X->alpha*X->rho2+csqr(X->rho2))* X->Q;

  for (i=0; i<2; i++)
    dh->chi[i] += 0.5* *kosci* X->T* dX->dS.chi[i]* 
      (15*csqr(X->alpha) + 10*X->alpha*X->rho2 + csqr(X->rho2))* X->R;

  dh->a += 0.5* *kosci* X->T* X->S* 
    (((30*X->alpha+10*X->rho2)*dX->dalpha+ 
      (10*X->alpha+2*X->rho2)*dX->drho2.a)*X->R+
     (15*csqr(X->alpha)+10*X->alpha*X->rho2+csqr(X->rho2))*dX->dR.a);
  for (i=0; i<3; i++)
    dh->b[i] += 0.5* *kosci* X->T* X->S* 
      ((10*X->alpha+2*X->rho2)*dX->drho2.b[i]* X->R+
	(15*csqr(X->alpha)+10*X->alpha*X->rho2+csqr(X->rho2))* dX->dR.b[i]);
}


void calcgradOsci2(const SlaterDet* Q, const SlaterDetAux* X, 
		   const gradSlaterDetAux* dX,
		   double omega, gradSlaterDet* dhosci)
{
  double kosci = mnucleon* omega*omega;
  gradOneBodyOperator gop_ob_osci2 = {opt: 1, par: &kosci, me: gob_osci2};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_osci2, dhosci);
}


void calcgradOsci4(const SlaterDet* Q, const SlaterDetAux* X, 
		   const gradSlaterDetAux* dX,
		   double kappa, gradSlaterDet* dhosci)
{
  double kosci = mnucleon*kappa;
  gradOneBodyOperator gop_ob_osci4 = {opt: 1, par: &kosci, me: gob_osci4};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_osci4, dhosci);
}


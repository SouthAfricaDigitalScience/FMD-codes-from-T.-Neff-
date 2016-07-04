/**

  \file Oscillator.h

  external quadratic and quartic oscillator


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"

#include "Oscillator.h"

#include "numerics/cmath.h"
#include "misc/physics.h"



static void ob_osci2(double* kosci,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* h)
{	
  *h += 0.5* *kosci* (3*X->alpha + X->rho2) * X->Q;
}


static void ob_osci4(double* kosci,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* h)
{
  *h += 0.5* *kosci* (15*csqr(X->alpha) + 10*X->alpha*X->rho2 + csqr(X->rho2))* X->Q;
}


void calcOsci2(const SlaterDet* Q, const SlaterDetAux* X, 
	       double omega, double* hosci)
{
  double kosci = mnucleon*omega*omega;

  OneBodyOperator op_ob_osci2 = {dim: 1, opt: 1, par: &kosci, me: ob_osci2};

  calcSlaterDetOBME(Q, X, &op_ob_osci2, hosci);
}


void calcOsci4(const SlaterDet* Q, const SlaterDetAux* X,
	       double kappa, double* hosci)
{
  double kosci = mnucleon*kappa;

  OneBodyOperator op_ob_osci4 = {dim: 1, opt: 1, par: &kosci, me: ob_osci4};

  calcSlaterDetOBME(Q, X, &op_ob_osci4, hosci);
}

/**

  \file NOscillator.c

  single-particle harmonic oscillator hamiltonian


  (c) 2005 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"

#include "Oscillator.h"

#include "numerics/cmath.h"
#include "misc/physics.h"

#define SQR(x) ((x)*(x))


static void ob_hosci(double* omega,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* hosci)
{
  double om = *omega;

  *hosci += (0.5/mass(G1->xi)* (3*X->lambda + X->pi2) +
	     0.5*mass(G1->xi)* SQR(om)* (3*X->alpha + X->rho2))* X->Q;
}


void calcHOsciHF(const SlaterDet* Q, const SlaterDetAux* X, double omega,
		 void* mes)
{
  OneBodyOperator op_ob_hosci = {dim: 1, opt: 1, par: &omega, me: ob_hosci};

  calcSlaterDetOBHFMEs(Q, X, &op_ob_hosci, mes);
}

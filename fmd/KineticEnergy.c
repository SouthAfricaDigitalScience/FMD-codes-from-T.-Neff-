/**

  \file KineticEnergy.c

  calculate kinetic energy


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"

#include "KineticEnergy.h"

#include "numerics/cmath.h"
#include "misc/physics.h"



static void ob_t(void* par,
		 const Gaussian* G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double* t)
{	
    *t += 0.5/mass(G1->xi)* (3*X->lambda + X->pi2) * X->Q;
}


void calcT(const SlaterDet* Q, const SlaterDetAux* X, double* t)
{
  OneBodyOperator op_ob_t = {dim: 1, opt: 1, par: NULL, me: ob_t};

  calcSlaterDetOBME(Q, X, &op_ob_t, t);
}


void calcTod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X, complex double* t)
{
  OneBodyOperator op_ob_t = {dim: 1, opt: 1, par: NULL, me: ob_t};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_t, t);
}


void calcTHF(const SlaterDet* Q, const SlaterDetAux* X, void* mes)
{
  OneBodyOperator op_ob_t = {dim: 1, opt: 1, par: NULL, me: ob_t};

  calcSlaterDetOBHFMEs(Q, X, &op_ob_t, mes);
}

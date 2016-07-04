/**

  \file gradKineticEnergy.c

  calculate gradient of kinetic energy


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "gradKineticEnergy.h"

#include "numerics/cmath.h"
#include "misc/physics.h"



static void gob_t(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, const gradGaussianAux* dX, 
		  complex double* t, gradGaussian* dt)
{	
  int i;
  double m;

  m = 0.5/mass(G1->xi);

  *t += m* (3*X->lambda + X->pi2)* X->Q;

  for (i=0; i<2; i++)
    dt->chi[i] += m* X->T* dX->dS.chi[i]* (3*X->lambda + X->pi2)* X->R;
  dt->a += m* X->T* X->S* ((3*dX->dlambda + dX->dpi2.a)* X->R +
			   (3*X->lambda + X->pi2)*dX->dR.a);
  for (i=0; i<3; i++)
    dt->b[i] += m* X->T* X->S* (dX->dpi2.b[i]* X->R +
				(3*X->lambda + X->pi2)* dX->dR.b[i]);
}



void calcgradT(const SlaterDet* Q, const SlaterDetAux* X, 
	       const gradSlaterDetAux* dX,
	       gradSlaterDet* dt)
{
  gradOneBodyOperator gop_ob_t = {opt: 1, par: NULL, me: gob_t};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_t, dt);
}


void calcgradTod(const SlaterDet* Q, const SlaterDet* Qp,
		 const SlaterDetAux* X, 
		 const gradSlaterDetAux* dX,
		 gradSlaterDet* dt)
{
  gradOneBodyOperator gop_ob_t = {opt: 1, par: NULL, me: gob_t};

  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_t, dt);
}



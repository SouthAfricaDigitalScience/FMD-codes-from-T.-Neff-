/**

  \file gradCenterofMass.c

  calculate gradient of center of mass energy


  (c) 2003 Thomas Neff

*/

#include <stdlib.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "gradCenterofMass.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


// we actually calculate the gradient of -Tcm !


static void gob_P2(SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* p2, gradGaussian* dp2)
{	
  int i;
  double M;

  M = 0.5/(Q->Z*mproton+Q->N*mneutron);

  *p2 -= M* (3*X->lambda + X->pi2) * X->Q;

  for (i=0; i<2; i++)
    dp2->chi[i] -= M* X->T* dX->dS.chi[i]* (3*X->lambda + X->pi2)* X->R;
  dp2->a -= M* X->T* X->S* ((3*dX->dlambda + dX->dpi2.a)* X->R +
			    (3*X->lambda + X->pi2)*dX->dR.a);
  for (i=0; i<3; i++)	
    dp2->b[i] -= M* X->T* X->S* (dX->dpi2.b[i]*X->R +
				 (3*X->lambda + X->pi2)* dX->dR.b[i]);
}


static void gtb_P2(SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24,
		   const gradGaussianAux* dX13,
		   complex double* p2, gradGaussian* dp2)
{	
  int i;
  double M;
  complex double pipi;

  M = 0.5/(Q->Z*mproton+Q->N*mneutron);

  pipi = cvec3mult(X13->pi, X24->pi);

  *p2 -= 2*M* pipi* X13->Q* X24->Q;

  for (i=0; i<2; i++)
    dp2->chi[i] -= 2*M* pipi* X13->T* dX13->dS.chi[i]* X13->R* X24->Q;
  dp2->a -= 2*M* (-X13->lambda*pipi* X13->R + pipi*dX13->dR.a)* 
    X13->T* X13->S* X24->Q;
  for (i=0; i<3; i++)
    dp2->b[i] -= 2*M* (dX13->dpi.b*X24->pi[i]*X13->R + pipi* dX13->dR.b[i])* 
      X13->T* X13->S* X24->Q;
}


void calcgradTCM(const SlaterDet* Q, const SlaterDetAux* X, 
		 const gradSlaterDetAux* dX,
		 gradSlaterDet* dtcm)
{
  gradOneBodyOperator gop_ob_P2 = {opt: 1, par: Q, me: gob_P2};
  gradTwoBodyOperator gop_tb_P2 = {opt: 1, par: Q, me: gtb_P2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_P2, dtcm);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_P2, dtcm);
}


void calcgradTCMod(const SlaterDet* Q, const SlaterDet* Qp, 
		   const SlaterDetAux* X, 
		   const gradSlaterDetAux* dX,
		   gradSlaterDet* dtcm)
{
  gradOneBodyOperator gop_ob_P2 = {opt: 1, par: Q, me: gob_P2};
  gradTwoBodyOperator gop_tb_P2 = {opt: 1, par: Q, me: gtb_P2};


  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_P2, dtcm);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_P2, dtcm);
}



/**

  \file ConstraintJ.c

  Constrain expectation values of Jx, Jy, Jz


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "ConstraintJ.h"

#include "numerics/cmath.h"


static void ob_j(void* par,
		 const Gaussian* G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double j[3])
{
  int i;

  for (i=0; i<3; i++)
    j[i] += X->T* ((X->rho[(i+1)%3]*X->pi[(i+2)%3]-
		    X->rho[(i+2)%3]*X->pi[(i+1)%3])* X->S* X->R +
	       0.5*X->sig[i]*X->R);
}


static void gob_j(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const GaussianAux* X, gradGaussianAux* dX,
		  complex double* jx, gradGaussian* djx)
{
  // to be rewritten

  complex double rhoxpi;
  int i;

  rhoxpi = X->rho[1]*X->pi[2] - X->rho[2]*X->pi[1];

  *jx += X->T* (rhoxpi*X->S + 0.5*X->sig[0])* X->R;

  for (i=0; i<2; i++)
    djx->chi[i] += X->T* (rhoxpi*dX->dS.chi[i] + 0.5*dX->dsig.chi[i][0])* X->R;
  djx->a += X->T* ((dX->drho.a[1]*X->pi[2] + X->rho[1]*dX->dpi.a[2] -
		    dX->drho.a[2]*X->pi[1] - X->rho[2]*dX->dpi.a[1])*X->S*X->R+
		   (rhoxpi*X->S + 0.5*X->sig[0])* dX->dR.a);
  djx->b[1] += X->T*(dX->drho.b*X->pi[2] - X->rho[2]*dX->dpi.b)*X->S*X->R;
  djx->b[2] += X->T*(X->rho[1]*dX->dpi.b - dX->drho.b*X->pi[1])*X->S*X->R;

  for (i=0; i<3; i++)
    djx->b[i] += X->T*(rhoxpi* X->S + 0.5*X->sig[0])*dX->dR.b[i]; 
}


#define SQR(x) (x)*(x)

void calcConstraintJ(const SlaterDet* Q, const SlaterDetAux* X, 
		     double* j2)
{
  double j[3];
  OneBodyOperator op_ob_j = {dim: 3, opt: 1, par: NULL, me: ob_j};

  calcSlaterDetOBME(Q, X, &op_ob_j, j);
  *j2 = SQR(j[0]) + SQR(j[1]) + SQR(j[2]);
}


void calcgradConstraintJ(const SlaterDet* Q, const SlaterDetAux* X, 
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* djx)
{
  // to be done

}

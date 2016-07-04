/**

  \file gradGaussian.c

  Derivatives of a Gaussian Wavepacket


  (c) 2003 Thomas Neff

*/

#include "Gaussian.h"
#include "gradGaussian.h"

#include <math.h>
#include <complex.h>

#include "numerics/cmath.h"



void calcgradGaussianAux(const Gaussian* G1, const Gaussian* G2, 
			 const GaussianAux* X, gradGaussianAux* dX)
{
  int i;

  dX->dsig.chi[0][0] = G2->chi[1];	dX->dsig.chi[1][0] = G2->chi[0];
  dX->dsig.chi[0][1] = -I*G2->chi[1];	dX->dsig.chi[1][1] = I*G2->chi[0];
  dX->dsig.chi[0][2] = G2->chi[0];	dX->dsig.chi[1][2] = -G2->chi[1];

  dX->dlambda = - csqr(X->lambda);
  dX->dalpha = X->lambda*(G2->a - X->alpha);

  for (i=0; i<3; i++) 
    dX->drho.a[i] = X->lambda*(G2->b[i] - X->rho[i]);
  dX->drho.b = G2->a*X->lambda;

  dX->drho2.a = 2*cvec3mult(dX->drho.a, X->rho);
  for (i=0; i<3; i++)
    dX->drho2.b[i] = 2*dX->drho.b*X->rho[i];
  
  for (i=0; i<3; i++)
    dX->dpi.a[i] = -X->lambda*X->pi[i];
  dX->dpi.b = I*X->lambda;  

  // dX->dpi2.a = 2*cvec3mult(dX->dpi.a, X->pi);
  dX->dpi2.a = -2*X->lambda*X->pi2;
  for (i=0; i<3; i++)
    dX->dpi2.b[i] = 2*dX->dpi.b*X->pi[i];
  
  for (i=0; i<2; i++) {
    dX->dS.chi[i] = G2->chi[i];
    dX->dQ.chi[i] = X->T* dX->dS.chi[i]* X->R;
  }
  dX->dR.a = 0.5*(3.0/conj(G1->a) - 3.0*X->lambda - X->pi2)*X->R;
  dX->dQ.a = X->T* X->S* dX->dR.a;
  for (i=0; i<3; i++) {
    dX->dR.b[i] = I* X->pi[i]* X->R;
    dX->dQ.b[i] = X->T* X->S* dX->dR.b[i];
  }
}

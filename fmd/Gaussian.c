/**

  \file Gaussian.c

  Parameters of a Gaussian Wavepacket


  (c) 2003 Thomas Neff

*/

#include "Gaussian.h"

#include <math.h>
#include <complex.h>

#include "misc/physics.h"
#include "numerics/cmath.h"
#include "numerics/rotationmatrices.h"


#define SQR(x) (x)*(x)

void moveGaussian(Gaussian* g, double d[3])
{
  int i;

  for (i=0; i<3; i++)
    g->b[i] += d[i];
}


void boostGaussian(Gaussian* g, double v[3])
{
  int i;

  for (i=0; i<2; i++)
    g->chi[i] *= cexp(I*mass(g->xi)*(v[0]*g->b[0]+v[1]*g->b[1]+v[2]*g->b[2])
		      -0.5*g->a*SQR(mass(g->xi))*(SQR(v[0])+SQR(v[1])+SQR(v[2])));
  for (i=0; i<3; i++)
    g->b[i] += I*mass(g->xi)*g->a*v[i];

}


void rotateGaussian(Gaussian* G, 
		    const double R3[3][3], const complex double R2[2][2])
{
  cm2mult(R2, G->chi);
  m3mult(R3, G->b);
}


void invertGaussian(Gaussian* G)
{
  int i;

  for (i=0; i<3; i++)
    G->b[i] *= -1;
}


void timerevertGaussian(Gaussian* G)
{
  int i;

  G->a = conj(G->a);
  for (i=0; i<3; i++)
    G->b[i] = conj(G->b[i]);

  complex double tmp = G->chi[0];
  G->chi[0] = -conj(G->chi[1]);
  G->chi[1] = conj(tmp);
}


void scaleGaussian(Gaussian* G, double kappa)
{
  int i;

  G->a *= kappa*kappa;	 
  for (i=0; i<3; i++)
    G->b[i] *= kappa;
}


void spinflipGaussian(Gaussian* G)
{
  complex double c;
  c = G->chi[0];
  G->chi[0] = conj(G->chi[1]);
  G->chi[1] = -conj(c);
}
    

void calcGaussianAux(const Gaussian* G1, const Gaussian* G2, GaussianAux* X)
{
  int i;

  X->sig[0] = conj(G1->chi[0])*G2->chi[1] + conj(G1->chi[1])*G2->chi[0];
  X->sig[1] = I*(conj(G1->chi[1])*G2->chi[0] - conj(G1->chi[0])*G2->chi[1]);
  X->sig[2] = conj(G1->chi[0])*G2->chi[0] - conj(G1->chi[1])*G2->chi[1];

  X->lambda = 1.0/(conj(G1->a)+G2->a);
  X->alpha = conj(G1->a)*G2->a*X->lambda;
  for (i=0; i<3; i++)
    X->rho[i] = X->lambda*(G2->a*conj(G1->b[i])+conj(G1->a)*G2->b[i]);
  X->rho2 = cvec3sqr(X->rho);
  for (i=0; i<3; i++)	
    X->pi[i] = X->lambda*I*(conj(G1->b[i]) - G2->b[i]);
  X->pi2 = cvec3sqr(X->pi);
  X->rhopi = cvec3mult(X->rho, X->pi);
  cvec3cross(X->rho, X->pi, X->rhoxpi);

  X->T = (1+G1->xi*G2->xi)/2;
  X->S = conj(G1->chi[0])*G2->chi[0] + conj(G1->chi[1])*G2->chi[1];
  X->R = cpow32(2*M_PI*X->alpha)*cexp(0.5*X->pi2/X->lambda);
  X->Q = X->T*X->S*X->R;
}

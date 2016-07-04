/**

  \file SpinDensities.c

  calculate one-body spin densities in coordinate and momentum space


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"

#include "SpinDensities.h"

#include "numerics/cmath.h"


typedef struct {
  int normal;
  int npoints;
  double range;
} densparameters;


static void ob_spindensx(densparameters* par,
			 const Gaussian* G1, const Gaussian* G2, 
			 const GaussianAux* X, complex double* dens)
{	
  int v = par->normal;
  int n = par->npoints;
  double xmax = par->range;
  int i,j,k;
  double x[3];
  int vx, vy;
  complex double rho;

  if (v==0) {
    vx = 1; vy = 2;
  } else if (v==1) {
    vx = 0; vy = 2;
  } else {
    vx = 0; vy = 1;
  }

  for (j=0; j<n; j++)
    for (i=0; i<n; i++) {

      x[v] = 0.0;
      x[vx] = -xmax+2*xmax*i/(n-1);
      x[vy] = -xmax+2*xmax*j/(n-1);

      rho =
	cexp(-0.5*(csqr(x[0] - conj(G1->b[0]))+csqr(x[1] - conj(G1->b[1]))+
		   csqr(x[2] - conj(G1->b[2])))/conj(G1->a)
	     -0.5*(csqr(x[0] - G2->b[0])+csqr(x[1] - G2->b[1])+
		   csqr(x[2] - G2->b[2]))/G2->a)*X->T;
      for (k=0; k<3; k++)
	dens[k+(i+j*n)*3] += 0.5*rho*X->sig[k];
    }
}


static void ob_spindensp(densparameters* par,
			 const Gaussian* G1, const Gaussian* G2, 
			 const GaussianAux* X, complex double* dens)
{	
  int v = par->normal;
  int n = par->npoints;
  double pmax = par->range;
  int i,j,k;
  double p[3], p2;
  int vx, vy;
  complex double rho;

  if (v==0) {
    vx = 1; vy = 2;
  } else if (v==1) {
    vx = 0; vy = 2;
  } else {
    vx = 0; vy = 1;
  }

  for (j=0; j<n; j++)
    for (i=0; i<n; i++) {

      p[v] = 0.0;
      p[vx] = -pmax+2*pmax*i/(n-1);
      p[vy] = -pmax+2*pmax*j/(n-1);
	
      p2 = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];

      rho= 
	cpow32(conj(G1->a)*G2->a)*
	cexp(-0.5*conj(G1->a)*p2 + 
	     I*(p[0]*conj(G1->b[0])+p[1]*conj(G1->b[1])+p[2]*conj(G1->b[2]))
	     -0.5*G2->a*p2 
	     -I*(p[0]*G2->b[0]+p[1]*G2->b[1]+p[2]*G2->b[2]))*X->T;
      
      for (k=0; k<3; k++)
	dens[k+(i+j*n)*3] += 0.5*rho*X->sig[k];
    }
}



void calcSpinDensitiesCoordinate(const SlaterDet* Q, const SlaterDetAux* X, 
				 int v, int n, double x, double* dens)
{
  densparameters P = { normal: v, npoints: n, range: x };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_spindensx};

  calcSlaterDetOBME(Q, X, &op_ob_dens, dens);
}


void calcSpinDensitiesMomentum(const SlaterDet* Q, const SlaterDetAux* X, 
			       int v, int n, double p, double* dens)
{
  densparameters P = { normal: v, npoints: n, range: p };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_spindensp};

  calcSlaterDetOBME(Q, X, &op_ob_dens, dens);
}



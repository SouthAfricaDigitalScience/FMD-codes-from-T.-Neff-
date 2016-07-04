/**

  \file Densities3d.c

  calculate one-body densities in coordinate and momentum space


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"

#include "Densities3d.h"

#include "numerics/cmath.h"


typedef struct {
  int npoints;
  double range;
} densparameters;


static void ob_densx(densparameters* par,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* dens)
{	
  int n = par->npoints;
  double xmax = par->range;
  int i,j,k;
  double x,y,z;

  for (k=0; k<n; k++) 
    for (j=0; j<n; j++)
      for (i=0; i<n; i++) {

	x = -xmax+2*xmax*i/(n-1);
	y = -xmax+2*xmax*j/(n-1);
	z = -xmax+2*xmax*k/(n-1);

	dens[i+j*n+k*n*n] += 
	  cexp(-0.5*(csqr(x - conj(G1->b[0]))+csqr(y - conj(G1->b[1]))+
		     csqr(z - conj(G1->b[2])))/conj(G1->a)
	       -0.5*(csqr(x - G2->b[0])+csqr(y - G2->b[1])+
		     csqr(z - G2->b[2]))/G2->a)*
	  X->T*X->S;
      }
}


static void ob_densp(densparameters* par,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* dens)
{	
  int n = par->npoints;
  double pmax = par->range;
  int i,j,k;
  double px,py,pz,p2;

  for (k=0; k<n; k++) 
    for (j=0; j<n; j++)
      for (i=0; i<n; i++) {

	px = -pmax+2*pmax*i/(n-1);
	py = -pmax+2*pmax*j/(n-1);
	pz = -pmax+2*pmax*k/(n-1);
	
	p2 = px*px+py*py+pz*pz;

	dens[i+j*n+k*n*n] += 
	  cpow32(conj(G1->a)*G2->a)*
	  cexp(-0.5*conj(G1->a)*p2 + 
	       I*(px*conj(G1->b[0])+py*conj(G1->b[1])+pz*conj(G1->b[2]))
	       -0.5*G2->a*p2 
	       -I*(px*G2->b[0]+py*G2->b[1]+pz*G2->b[2]))*
	  X->T*X->S;
      }
}



void calcDensitiesCoordinate3d(const SlaterDet* Q, const SlaterDetAux* X, 
			       int n, double x, double* dens)
{
  densparameters P = { npoints: n, range: x };
  OneBodyOperator op_ob_dens = {dim: n*n*n, opt: 1, par: &P, me: ob_densx};

  calcSlaterDetOBME(Q, X, &op_ob_dens, dens);
}


void calcDensitiesMomentum3d(const SlaterDet* Q, const SlaterDetAux* X, 
			     int n, double p, double* dens)
{
  densparameters P = { npoints: n, range: p };
  OneBodyOperator op_ob_dens = {dim: n*n*n, opt: 1, par: &P, me: ob_densp};

  calcSlaterDetOBME(Q, X, &op_ob_dens, dens);
}



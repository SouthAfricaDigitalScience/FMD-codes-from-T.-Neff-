/**

  \file Densities.c

  calculate one-body densities in coordinate and momentum space


  (c) 2003 Thomas Neff

*/

#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "SlaterDet.h"

#include "Densities.h"

#include "numerics/cmath.h"
#include "numerics/zcw.h"


typedef struct {
  int normal;
  int npoints;
  double range;
} densparameters;


typedef struct {
  int izcw;
  int npoints;
  double range;
} rdensparameters;


static void ob_densx(densparameters* par,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* dens)
{	
  int v = par->normal;
  int n = par->npoints;
  double xmax = par->range;
  int i,j;
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
		   csqr(x[2] - G2->b[2]))/G2->a)* X->S;

      dens[i+j*n      ] += rho* X->T;
      if (G1->xi == +1)
	dens[i+j*n+  n*n] += rho* X->T;
      if (G1->xi == -1)
	dens[i+j*n+2*n*n] += rho* X->T;
    }
}


static void ob_rdensx(rdensparameters* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const GaussianAux* X, complex double* dens)
{	
  int izcw = par->izcw;
  int n = par->npoints;
  double rmax = par->range;
  int i,j;
  double r, alpha, beta, x[3];
  complex double rho;

  int nzcw = nangles2(izcw);
  double weight = 1.0/nzcw;

  for (i=0; i<n; i++) {
    r = rmax* i/(n-1);

    for (j=0; j<nzcw; j++) {
      getangles2(izcw, j, &alpha, &beta);

      x[0] = r*cos(alpha)*sin(beta);
      x[1] = r*sin(alpha)*sin(beta);
      x[2] = r*cos(beta);

      rho =
	cexp(-0.5*(csqr(x[0] - conj(G1->b[0]))+csqr(x[1] - conj(G1->b[1]))+
		   csqr(x[2] - conj(G1->b[2])))/conj(G1->a)
	     -0.5*(csqr(x[0] - G2->b[0])+csqr(x[1] - G2->b[1])+
		   csqr(x[2] - G2->b[2]))/G2->a)* X->S;

      dens[i   ] += weight* rho* X->T;
      if (G1->xi == 1)
	dens[i+  n] += weight* rho* X->T;
      if (G1->xi == -1)
	dens[i+2*n] += weight* rho* X->T;
    }
  }	
}


static void ob_densp(densparameters* par,
		     const Gaussian* G1, const Gaussian* G2, 
		     const GaussianAux* X, complex double* dens)
{	
  int v = par->normal;
  int n = par->npoints;
  double pmax = par->range;
  int i,j;
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

      rho = 
	cpow32(conj(G1->a)*G2->a)*
	cexp(-0.5*conj(G1->a)*p2 + 
	     I*(p[0]*conj(G1->b[0])+p[1]*conj(G1->b[1])+p[2]*conj(G1->b[2]))
	     -0.5*G2->a*p2 
	     -I*(p[0]*G2->b[0]+p[1]*G2->b[1]+p[2]*G2->b[2]))* X->S;
      
      dens[i+j*n      ] += rho* X->T;
      if (G1->xi == +1)
	dens[i+j*n+  n*n] += rho* X->T;
      if (G1->xi == -1)
	dens[i+j*n+2*n*n] += rho* X->T;
    }
}


static void ob_rdensp(rdensparameters* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const GaussianAux* X, complex double* dens)
{	
  int izcw = par->izcw;
  int n = par->npoints;
  double pmax = par->range;
  int i,j;
  double pr, alpha, beta, p[3], p2;
  complex double rho;

  int nzcw = nangles2(izcw);
  double weight = 1.0/nzcw;

  for (i=0; i<n; i++) {
    pr = pmax* i/(n-1);

    for (j=0; j<nzcw; j++) {
      getangles2(izcw, j, &alpha, &beta);

      p[0] = pr*cos(alpha)*sin(beta);
      p[1] = pr*sin(alpha)*sin(beta);
      p[2] = pr*cos(beta);

      p2 = pr*pr;

      rho = 
	cpow32(conj(G1->a)*G2->a)*
	cexp(-0.5*conj(G1->a)*p2 + 
	     I*(p[0]*conj(G1->b[0])+p[1]*conj(G1->b[1])+p[2]*conj(G1->b[2]))
	     -0.5*G2->a*p2 
	     -I*(p[0]*G2->b[0]+p[1]*G2->b[1]+p[2]*G2->b[2]))* X->S;

      dens[i   ] += weight* rho* X->T;
      if (G1->xi == 1)
	dens[i+  n] += weight* rho* X->T;
      if (G1->xi == -1)
	dens[i+2*n] += weight* rho* X->T;
    }
  }
}


void calcDensitiesCoordinate(const SlaterDet* Q,
			     int v, int n, double x, double* dens)
{
  SlaterDetAux X;
  initSlaterDetAux(Q, &X);
  calcSlaterDetAux(Q, &X);

  densparameters P = { normal: v, npoints: n, range: x };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_densx};

  calcSlaterDetOBME(Q, &X, &op_ob_dens, dens);

  freeSlaterDetAux(&X);
}


void calcDensitiesCoordinateParity(const SlaterDet* Q, int par,
				   int v, int n, double x, double* dens)
{
  double nlr;
  complex double ovllr;
  complex double denslr[3*n*n];

  SlaterDet Ql, Qr;
  initSlaterDet(Q, &Ql);
  initSlaterDet(Q, &Qr);

  SlaterDetAux X;
  initSlaterDetAux(Q, &X);

  densparameters P = { normal: v, npoints: n, range: x };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_densx};

  int ipl, ipr, k;

  nlr = 0;
  copySlaterDet(Q, &Ql);
  for (ipl=0; ipl<=1; ipl++) {

    if (ipl)
      invertSlaterDet(&Ql);

    copySlaterDet(Q, &Qr);
    for (ipr=0; ipr<=1; ipr++) {

      if (ipr)
        invertSlaterDet(&Qr);

      calcSlaterDetAuxod(&Ql, &Qr, &X);
      ovllr = X.ovlap;

      calcSlaterDetOBMEod(&Ql, &Qr, &X, &op_ob_dens, denslr);

      nlr += ((par == -1 && (ipl+ipr)%2) ? -1 : 1)* ovllr;

      for (k=0; k<3*n*n; k++)
        dens[k] += ((par == -1 && (ipl+ipr)%2) ? -1 : 1)* denslr[k];

    }

  }

  for (k=0; k<3*n*n; k++)
    dens[k] = dens[k]/nlr;

  freeSlaterDet(&Ql); freeSlaterDet(&Qr);
  freeSlaterDetAux(&X);
}


void calcRadialDensitiesCoordinate(const SlaterDet* Q,
				   int izcw, int n, double r, double* dens)
{
  SlaterDetAux X;
  initSlaterDetAux(Q, &X);
  calcSlaterDetAux(Q, &X);

  rdensparameters P = { izcw: izcw, npoints: n, range: r };
  OneBodyOperator op_ob_rdens = {dim: 3*n, opt: 1, par: &P, me: ob_rdensx};

  calcSlaterDetOBME(Q, &X, &op_ob_rdens, dens);

  freeSlaterDetAux(&X);
}


void calcDensitiesMomentum(const SlaterDet* Q,
			   int v, int n, double p, double* dens)
{
  SlaterDetAux X;
  initSlaterDetAux(Q, &X);
  calcSlaterDetAux(Q, &X);

  densparameters P = { normal: v, npoints: n, range: p };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_densp};

  calcSlaterDetOBME(Q, &X, &op_ob_dens, dens);

  freeSlaterDetAux(&X);
}


void calcRadialDensitiesMomentum(const SlaterDet* Q,
				 int izcw, int n, double p, double* dens)
{
  SlaterDetAux X;
  initSlaterDetAux(Q, &X);
  calcSlaterDetAux(Q, &X);

  rdensparameters P = { izcw: izcw, npoints: n, range: p };
  OneBodyOperator op_ob_rdens = {dim: 3*n, opt: 1, par: &P, me: ob_rdensp};

  calcSlaterDetOBME(Q, &X, &op_ob_rdens, dens);

  freeSlaterDetAux(&X);
}


void calcDensitiesCoordinateHF(const SlaterDet* Q, const SlaterDetAux* X, 
                               const int v, const int n, const double x, 
                               void* mes)
{
  densparameters P = { normal: v, npoints: n, range: x };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_densx};
	
  calcSlaterDetOBHFMEs(Q, X, &op_ob_dens, mes);
}


void calcDensitiesMomentumHF(const SlaterDet* Q, const SlaterDetAux* X, 
                             const int v, const int n, const double x, 
                             void* mes)
{
  densparameters P = { normal: v, npoints: n, range: x };
  OneBodyOperator op_ob_dens = {dim: 3*n*n, opt: 1, par: &P, me: ob_densp};
	
  calcSlaterDetOBHFMEs(Q, X, &op_ob_dens, mes);
}

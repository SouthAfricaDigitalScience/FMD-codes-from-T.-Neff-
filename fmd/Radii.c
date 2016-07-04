/**

  \file Radii.c

  calculate radii


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Radii.h"

#include "numerics/cmath.h"


static void ob_radii2(const SlaterDet* Q,
	  const Gaussian*G1, const Gaussian* G2, 
	  const GaussianAux* X, complex double r2[3])
{
  complex double mr2;
  int A=Q->A, Z=Q->Z, N=Q->N;

  mr2 = (3.0*X->alpha + X->rho2)* X->Q;
 
  r2[0] += 1.0/A*(1.0-1.0/A)* mr2;
  r2[1] += (1.0/(A*A) + (1+G1->xi)/2 *(1.0/Z-2.0/(A*Z)))* mr2;
  r2[2] += (1.0/(A*A) + (1-G1->xi)/2 *(1.0/N-2.0/(A*N)))* mr2;
}


static void tb_radii2(const SlaterDet* Q,
	   const Gaussian* G1, const Gaussian* G2, 
	   const Gaussian* G3, const Gaussian* G4, 
	   const GaussianAux* X13, const GaussianAux* X24, 
	   complex double r2[3])
{	
  complex double mr2;
  int A=Q->A, Z=Q->Z, N=Q->N;

  mr2 = cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;

  r2[0] += -2.0/(A*A)* mr2;
  r2[1] +=  2.0*(1.0/(A*A) - 
		 ((1+G1->xi)/2+(1+G2->xi)/2) *1.0/(A*Z))* mr2;
  r2[2] +=  2.0*(1.0/(A*A) - 
		 ((1-G1->xi)/2+(1-G2->xi)/2) *1.0/(A*N))* mr2;
}


void calcRadii2(const SlaterDet* Q, const SlaterDetAux* X, 
		double* r2mass, double* r2proton, double* r2neutron)
{
  double r2one[3], r2two[3];
  OneBodyOperator op_ob_radii2 = {dim: 3, opt: 1, par: Q, me: ob_radii2};
  TwoBodyOperator op_tb_radii2 = {dim: 3, opt: 1, par: Q, me: tb_radii2};


  calcSlaterDetOBME(Q, X, &op_ob_radii2, r2one);
  calcSlaterDetTBME(Q, X, &op_tb_radii2, r2two);

  *r2mass = r2one[0] + r2two[0];
  *r2proton = r2one[1] + r2two[1];
  *r2neutron = r2one[2] + r2two[2];
}


void calcRadii2od(const SlaterDet* Q, const SlaterDet* Qp,
		  const SlaterDetAux* X, 
		  complex double* r2mass, 
		  complex double* r2proton, complex double* r2neutron)
{
  complex double r2one[3], r2two[3];
  OneBodyOperator op_ob_radii2 = {dim: 3, opt: 1, par: Q, me: ob_radii2};
  TwoBodyOperator op_tb_radii2 = {dim: 3, opt: 1, par: Q, me: tb_radii2};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_radii2, r2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_radii2, r2two);

  *r2mass = r2one[0] + r2two[0];
  *r2proton = r2one[1] + r2two[1];
  *r2neutron = r2one[2] + r2two[2];
}


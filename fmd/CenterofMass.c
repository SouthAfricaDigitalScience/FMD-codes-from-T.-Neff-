/**

  \file CenterofMass.c

  calculate observables related to center of mass motion


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "CenterofMass.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


#define SQR(x)	((x)*(x))

static void ob_mx(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double x[3])
{
  int i;
  
  for (i=0; i<3; i++)
    x[i] += mass(G1->xi)*X->rho[i]* X->Q;
}


static void ob_p(void* par,
		 const Gaussian*G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double p[3])
{	
  int i;
  
  for (i=0; i<3; i++)
    p[i] += X->pi[i] * X->Q;
}


static void ob_tcm(SlaterDet* Q,
		   const Gaussian*G1, const Gaussian* G2, 
		   const GaussianAux* X, complex double* p2)
{	
  double M;
  M = 0.5/(Q->Z*mproton+Q->N*mneutron);
   
  *p2 += M* (3*X->lambda + X->pi2) * X->Q;
}


static void tb_tcm(SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24, 
		   complex double* p2)
{	
  double M;
  M = 0.5/(Q->Z*mproton+Q->N*mneutron);

  *p2 += M* 2.0*cvec3mult(X13->pi, X24->pi)* X13->Q* X24->Q;
}



void calcCMPosition(const SlaterDet* Q, const SlaterDetAux* X, double xcm[3])
{
  int i;
  double M;
  OneBodyOperator op_ob_mx = {dim: 3, opt: 1, par: NULL, me: ob_mx};

  calcSlaterDetOBME(Q, X, &op_ob_mx, xcm);

  M = Q->Z*mproton+Q->N*mneutron;
  for (i=0; i<3; i++)
    xcm[i] /= M;
}


void calcCMVelocity(const SlaterDet* Q, const SlaterDetAux* X, double vcm[3])
{
  int i;
  double M;
  OneBodyOperator op_ob_p = {dim: 3, opt: 1, par: NULL, me: ob_p};

  calcSlaterDetOBME(Q, X, &op_ob_p, vcm);

  M = Q->Z*mproton+Q->N*mneutron;
  for (i=0; i<3; i++)
    vcm[i] /= M;
}


void calcTCM(const SlaterDet* Q, const SlaterDetAux* X, double* tcm)
{
  double tcmone, tcmtwo;

  OneBodyOperator op_ob_tcm = {dim: 1, opt: 1, par: Q, me: ob_tcm};
  TwoBodyOperator op_tb_tcm = {dim: 1, opt: 1, par: Q, me: tb_tcm};

  calcSlaterDetOBME(Q, X, &op_ob_tcm, &tcmone);
  calcSlaterDetTBME(Q, X, &op_tb_tcm, &tcmtwo);

  *tcm = tcmone+tcmtwo;
}


void calcTCMod(const SlaterDet* Q, const SlaterDet* Qp,
	       const SlaterDetAux* X, complex double* tcm)
{
  complex double tcmone, tcmtwo;

  OneBodyOperator op_ob_tcm = {dim: 1, opt: 1, par: Q, me: ob_tcm};
  TwoBodyOperator op_tb_tcm = {dim: 1, opt: 1, par: Q, me: tb_tcm};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_tcm, &tcmone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_tcm, &tcmtwo);

  *tcm = tcmone+tcmtwo;
}


void calcTCMHF(const SlaterDet* Q, const SlaterDetAux* X, void* mes)
{
  OneBodyOperator op_ob_tcm = {dim: 1, opt: 1, par: Q, me: ob_tcm};
  TwoBodyOperator op_tb_tcm = {dim: 1, opt: 1, par: Q, me: tb_tcm};

  int A=Q->A;
  complex double tcmone[A*A];
  complex double tcmtwo[A*A];

  calcSlaterDetOBHFMEs(Q, X, &op_ob_tcm, tcmone);
  calcSlaterDetTBHFMEs(Q, X, &op_tb_tcm, tcmtwo);

  complex double *tcm = mes;
  int k,l;
  
  for (l=0; l<A; l++)
    for (k=0; k<A; k++)
      tcm[k+l*A] = tcmone[k+l*A] + tcmtwo[k+l*A];
}

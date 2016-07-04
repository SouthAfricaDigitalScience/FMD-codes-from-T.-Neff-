/**

  \file gradCMSlaterDet.c

  gradients of matrix elements that that depend
  on the expectation values of Center of Mass coordinate and momentum


  (c) 2003 Thomas Neff

*/



#include "Gaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "gradCMSlaterDet.h"

#include "misc/physics.h"


typedef struct {
  int A; int Z; int N;
  double X[3], V[3];
} gcmpara;


void gob_dcm(const gcmpara* para,
	     const Gaussian* G1, const Gaussian* G2,
	     const GaussianAux* X, const gradGaussianAux* dX,
	     complex double* dcm, gradGaussian* ddcm)
{
  int i;
  double M;
  double* dXcm = para->X;
  double* dVcm = para->V;
  complex double meXVcm = 0.0;
  
  M=para->Z*mproton+para->N*mneutron;

  for (i=0; i<3; i++) {
    meXVcm += dXcm[i]*mass(G1->xi)/M* X->rho[i]; 
    meXVcm += dVcm[i]/M* X->pi[i]; 
  }

  *dcm += meXVcm* X->Q;

  for (i=0; i<2; i++)
    ddcm->chi[i] += meXVcm* dX->dQ.chi[i];
  for (i=0; i<3; i++)
    ddcm->a += (dXcm[i]*mass(G1->xi)/M* dX->drho.a[i] +
	       dVcm[i]/M* dX->dpi.a[i])* X->Q;
  ddcm->a += meXVcm* dX->dQ.a;
  for (i=0; i<3; i++)
    ddcm->b[i] += (dXcm[i]*mass(G1->xi)/M* dX->drho.b +
		  dVcm[i]/M *dX->dpi.b)* X->Q + meXVcm* dX->dQ.b[i];
}


void calcgradCMSlaterDetOBME(const SlaterDet* Q, const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const OneBodyOperator* op,
			     gradSlaterDet* grad)
{
  int i;
  double dcm[6];
  calcSlaterDetOBME(Q, X, op, dcm);

  gcmpara para;
  para.A = Q->A; para.Z = Q->Z; para.N = Q->N;
  for (i=0; i<3; i++)
    para.X[i] = dcm[i];
  for (i=0; i<3; i++)
    para.V[i] = dcm[i+3];

  gradOneBodyOperator gop_ob_dcm = {opt: op->opt, par: &para, me: gob_dcm}; 

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_dcm, grad);
}


void calcgradCMSlaterDetTBME(const SlaterDet* Q, const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const TwoBodyOperator* op,
			     gradSlaterDet* grad)
{
  int i;
  double dcm[6];
  calcSlaterDetTBME(Q, X, op, dcm);

  gcmpara para;
  para.A = Q->A; para.Z = Q->Z; para.N = Q->N;
  for (i=0; i<3; i++)
    para.X[i] = dcm[i];
  for (i=0; i<3; i++)
    para.V[i] = dcm[i+3];

  gradOneBodyOperator gop_ob_dcm = {opt: op->opt, par: &para, me: gob_dcm}; 

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_dcm, grad);
}

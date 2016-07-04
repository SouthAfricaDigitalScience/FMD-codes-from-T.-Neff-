/**

  \file ConstraintMainAxes.c

  Constrain Main Axes of Quadrupole tensor


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "ConstraintMainAxes.h"

#include "numerics/cmath.h"


typedef struct {
  int A;
  double q[3];
} quadrupolepara;


static void ob_qdiag(const SlaterDet* Q,
		     const Gaussian*G1, const Gaussian* G2, 
		     const GaussianAux* X, 
		     complex double q[3])
{
  int A=Q->A;
  int i;

  for (i=0; i<3; i++)
    q[i] += 3*(1.0-1.0/A)* (X->rho[i]*X->rho[i]-X->rho2)* X->Q;
}


static void ob_qoffdiag(const SlaterDet* Q,
			const Gaussian*G1, const Gaussian* G2, 
			const GaussianAux* X, 
			complex double q[3])
{
  int A=Q->A;
  int i;

  for (i=0; i<3; i++)
    q[i] += 3*(1.0-1.0/A)*X->rho[(i+1)%3]*X->rho[(i+2)%3]*X->Q;
}


static void tb_qdiag(const SlaterDet* Q,
		     const Gaussian* G1, const Gaussian* G2, 
		     const Gaussian* G3, const Gaussian* G4, 
		     const GaussianAux* X13, const GaussianAux* X24, 
		     complex double q[3])
{
  int A=Q->A;
  int i;

  for (i=0; i<3; i++)
    q[i] += (-6.0/A* X13->rho[i]*X24->rho[i] +
	      2.0/A* cvec3mult(X13->rho, X24->rho))* X13->Q* X24->Q;
}


static void tb_qoffdiag(const SlaterDet* Q,
			const Gaussian* G1, const Gaussian* G2, 
			const Gaussian* G3, const Gaussian* G4, 
			const GaussianAux* X13, const GaussianAux* X24, 
			complex double q[3])
{
  int A=Q->A;
  int i;

  for (i=0; i<3; i++)
    q[i] += -3.0/A*(X13->rho[(i+1)%3]*X24->rho[(i+2)%3]+
		    X13->rho[(i+2)%3]*X24->rho[(i+1)%3])* X13->Q* X24->Q;
}


static void gob_mainaxes(quadrupolepara* qpar,
			 const Gaussian* G1, const Gaussian* G2,
			 const GaussianAux* X, const gradGaussianAux* dX,
			 complex double* q2, gradGaussian* dq2)
{
  int A=qpar->A;
  double* q=qpar->q;
  int i, k;
  complex double merr;

  for (i=0; i<3; i++) {
    
    merr = 3*(1.0-1.0/A)*X->rho[(i+1)%3]*X->rho[(i+2)%3];
      
    *q2 += 2*q[i]* merr *X->Q;

    for (k=0; k<2; k++)
      dq2->chi[k] +=  2*q[i]* merr* dX->dQ.chi[k];
    dq2->a += 2*q[i]*
      (3*(1.0-1.0/A)*(dX->drho.a[(i+1)%3]*X->rho[(i+2)%3]+
		      X->rho[(i+1)%3]*dX->drho.a[(i+2)%3])*X->Q+
       merr* dX->dQ.a);
    dq2->b[(i+1)%3] += 2*q[i]*
      3*(1.0-1.0/A)*dX->drho.b*X->rho[(i+2)%3]*X->Q;
    dq2->b[(i+2)%3] += 2*q[i]*
      3*(1.0-1.0/A)*X->rho[(i+1)%3]*dX->drho.b*X->Q;
    for (k=0; k<3; k++) 
      dq2->b[k] += 2*q[i]* merr* dX->dQ.b[k];
  }
}


static void gtb_mainaxes(quadrupolepara* qpar,
			 const Gaussian* G1, const Gaussian* G2,
			 const Gaussian* G3, const Gaussian* G4,
			 const GaussianAux* X13, const GaussianAux* X24,
			 const gradGaussianAux* dX13,
			 complex double* q2, gradGaussian* dq2)
{
  int A=qpar->A;
  double* q=qpar->q;
  int i, k;
  complex double merr;

  for (i=0; i<3; i++) {

    merr = -3.0/A* (X13->rho[(i+1)%3]*X24->rho[(i+2)%3]+
		    X13->rho[(i+2)%3]*X24->rho[(i+1)%3]);
      
    *q2 += 2*q[i]* merr* X13->Q* X24->Q;

    for (k=0; k<2; k++)
      dq2->chi[k] +=  2*q[i]* merr* dX13->dQ.chi[k]* X24->Q;
    dq2->a += 2*q[i]*
      (-3.0/A*(dX13->drho.a[(i+1)%3]*X24->rho[(i+2)%3] +
	       dX13->drho.a[(i+2)%3]*X24->rho[(i+1)%3])* X13->Q* X24->Q + 
       merr* dX13->dQ.a* X24->Q);
    dq2->b[(i+1)%3] += 2*q[i]*
      (-3.0/A)*(dX13->drho.b* X24->rho[(i+2)%3] +
		dX13->drho.b* X24->rho[(i+3)%3])* X13->Q* X24->Q;
    for (k=0; k<3; k++) 
	dq2->b[k] += 2*q[i]* merr* dX13->dQ.b[k]* X24->Q;
  }
}


static void gob_qbeta2(quadrupolepara* qpar,
		       const Gaussian* G1, const Gaussian* G2,
		       const GaussianAux* X, const gradGaussianAux* dX,
		       complex double* beta2, gradGaussian* dbeta2)
{
  int A=qpar->A;
  double* q=qpar->q;
  int i, k;
  complex double merr;

  for (i=0; i<3; i++) {
    
    merr = (1.0-1.0/A)*(3*X->rho[i]*X->rho[i] - X->rho2);
      
    *beta2 += 2*q[i]* merr *X->Q;

    for (k=0; k<2; k++)
      dq2->chi[k] +=  2*q[i]* merr* dX->dQ.chi[k];
    dq2->a += 2*q[i]*
      ((1.0-1.0/A)*(3*dX->drho.a[i]*X->rho[i] - dX->drho2.a)*X->Q+
       merr* dX->dQ.a);
    dq2->b[i] += 2*q[i]*
      (1.0-1.0/A)*(6*dX->drho.b*X->rho[i] - dX->drho2.b[i])*X->Q;
    for (k=0; k<3; k++) 
      dq2->b[k] += 2*q[i]* merr* dX->dQ.b[k];
  }
}


static void gtb_beta2(quadrupolepara* qpar,
		      const Gaussian* G1, const Gaussian* G2,
		      const Gaussian* G3, const Gaussian* G4,
		      const GaussianAux* X13, const GaussianAux* X24,
		      const gradGaussianAux* dX13,
		      complex double* beta2, gradGaussian* dbeta2)
{
  int A=qpar->A;
  double* q=qpar->q;
  int i, k;
  complex double merr;

  for (i=0; i<3; i++) {

    merr = -6.0/A* X13->rho[i]*X24->rho[i] + 
            2.0/A* cvec3mult(X13->rho, X24->rho);
      
    *beta2 += 2*q[i]* merr* X13->Q* X24->Q;

    for (k=0; k<2; k++)
      dq2->chi[k] +=  2*q[i]* merr* dX13->dQ.chi[k]* X24->Q;
    dq2->a += 2*q[i]*
      (-6.0/A*dX13->drho.a[i]*X24->rho[i] +
        2.0/A*cvec3mult(dX13->drho.a, X24->rho))* X13->Q* X24->Q + 
       merr* dX13->dQ.a* X24->Q);
    dq2->b[(i+1)%3] += 2*q[i]*
      (-3.0/A)*(dX13->drho.b* X24->rho[(i+2)%3] +
		dX13->drho.b* X24->rho[(i+3)%3])* X13->Q* X24->Q;
    for (k=0; k<3; k++) 
	dq2->b[k] += 2*q[i]* merr* dX13->dQ.b[k]* X24->Q;
  }
}



#define SQR(x) (x)*(x)

void calcConstraintMainAxes(const SlaterDet* Q, const SlaterDetAux* X, 
			    double* q2)
{
  double qijone[3], qijtwo[3], qsq = 0.0;
  int i;
  OneBodyOperator op_ob_qij = {dim: 3, opt: 1, par: Q, me: ob_qoffdiag};
  TwoBodyOperator op_tb_qij = {dim: 3, opt: 1, par: Q, me: tb_qoffdiag};


  calcSlaterDetOBME(Q, X, &op_ob_qij, qijone);
  calcSlaterDetTBME(Q, X, &op_tb_qij, qijtwo);

  for (i=0; i<3; i++)
    qsq += SQR(qijone[i]+qijtwo[i]);

  *q2 = qsq;
}


void calcgradConstraintMainAxes(const SlaterDet* Q, const SlaterDetAux* X,
				const gradSlaterDetAux* dX,
				gradSlaterDet* dq2)
{
  double qijone[3], qijtwo[3];
  int i;
  quadrupolepara qpar;

  OneBodyOperator op_ob_qij = {dim: 3, opt: 1, par: Q, me: ob_offdiag};
  TwoBodyOperator op_tb_qij = {dim: 3, opt: 1, par: Q, me: tb_offdiag};


  calcSlaterDetOBME(Q, X, &op_ob_qij, qijone);
  calcSlaterDetTBME(Q, X, &op_tb_qij, qijtwo);

  q2par.A = Q->A;
  for (i=0; i<3; i++)
    q2par.q[i] = qijone[i]+qijtwo[i];

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_mainaxes};
  gradTwoBodyOperator gop_tb_q2 = {opt: 1, par: &q2par, me: gtb_mainaxes};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_q2, dq2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_q2, dq2);
}


void calcConstraintQbeta2(const SlaterDet* Q, const SlaterDetAux* X, 
			  double* beta2)
{
  double qone[3], qtwo[3], qsq = 0.0;
  int i;
  OneBodyOperator op_ob_q = {dim: 3, opt: 1, par: Q, me: ob_qdiag};
  TwoBodyOperator op_tb_q = {dim: 3, opt: 1, par: Q, me: tb_qdiag};


  calcSlaterDetOBME(Q, X, &op_ob_q, qone);
  calcSlaterDetTBME(Q, X, &op_tb_q, qtwo);

  for (i=0; i<3; i++)
    qsq += SQR(qone[i]+qtwo[i]);

  *beta = qsq;
}


void calcgradConstraintQbeta2(const SlaterDet* Q, const SlaterDetAux* X,
			      const gradSlaterDetAux* dX,
			      gradSlaterDet* dbeta2)
{
  double qone[3], qtwo[3];
  int i;
  quadrupolepara qpara;

  OneBodyOperator op_ob_q = {dim: 3, opt: 1, par: Q, me: ob_offdiag};
  TwoBodyOperator op_tb_q = {dim: 3, opt: 1, par: Q, me: tb_offdiag};


  calcSlaterDetOBME(Q, X, &op_ob_q, qone);
  calcSlaterDetTBME(Q, X, &op_tb_q, qtwo);

  qpara.A = Q->A;
  for (i=0; i<3; i++)
    qpara.q[i] = qone[i]+qtwo[i];

  gradOneBodyOperator gop_ob_qbeta2 = {opt: 1, par: &qpara, me: gob_qbeta2};
  gradTwoBodyOperator gop_tb_qbeta2 = {opt: 1, par: &qpara, me: gtb_qbeta2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_qbeta2, dbeta2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_qbeta2, dbeta2);
}

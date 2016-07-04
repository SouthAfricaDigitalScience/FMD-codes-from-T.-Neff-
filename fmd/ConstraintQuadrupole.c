/**

  \file ConstraintQuadrupole.c

  Constrain Mass Quadrupole
  Galilei invariance by considering expectation value of center of mass


  (c) 2003 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "./Gaussian.h"
#include "./gradGaussian.h"
#include "./SlaterDet.h"
#include "./gradSlaterDet.h"
#include "./CenterofMass.h"
#include "./gradCMSlaterDet.h"

#include "./Constraint.h"
#include "./ConstraintQuadrupole.h"

#include "../numerics/cmath.h"


Constraint ConstraintQ2 = {
  val : 0.0,
  label : "sQ2",
  me : calcConstraintQ2,
  gradme : calcgradConstraintQ2,
  output : outputConstraintQ2
};

Constraintod ConstraintQ2od = {
  val : 0.0,
  label : "sQ2",
  dim : 9,
  gradneedsme : 1,
  me : calcConstraintQ2od,
  gradme : calcgradConstraintQ2od,
  output : outputConstraintQ2
};

Constraint ConstraintEQ2 = {
  val : 0.0,
  label : "sEQ2",
  me : calcConstraintEQ2,
  gradme : calcgradConstraintEQ2,
  output : outputConstraintEQ2
};

Constraintod ConstraintEQ2od = {
  val : 0.0,
  label : "sEQ2",
  dim : 9,
  gradneedsme : 1,
  me : calcConstraintEQ2od,
  gradme : calcgradConstraintEQ2od,
  output : outputConstraintEQ2
};

Constraint ConstraintNQ2 = {
  val : 0.0,
  label : "sNQ2",
  me : calcConstraintNQ2,
  gradme : calcgradConstraintNQ2,
  output : outputConstraintNQ2
};

Constraintod ConstraintNQ2od = {
  val : 0.0,
  label : "sNQ2",
  dim : 9,
  gradneedsme : 1,
  me : calcConstraintNQ2od,
  gradme : calcgradConstraintNQ2od,
  output : outputConstraintNQ2
};

Constraint ConstraintDetQ = {
  val : 0.0,
  label : "cdetQ",
  me : calcConstraintDetQ,
  gradme : calcgradConstraintDetQ,
  output : outputConstraintDetQ
};

Constraint ConstraintQ2offdiagonal = {
  val : 0.0,
  label : "sQ2off",
  me : calcConstraintQ2offdiagonal,
  gradme : calcgradConstraintQ2offdiagonal,
  output : outputConstraintQ2offdiagonal
};


typedef struct {
  double Xcm[3];
  double q[9];
} quadrupolepara;


static void ob_quadrupole(const double* Xcm,
			  const Gaussian*G1, const Gaussian* G2, 
			  const GaussianAux* X, 
			  complex double q[9])
{
  int i,j;
  complex double rho2cm = 0.0;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++)
      q[i+j*3] += 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*X->Q;

    q[j+j*3] += - rho2cm* X->Q; 
  }
}


static void ob_equadrupole(const double* Xcm,
			   const Gaussian* G1, const Gaussian* G2, 
			   const GaussianAux* X, 
			   complex double q[9])
{
  int i,j;
  complex double rho2cm = 0.0;

  if (G1->xi == 1) {

    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    for (j=0; j<3; j++) {
      for (i=0; i<3; i++)
	q[i+j*3] += 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*X->Q;

      q[j+j*3] += - rho2cm* X->Q; 
    }
  }	
}


static void ob_nquadrupole(const double* Xcm,
			   const Gaussian* G1, const Gaussian* G2, 
			   const GaussianAux* X, 
			   complex double q[9])
{
  int i,j;
  complex double rho2cm = 0.0;

  if (G1->xi == -1) {

    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    for (j=0; j<3; j++) {
      for (i=0; i<3; i++)
	q[i+j*3] += 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*X->Q;

      q[j+j*3] += - rho2cm* X->Q; 
    }
  }	
}


static void ob_quadrupoleoffd(const double* Xcm,
			      const Gaussian*G1, const Gaussian* G2, 
			      const GaussianAux* X, 
			      complex double q[9])
{
  int i,j;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++)
      if (i!=j)
	q[i+j*3] += 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j])*X->Q;
  }
}


static void gob_quadrupole(quadrupolepara* q2par,
			   const Gaussian* G1, const Gaussian* G2,
			   const GaussianAux* X, const gradGaussianAux* dX,
			   complex double* q2, gradGaussian* dq2)
{
  double* Xcm = q2par->Xcm;
  double* q=q2par->q;
  int i,j, k;
  complex double merr;
  complex double rho2cm = 0.0;
  gradGaussian drho2cm;

  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);

  drho2cm.a = 0.0;
  for (i=0; i<3; i++)
    drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
  for (i=0; i<3; i++)
    drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {

      merr = 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j]);
      
      *q2 += q[i+j*3]* merr *X->Q;

      for (k=0; k<2; k++)
	dq2->chi[k] +=  q[i+j*3]* merr* dX->dQ.chi[k];
      dq2->a += q[i+j*3]*
	(3*(dX->drho.a[i]*(X->rho[j]-Xcm[j])+(X->rho[i]-Xcm[i])*dX->drho.a[j])*X->Q+
	 merr* dX->dQ.a);
      dq2->b[i] += q[i+j*3]*
	3*dX->drho.b*(X->rho[j]-Xcm[j])*X->Q;
      dq2->b[j] += q[i+j*3]*
	3*(X->rho[i]-Xcm[i])*dX->drho.b*X->Q;
      for (k=0; k<3; k++) 
	dq2->b[k] += q[i+j*3]* merr* dX->dQ.b[k];
    }

    merr = rho2cm;

    *q2 += - q[j+j*3]* merr* X->Q;

    for (k=0; k<2; k++)
      dq2->chi[k] += - q[j+j*3]* merr* dX->dQ.chi[k];
    dq2->a += - q[j+j*3]* (drho2cm.a* X->Q + merr* dX->dQ.a);
    for (k=0; k<3; k++)
      dq2->b[k] += - q[j+j*3]* 
	(drho2cm.b[k]*X->Q + merr* dX->dQ.b[k]);
  }
}


static void gob_equadrupole(quadrupolepara* q2par,
			    const Gaussian* G1, const Gaussian* G2,
			    const GaussianAux* X, const gradGaussianAux* dX,
			    complex double* q2, gradGaussian* dq2)
{	
  double* Xcm = q2par->Xcm;
  double* q=q2par->q;
  int i,j, k;
  complex double merr;
  complex double rho2cm = 0.0;
  gradGaussian drho2cm;

  if (G1->xi == 1) {

    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    drho2cm.a = 0.0;
    for (i=0; i<3; i++)
      drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
    for (i=0; i<3; i++)
      drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

    for (j=0; j<3; j++) {
      for (i=0; i<3; i++) {

	merr = 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j]);
      
	*q2 += q[i+j*3]* merr *X->Q;

	for (k=0; k<2; k++)
	  dq2->chi[k] +=  q[i+j*3]* merr* dX->dQ.chi[k];
	dq2->a += q[i+j*3]*
	  (3*(dX->drho.a[i]*(X->rho[j]-Xcm[j])+(X->rho[i]-Xcm[i])*dX->drho.a[j])*X->Q+
	   merr* dX->dQ.a);
	dq2->b[i] += q[i+j*3]*
	  3*dX->drho.b*(X->rho[j]-Xcm[j])*X->Q;
	dq2->b[j] += q[i+j*3]*
	  3*(X->rho[i]-Xcm[i])*dX->drho.b*X->Q;
	for (k=0; k<3; k++) 
	  dq2->b[k] += q[i+j*3]* merr* dX->dQ.b[k];
      }

      merr = rho2cm;

      *q2 += - q[j+j*3]* merr* X->Q;

      for (k=0; k<2; k++)
	dq2->chi[k] += - q[j+j*3]* merr* dX->dQ.chi[k];
      dq2->a += - q[j+j*3]* (drho2cm.a* X->Q + merr* dX->dQ.a);
      for (k=0; k<3; k++)
	dq2->b[k] += - q[j+j*3]* 
	  (drho2cm.b[k]*X->Q + merr* dX->dQ.b[k]);
    }
  }
}	


static void gob_nquadrupole(quadrupolepara* q2par,
			    const Gaussian* G1, const Gaussian* G2,
			    const GaussianAux* X, const gradGaussianAux* dX,
			    complex double* q2, gradGaussian* dq2)
{	
  double* Xcm = q2par->Xcm;
  double* q=q2par->q;
  int i,j, k;
  complex double merr;
  complex double rho2cm = 0.0;
  gradGaussian drho2cm;

  if (G1->xi == -1) {

    for (i=0; i<3; i++)
      rho2cm += csqr(X->rho[i]-Xcm[i]);

    drho2cm.a = 0.0;
    for (i=0; i<3; i++)
      drho2cm.a += 2*(X->rho[i]-Xcm[i])*dX->drho.a[i];
    for (i=0; i<3; i++)
      drho2cm.b[i] = 2*(X->rho[i]-Xcm[i])*dX->drho.b;

    for (j=0; j<3; j++) {
      for (i=0; i<3; i++) {

	merr = 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j]);
      
	*q2 += q[i+j*3]* merr *X->Q;

	for (k=0; k<2; k++)
	  dq2->chi[k] +=  q[i+j*3]* merr* dX->dQ.chi[k];
	dq2->a += q[i+j*3]*
	  (3*(dX->drho.a[i]*(X->rho[j]-Xcm[j])+(X->rho[i]-Xcm[i])*dX->drho.a[j])*X->Q+
	   merr* dX->dQ.a);
	dq2->b[i] += q[i+j*3]*
	  3*dX->drho.b*(X->rho[j]-Xcm[j])*X->Q;
	dq2->b[j] += q[i+j*3]*
	  3*(X->rho[i]-Xcm[i])*dX->drho.b*X->Q;
	for (k=0; k<3; k++) 
	  dq2->b[k] += q[i+j*3]* merr* dX->dQ.b[k];
      }

      merr = rho2cm;

      *q2 += - q[j+j*3]* merr* X->Q;

      for (k=0; k<2; k++)
	dq2->chi[k] += - q[j+j*3]* merr* dX->dQ.chi[k];
      dq2->a += - q[j+j*3]* (drho2cm.a* X->Q + merr* dX->dQ.a);
      for (k=0; k<3; k++)
	dq2->b[k] += - q[j+j*3]* 
	  (drho2cm.b[k]*X->Q + merr* dX->dQ.b[k]);
    }
  }
}	


static void gob_quadrupoleoffd(quadrupolepara* q2par,
			       const Gaussian* G1, const Gaussian* G2,
			       const GaussianAux* X, const gradGaussianAux* dX,
			       complex double* q2, gradGaussian* dq2)
{
  double* Xcm = q2par->Xcm;
  double* q=q2par->q;
  int i,j, k;
  complex double merr;

  for (j=0; j<3; j++)
    for (i=0; i<3; i++)
      if (i != j) {

	merr = 3*(X->rho[i]-Xcm[i])*(X->rho[j]-Xcm[j]);
      
	*q2 += q[i+j*3]* merr *X->Q;

	for (k=0; k<2; k++)
	  dq2->chi[k] +=  q[i+j*3]* merr* dX->dQ.chi[k];
	dq2->a += q[i+j*3]*
	  (3*(dX->drho.a[i]*(X->rho[j]-Xcm[j])+(X->rho[i]-Xcm[i])*dX->drho.a[j])*X->Q+
	   merr* dX->dQ.a);
	dq2->b[i] += q[i+j*3]*
	  3*dX->drho.b*(X->rho[j]-Xcm[j])*X->Q;
	dq2->b[j] += q[i+j*3]*
	  3*(X->rho[i]-Xcm[i])*dX->drho.b*X->Q;
	for (k=0; k<3; k++) 
	  dq2->b[k] += q[i+j*3]* merr* dX->dQ.b[k];
      }
}


static void gcmob_equadrupole(const quadrupolepara* P,
			      const Gaussian* G1, const Gaussian* G2,
			      const GaussianAux* X,
			      gradCM* dCMq2)
{
  double* Xcm = P->Xcm;
  double* q = P->q;
  int i,j;

  if (G1->xi == 1) {

    for (j=0; j<3; j++) {
      for (i=0; i<3; i++) {
	dCMq2->X[i] -= q[i+j*3]* 3*(X->rho[j]-Xcm[j])* X->Q;
       	dCMq2->X[j] -= q[i+j*3]* 3*(X->rho[i]-Xcm[i])* X->Q;
      }
      dCMq2->X[j] -= q[j+j*3]* 2*(X->rho[j]-Xcm[j])* X->Q;
    }
  }
}


static void gcmob_nquadrupole(const quadrupolepara* P,
			      const Gaussian* G1, const Gaussian* G2,
			      const GaussianAux* X,
			      gradCM* dCMq2)
{
  double* Xcm = P->Xcm;
  double* q = P->q;
  int i,j;

  if (G1->xi == -1) {

    for (j=0; j<3; j++) {
      for (i=0; i<3; i++) {
	dCMq2->X[i] -= q[i+j*3]* 3*(X->rho[j]-Xcm[j])* X->Q;
       	dCMq2->X[j] -= q[i+j*3]* 3*(X->rho[i]-Xcm[i])* X->Q;
      }
      dCMq2->X[j] -= q[j+j*3]* 2*(X->rho[j]-Xcm[j])* X->Q;
    }
  }
}


#define SQR(x) (x)*(x)

void calcConstraintQ2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* q2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);

  *q2 = sqrt(qsq);
}


void calcConstraintQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X, 
			complex double *q)
{
  double Xcm[3] = { 0.0, 0.0, 0.0 };

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_q, q);
}


void calcConstraintEQ2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* q2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_equadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);

  *q2 = sqrt(qsq);
}


void calcConstraintEQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X, 
			 complex double *q)
{
  double Xcm[3] = { 0.0, 0.0, 0.0 };

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_equadrupole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_q, q);
}


void calcConstraintNQ2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* q2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_nquadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);

  *q2 = sqrt(qsq);
}


void calcConstraintNQ2od(const SlaterDet* Q, const SlaterDet* Qp, 
			 const SlaterDetAux* X, 
			 complex double *q)
{
  double Xcm[3] = { 0.0, 0.0, 0.0 };

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_nquadrupole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_q, q);
}


void calcgradConstraintQ2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dq2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;

  int i;
  quadrupolepara q2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);
  
  for (i=0; i<3; i++)
    q2par.Xcm[i] = Xcm[i];
  for (i=0; i<9; i++)
    q2par.q[i] = q[i]/sqrt(qsq);

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_quadrupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_q2, dq2);
}


void calcgradConstraintQ2od(const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    const double* q,
			    gradSlaterDet* dq2)
{
  double Xcm[3] = { 0.0, 0.0, 0.0 };

  int i;
  quadrupolepara q2par;

  for (i=0; i<9; i++)
    q2par.q[i] = q[i];

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_quadrupole};

  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_q2, dq2);
}


void calcgradConstraintEQ2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dq2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;

  int i;
  quadrupolepara q2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_equadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);
  
  for (i=0; i<3; i++)
    q2par.Xcm[i] = Xcm[i];
  for (i=0; i<9; i++)
    q2par.q[i] = q[i]/sqrt(qsq);

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_equadrupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_q2, dq2);

  OneBodyOperator gcmop_ob_q2 = {dim: 6, opt: 1, par: &q2par, me: gcmob_equadrupole};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_q2, dq2);
  dq2->val = sqrt(qsq);
}


void calcgradConstraintEQ2od(const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* q,
			     gradSlaterDet* dq2)
{
  double Xcm[3] = { 0.0, 0.0, 0.0 };

  int i;
  quadrupolepara q2par;

  for (i=0; i<9; i++)
    q2par.q[i] = q[i];

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_equadrupole};

  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_q2, dq2);
}


void calcgradConstraintNQ2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dq2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;

  int i;
  quadrupolepara q2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_nquadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);
  
  for (i=0; i<3; i++)
    q2par.Xcm[i] = Xcm[i];
  for (i=0; i<9; i++)
    q2par.q[i] = q[i]/sqrt(qsq);

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_nquadrupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_q2, dq2);

  OneBodyOperator gcmop_ob_q2 = {dim: 6, opt: 1, par: &q2par, me: gcmob_nquadrupole};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_q2, dq2);
  dq2->val = sqrt(qsq);
}


void calcgradConstraintNQ2od(const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* q,
			     gradSlaterDet* dq2)
{
  double Xcm[3] = { 0.0, 0.0, 0.0 };

  int i;
  quadrupolepara q2par;

  for (i=0; i<9; i++)
    q2par.q[i] = q[i];

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_nquadrupole};

  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_q2, dq2);
}


double outputConstraintQ2(double val)
{
  return val;
}


double outputConstraintEQ2(double val)
{
  return val;
}


double outputConstraintNQ2(double val)
{
  return val;
}


static double det(double R[3][3])
{
  return (R[0][0]*R[1][1]*R[2][2] + R[0][1]*R[1][2]*R[2][0] + 
          R[0][2]*R[1][0]*R[2][1] - R[0][0]*R[1][2]*R[2][1] - 
          R[0][1]*R[1][0]*R[2][2] - R[0][2]*R[1][1]*R[2][0]);
}


void calcConstraintDetQ(const SlaterDet* Q, const SlaterDetAux* X, 
			double* detq)
{
  double Xcm[3];
  double q[9];

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  *detq = cbrt(det(q));
}

static void adjunct(const double q[3][3], double a[3][3])
{
  int i,j;

  for (j=0; j<3; j++){
    for (i=0; i<3; i++){
      a[i][j] = 
	q[(j+1)%3][(i+1)%3]*q[(j+2)%3][(i+2)%3]-
	q[(j+1)%3][(i+2)%3]*q[(j+2)%3][(i+1)%3];
      printf("%f\n", a[i][j]);
                       }
                      }
}


void calcgradConstraintDetQ(const SlaterDet* Q, const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* ddetq)
{
  double Xcm[3];
  double q[9], detq, adjq[9];

  int i;
  quadrupolepara q2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<3; i++)
    q2par.Xcm[i] = Xcm[i];

  detq = det(q);
  adjunct(q, adjq);

  for (i=0; i<9; i++)
    q2par.q[i] = 1.0/(3*SQR(cbrt(detq)))*adjq[i];

  gradOneBodyOperator gop_ob_detq = {opt: 1, par: &q2par, me: gob_quadrupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_detq, ddetq);
  ddetq->val = cbrt(detq);
}


double outputConstraintDetQ(double val)
{
  return 3*val;
}



void calcConstraintQ2offdiagonal(const SlaterDet* Q, const SlaterDetAux* X, 
				 double* q2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;
  int i;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupoleoffd};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);

  *q2 = sqrt(qsq);
}


void calcgradConstraintQ2offdiagonal(const SlaterDet* Q, const SlaterDetAux* X,
				     const gradSlaterDetAux* dX,
				     gradSlaterDet* dq2)
{
  double Xcm[3];
  double q[9], qsq = 0.0;

  int i;
  quadrupolepara q2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_quadrupoleoffd};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<9; i++)
    qsq += SQR(q[i]);
  
  for (i=0; i<3; i++)
    q2par.Xcm[i] = Xcm[i];
  for (i=0; i<9; i++)
    q2par.q[i] = q[i]/sqrt(qsq);

  gradOneBodyOperator gop_ob_q2 = {opt: 1, par: &q2par, me: gob_quadrupoleoffd};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_q2, dq2);
}


double outputConstraintQ2offdiagonal(double val)
{
  return val;
}

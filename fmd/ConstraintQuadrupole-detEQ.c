/**

  \file ConstraintQuadrupole-detEQ.c

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
#include "./ConstraintQuadrupole-detEQ.h"

#include "../numerics/cmath.h"




Constraint ConstraintDetEQ = {
  val : 0.0,
  label : "cdetEQ",
  me : calcConstraintDetEQ,
  gradme : calcgradConstraintDetEQ,
  output : outputConstraintDetEQ
};


typedef struct {
  double Xcm[3];
  double q[9];
} quadrupolepara;




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


#define SQR(x) (x)*(x)


static double det(double R[3][3])
{
  return (R[0][0]*R[1][1]*R[2][2] + R[0][1]*R[1][2]*R[2][0] + 
          R[0][2]*R[1][0]*R[2][1] - R[0][0]*R[1][2]*R[2][1] - 
          R[0][1]*R[1][0]*R[2][2] - R[0][2]*R[1][1]*R[2][0]);
}


void calcConstraintDetEQ(const SlaterDet* Q, const SlaterDetAux* X, 
			double* deteq)
{
  double Xcm[3];
  double q[9];

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_equadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  *deteq = cbrt(det(q));
}

/** Print out values for quadrupole array */

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


void calcgradConstraintDetEQ(const SlaterDet* Q, const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* ddeteq)
{
  double Xcm[3];
  double q[9], deteq, adjq[9];

  int i;
  quadrupolepara q2par;

  calcCMPosition(Q, X, Xcm);

  OneBodyOperator op_ob_q = {dim: 9, opt: 1, par: Xcm, me: ob_equadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_q, q);

  for (i=0; i<3; i++)
    q2par.Xcm[i] = Xcm[i];

  deteq = det(q);
  adjunct(q, adjq);

  for (i=0; i<9; i++)
    q2par.q[i] = 1.0/(3*SQR(cbrt(deteq)))*adjq[i];

  gradOneBodyOperator gop_ob_detq = {opt: 1, par: &q2par, me: gob_equadrupole};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_detq, ddeteq);
  ddeteq->val = cbrt(deteq);
}

double outputConstraintDetEQ(double val)
{
  return 3*val;
}




/**

  \file ConstraintR2.c

  Constrain Mass/Charge Radius


  (c) 2003 Thomas Neff

*/

#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Constraint.h"
#include "ConstraintR2.h"

#include "numerics/cmath.h"


Constraint ConstraintR2 = {
  val : 0.0,
  label : "R",
  me : calcConstraintR2,
  gradme : calcgradConstraintR2,
  output : outputConstraintR2
};

Constraintod ConstraintR2od = {
  val : 0.0,
  label : "R",
  dim : 1,
  gradneedsme : 0,
  me : calcConstraintR2od,
  gradme : calcgradConstraintR2od,
  output : outputConstraintR2
};

Constraint ConstraintER2 = {
  val : 0.0,
  label : "ER",
  me : calcConstraintER2,
  gradme : calcgradConstraintER2,
  output : outputConstraintER2
};

Constraintod ConstraintER2od = {
  val : 0.0,
  label : "ER",
  dim : 1,
  gradneedsme : 0,
  me : calcConstraintER2od,
  gradme : calcgradConstraintER2od,
  output : outputConstraintR2
};

Constraint ConstraintNR2 = {
  val : 0.0,
  label : "NR",
  me : calcConstraintNR2,
  gradme : calcgradConstraintNR2,
  output : outputConstraintNR2
};

Constraintod ConstraintNR2od = {
  val : 0.0,
  label : "NR",
  dim : 1,
  gradneedsme : 0,
  me : calcConstraintNR2od,
  gradme : calcgradConstraintNR2od,
  output : outputConstraintR2
};

// matrix elements as in Radii.c

static void ob_r2(const SlaterDet* Q,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* r2)
{
  int A=Q->A;

  *r2 += 1.0/A*(1.0-1.0/A)* (3*X->alpha + X->rho2)* X->Q;
}

static void ob_er2(const SlaterDet* Q,
		   const Gaussian*G1, const Gaussian* G2, 
		   const GaussianAux* X, complex double* r2)
{
  int A=Q->A;
  int Z=Q->Z;

  // radius with respec to mass cm
  // *r2 += (1.0/(A*A) + (1+G1->xi)/2 *(1.0/Z-2.0/(A*Z)))*         
  //   (3*X->alpha + X->rho2)* X->Q;

  // radius with respect to proton cm
  *r2 += 1.0/Z*(1.0-1.0/Z)* (1+G1->xi)/2* (3*X->alpha + X->rho2)* X->Q;
}

static void ob_nr2(const SlaterDet* Q,
		   const Gaussian*G1, const Gaussian* G2, 
		   const GaussianAux* X, complex double* r2)
{
  int A=Q->A;
  int N=Q->N;

  // radius with respec to mass cm
  // *r2 += (1.0/(A*A) + (1-G1->xi)/2 *(1.0/N-2.0/(A*N)))* 
  //   (3*X->alpha + X->rho2)* X->Q;

  // radius with respect to neutron cm
  *r2 += 1.0/N*(1.0-1.0/N)* (1-G1->xi)/2* (3*X->alpha + X->rho2)* X->Q;
}


static void tb_r2(const SlaterDet* Q,
		  const Gaussian* G1, const Gaussian* G2, 
		  const Gaussian* G3, const Gaussian* G4, 
		  const GaussianAux* X13, const GaussianAux* X24, 
		  complex double* r2)
{	
  int A=Q->A;

  *r2 += -2.0/(A*A)* cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;
}

static void tb_er2(const SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24, 
		   complex double* r2)
{	
  int A=Q->A;
  int Z=Q->Z;

  // radius with respect to mass cm
  // *r2 += 2.0*(1.0/(A*A) - ((1+G1->xi)/2+(1+G2->xi)/2) *1.0/(A*Z))*
  //   cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;

  // radius with respect to proton cm
  *r2 += -2.0/(Z*Z)* (1+G1->xi)/2* (1+G2->xi)/2*
    cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;
}

static void tb_nr2(const SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24, 
		   complex double* r2)
{	
  int A=Q->A;
  int N=Q->N;

  // radius with respect to mass cm
  // *r2 += 2.0*(1.0/(A*A) - ((1-G1->xi)/2+(1-G2->xi)/2) *1.0/(A*N))*
  //   cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;

  // radius with respect to neutron cm
  *r2 += -2.0/(N*N)* (1-G1->xi)/2* (1-G2->xi)/2*
    cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;
}


static void gob_r2(const SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2,
		   const GaussianAux* X, const gradGaussianAux* dX,
		   complex double* r2, gradGaussian* dr2)
{
  int A=Q->A;
  complex double mr2;
  double factor;
  int i;

  mr2 = 3*X->alpha + X->rho2;
  factor = 1.0/A*(1.0-1.0/A);

  *r2 += factor* mr2* X->Q;

  for (i=0; i<2; i++)
    dr2->chi[i] += factor* mr2* dX->dQ.chi[i];
  dr2->a += factor*
    ((3*dX->dalpha + dX->drho2.a)* X->Q + mr2* dX->dQ.a);
  for (i=0; i<3; i++) 
    dr2->b[i] += factor*
      (dX->drho2.b[i]* X->Q + mr2* dX->dQ.b[i]);
}

static void gob_er2(const SlaterDet* Q,
		    const Gaussian* G1, const Gaussian* G2,
		    const GaussianAux* X, const gradGaussianAux* dX,
		    complex double* r2, gradGaussian* dr2)
{
  int A=Q->A;
  int Z=Q->Z;
  complex double mr2;
  double factor;
  int i;

  mr2 = 3*X->alpha + X->rho2;
  // factor = (1.0/(A*A) + (1+G1->xi)/2 *(1.0/Z-2.0/(A*Z)));
  factor = 1.0/Z*(1.0-1.0/Z)* (1+G1->xi)/2;

  *r2 += factor* mr2* X->Q;

  for (i=0; i<2; i++)
    dr2->chi[i] += factor* mr2* dX->dQ.chi[i];
  dr2->a += factor*
    ((3*dX->dalpha + dX->drho2.a)* X->Q + mr2* dX->dQ.a);
  for (i=0; i<3; i++) 
    dr2->b[i] += factor*
      (dX->drho2.b[i]* X->Q + mr2* dX->dQ.b[i]);
}

static void gob_nr2(const SlaterDet* Q,
		    const Gaussian* G1, const Gaussian* G2,
		    const GaussianAux* X, const gradGaussianAux* dX,
		    complex double* r2, gradGaussian* dr2)
{
  int A=Q->A;
  int N=Q->N;
  complex double mr2;
  double factor;
  int i;

  mr2 = 3*X->alpha + X->rho2;
  // factor = (1.0/(A*A) + (1-G1->xi)/2 *(1.0/N-2.0/(A*N)));
  factor = 1.0/N*(1.0-1.0/N)* (1-G1->xi)/2;

  *r2 += factor* mr2* X->Q;

  for (i=0; i<2; i++)
    dr2->chi[i] += factor* mr2* dX->dQ.chi[i];
  dr2->a += factor*
    ((3*dX->dalpha + dX->drho2.a)* X->Q + mr2* dX->dQ.a);
  for (i=0; i<3; i++) 
    dr2->b[i] += factor*
      (dX->drho2.b[i]* X->Q + mr2* dX->dQ.b[i]);
}


static void gtb_r2(const SlaterDet* Q,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24,
		   const gradGaussianAux* dX13,
		   complex double* r2, gradGaussian* dr2)
{
  int A=Q->A;
  int i;
  complex double mr2;
  double factor;

  mr2 = cvec3mult(X13->rho, X24->rho);
  factor = -2.0/(A*A);

  *r2 += factor* mr2* X13->Q*X24->Q;

  for (i=0; i<2; i++)
    dr2->chi[i] += factor* mr2* dX13->dQ.chi[i]* X24->Q; 

  dr2->a += factor*
    (cvec3mult(dX13->drho.a, X24->rho)* X13->Q*X24->Q +
     mr2* dX13->dQ.a* X24->Q);
  for (i=0; i<3; i++)
    dr2->b[i] += factor*
      (dX13->drho.b* X24->rho[i]* X13->Q* X24->Q +
       mr2* dX13->dQ.b[i]* X24->Q);
}

static void gtb_er2(const SlaterDet* Q,
		    const Gaussian* G1, const Gaussian* G2, 
		    const Gaussian* G3, const Gaussian* G4, 
		    const GaussianAux* X13, const GaussianAux* X24,
		    const gradGaussianAux* dX13,
		    complex double* r2, gradGaussian* dr2)
{
  int A=Q->A;
  int Z=Q->Z;
  int i;
  complex double mr2;
  double factor;

  mr2 = cvec3mult(X13->rho, X24->rho);
  // factor = 2.0*(1.0/(A*A) - ((1+G1->xi)/2+(1+G2->xi)/2) *1.0/(A*Z));
  factor = -2.0/(Z*Z)* (1+G1->xi)/2* (1+G2->xi)/2;

  *r2 += factor* mr2* X13->Q*X24->Q;

  for (i=0; i<2; i++)
    dr2->chi[i] += factor* mr2* dX13->dQ.chi[i]* X24->Q; 

  dr2->a += factor*
    (cvec3mult(dX13->drho.a, X24->rho)* X13->Q*X24->Q +
     mr2* dX13->dQ.a* X24->Q);
  for (i=0; i<3; i++)
    dr2->b[i] += factor*
      (dX13->drho.b* X24->rho[i]* X13->Q* X24->Q +
       mr2* dX13->dQ.b[i]* X24->Q);
}

static void gtb_nr2(const SlaterDet* Q,
		    const Gaussian* G1, const Gaussian* G2, 
		    const Gaussian* G3, const Gaussian* G4, 
		    const GaussianAux* X13, const GaussianAux* X24,
		    const gradGaussianAux* dX13,
		    complex double* r2, gradGaussian* dr2)
{
  int A=Q->A;
  int N=Q->N;
  int i;
  complex double mr2;
  double factor;

  mr2 = cvec3mult(X13->rho, X24->rho);
  // factor = 2.0*(1.0/(A*A) - ((1-G1->xi)/2+(1-G2->xi)/2) *1.0/(A*N));
  factor = -2.0/(N*N)* (1+G1->xi)/2* (1+G2->xi)/2;

  *r2 += factor* mr2* X13->Q*X24->Q;

  for (i=0; i<2; i++)
    dr2->chi[i] += factor* mr2* dX13->dQ.chi[i]* X24->Q; 

  dr2->a += factor*
    (cvec3mult(dX13->drho.a, X24->rho)* X13->Q*X24->Q +
     mr2* dX13->dQ.a* X24->Q);
  for (i=0; i<3; i++)
    dr2->b[i] += factor*
      (dX13->drho.b* X24->rho[i]* X13->Q* X24->Q +
       mr2* dX13->dQ.b[i]* X24->Q);
}


void calcConstraintR2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* r2)
{
  double r2one, r2two;
  OneBodyOperator op_ob_r2 = {dim: 1, opt: 1, par: Q, me: ob_r2};
  TwoBodyOperator op_tb_r2 = {dim: 1, opt: 1, par: Q, me: tb_r2};


  calcSlaterDetOBME(Q, X, &op_ob_r2, &r2one);
  calcSlaterDetTBME(Q, X, &op_tb_r2, &r2two);

  *r2 = r2one + r2two;
}


void calcgradConstraintR2(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dr2)
{
  gradOneBodyOperator gop_ob_r2 = {opt: 1, par: Q, me: gob_r2};
  gradTwoBodyOperator gop_tb_r2 = {opt: 1, par: Q, me: gtb_r2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_r2, dr2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_r2, dr2);
}


void calcConstraintR2od(const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X, 
			complex double* r2)
{
  complex double r2one, r2two;
  OneBodyOperator op_ob_r2 = {dim: 1, opt: 1, par: Q, me: ob_r2};
  TwoBodyOperator op_tb_r2 = {dim: 1, opt: 1, par: Q, me: tb_r2};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_r2, &r2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_r2, &r2two);

  *r2 = r2one + r2two;
}


void calcgradConstraintR2od(const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X,
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dr2)
{
  gradOneBodyOperator gop_ob_r2 = {opt: 1, par: Q, me: gob_r2};
  gradTwoBodyOperator gop_tb_r2 = {opt: 1, par: Q, me: gtb_r2};


  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_r2, dr2);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_r2, dr2);
}


double outputConstraintR2(double val)
{
  return(sqrt(val));
}


void calcConstraintER2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* r2)
{
  double r2one, r2two;
  OneBodyOperator op_ob_r2 = {dim: 1, opt: 1, par: Q, me: ob_er2};
  TwoBodyOperator op_tb_r2 = {dim: 1, opt: 1, par: Q, me: tb_er2};


  calcSlaterDetOBME(Q, X, &op_ob_r2, &r2one);
  calcSlaterDetTBME(Q, X, &op_tb_r2, &r2two);

  *r2 = r2one + r2two;
}


void calcgradConstraintER2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dr2)
{
  gradOneBodyOperator gop_ob_r2 = {opt: 1, par: Q, me: gob_er2};
  gradTwoBodyOperator gop_tb_r2 = {opt: 1, par: Q, me: gtb_er2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_r2, dr2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_r2, dr2);
}


void calcConstraintER2od(const SlaterDet* Q, const SlaterDet* Qp,
			 const SlaterDetAux* X, 
			 complex double* r2)
{
  complex double r2one, r2two;
  OneBodyOperator op_ob_r2 = {dim: 1, opt: 1, par: Q, me: ob_er2};
  TwoBodyOperator op_tb_r2 = {dim: 1, opt: 1, par: Q, me: tb_er2};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_r2, &r2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_r2, &r2two);

  *r2 = r2one + r2two;
}


void calcgradConstraintER2od(const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* dummy,
			     gradSlaterDet* dr2)
{
  gradOneBodyOperator gop_ob_r2 = {opt: 1, par: Q, me: gob_er2};
  gradTwoBodyOperator gop_tb_r2 = {opt: 1, par: Q, me: gtb_er2};


  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_r2, dr2);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_r2, dr2);
}


double outputConstraintER2(double val)
{
  if (val > 0.0)
    return(sqrt(val));
  else
    return(-sqrt(-val));
}


void calcConstraintNR2(const SlaterDet* Q, const SlaterDetAux* X, 
		       double* r2)
{
  double r2one, r2two;
  OneBodyOperator op_ob_r2 = {dim: 1, opt: 1, par: Q, me: ob_nr2};
  TwoBodyOperator op_tb_r2 = {dim: 1, opt: 1, par: Q, me: tb_nr2};


  calcSlaterDetOBME(Q, X, &op_ob_r2, &r2one);
  calcSlaterDetTBME(Q, X, &op_tb_r2, &r2two);

  *r2 = r2one + r2two;
}


void calcgradConstraintNR2(const SlaterDet* Q, const SlaterDetAux* X,
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dr2)
{
  gradOneBodyOperator gop_ob_r2 = {opt: 1, par: Q, me: gob_nr2};
  gradTwoBodyOperator gop_tb_r2 = {opt: 1, par: Q, me: gtb_nr2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_r2, dr2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_r2, dr2);
}


void calcConstraintNR2od(const SlaterDet* Q, const SlaterDet* Qp,
			 const SlaterDetAux* X, 
			 complex double* r2)
{
  complex double r2one, r2two;
  OneBodyOperator op_ob_r2 = {dim: 1, opt: 1, par: Q, me: ob_nr2};
  TwoBodyOperator op_tb_r2 = {dim: 1, opt: 1, par: Q, me: tb_nr2};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_r2, &r2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_r2, &r2two);

  *r2 = r2one + r2two;
}


void calcgradConstraintNR2od(const SlaterDet* Q, const SlaterDet* Qp,
			     const SlaterDetAux* X,
			     const gradSlaterDetAux* dX,
			     const double* dummy,
			     gradSlaterDet* dr2)
{
  gradOneBodyOperator gop_ob_r2 = {opt: 1, par: Q, me: gob_nr2};
  gradTwoBodyOperator gop_tb_r2 = {opt: 1, par: Q, me: gtb_nr2};


  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_r2, dr2);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_r2, dr2);
}


double outputConstraintNR2(double val)
{
  if (val > 0.0)
    return(sqrt(val));
  else
    return(-sqrt(-val));
}

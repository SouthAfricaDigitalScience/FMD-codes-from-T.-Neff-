/**

  \file ConstraintT2.c

  Constrain expectation value of T2


  (c) 2004 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Constraint.h"
#include "ConstraintT2.h"

#include "numerics/cmath.h"


Constraint ConstraintT2 = {
  val : 0.0,
  label : "T2",
  me : calcConstraintT2,
  gradme : calcgradConstraintT2,
  output : outputConstraintT2
};

Constraint ConstraintPT2 = {
  val : 0.0,
  label : "PT2",
  me : calcConstraintPT2,
  gradme : calcgradConstraintPT2,
  output : outputConstraintT2
};

Constraint ConstraintNT2 = {
  val : 0.0,
  label : "NT2",
  me : calcConstraintNT2,
  gradme : calcgradConstraintNT2,
  output : outputConstraintT2
};


Constraintod ConstraintT2od = {
  val : 0.0,
  label : "T2",
  dim : 1,
  gradneedsme : 0,
  me : calcConstraintT2od,
  gradme : calcgradConstraintT2od,
  output : outputConstraintT2
};


static void ob_t2(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* t2)
{
  *t2 += 0.75*X->Q;
}

static void ob_pt2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const GaussianAux* X, complex double* t2)
{
  if (G1->xi == +1)
    *t2 += 0.75*X->Q;
}

static void ob_nt2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const GaussianAux* X, complex double* t2)
{
  if (G1->xi == -1)
    *t2 += 0.75*X->Q;
}

static void tb_t2(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const Gaussian* G3, const Gaussian* G4, 
		  const GaussianAux* X13, const GaussianAux* X24, 
		  complex double* t2)
{	
  *t2 += 0.5*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
	   (G1->xi*G4->xi+G2->xi*G3->xi)/2)* 
    X13->S* X24->S* X13->R* X24->R;
}

static void tb_pt2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const Gaussian* G3, const Gaussian* G4, 
                   const GaussianAux* X13, const GaussianAux* X24, 
                   complex double* t2)
{	
  if (G1->xi == +1 && G2->xi == +1)
    *t2 += 0.5*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
                (G1->xi*G4->xi+G2->xi*G3->xi)/2)* 
      X13->S* X24->S* X13->R* X24->R;
}

static void tb_nt2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const Gaussian* G3, const Gaussian* G4, 
                   const GaussianAux* X13, const GaussianAux* X24, 
                   complex double* t2)
{	
  if (G1->xi == -1 && G2->xi == -1)
    *t2 += 0.5*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
                (G1->xi*G4->xi+G2->xi*G3->xi)/2)* 
      X13->S* X24->S* X13->R* X24->R;
}


static void gob_t2(void* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX,
		   complex double* t2, gradGaussian* dt2)
{
  int i;
  complex double t2m;

  t2m = 0.75*X->T;

  *t2 += t2m*X->S*X->R;

  for (i=0; i<2; i++)
    dt2->chi[i] += t2m*dX->dS.chi[i]*X->R;

  dt2->a += t2m* X->S* dX->dR.a;

  for (i=0; i<3; i++)
    dt2->b[i] += t2m* X->S* dX->dR.b[i];
}

static void gob_pt2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const GaussianAux* X, const gradGaussianAux* dX,
                    complex double* t2, gradGaussian* dt2)
{
  int i;
  complex double t2m;

  if (G1->xi == +1) {
    t2m = 0.75*X->T;

    *t2 += t2m*X->S*X->R;

    for (i=0; i<2; i++)
      dt2->chi[i] += t2m*dX->dS.chi[i]*X->R;

    dt2->a += t2m* X->S* dX->dR.a;

    for (i=0; i<3; i++)
      dt2->b[i] += t2m* X->S* dX->dR.b[i];
  }
}

static void gob_nt2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const GaussianAux* X, const gradGaussianAux* dX,
                    complex double* t2, gradGaussian* dt2)
{
  int i;
  complex double t2m;

  if (G1->xi == -1) {
    t2m = 0.75*X->T;

    *t2 += t2m*X->S*X->R;

    for (i=0; i<2; i++)
      dt2->chi[i] += t2m*dX->dS.chi[i]*X->R;

    dt2->a += t2m* X->S* dX->dR.a;

    for (i=0; i<3; i++)
      dt2->b[i] += t2m* X->S* dX->dR.b[i];
  }
}


static void gtb_t2(void* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24,
		   const gradGaussianAux* dX13,
		   complex double* t2, gradGaussian* dt2)
{	
  complex double t2m;
  int i;

  t2m = 0.25*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
	      (G1->xi*G4->xi+G2->xi*G3->xi)/2);

  *t2 += 2* t2m* X13->S*X24->S* X13->R*X24->R;
  
  for (i=0; i<2; i++)
    dt2->chi[i] += 2* t2m* dX13->dS.chi[i]*X24->S* X13->R*X24->R;

  dt2->a += 2* t2m* X13->S*X24->S* dX13->dR.a*X24->R;

  for (i=0; i<3; i++)
    dt2->b[i] += 2* t2m* X13->S*X24->S* dX13->dR.b[i]*X24->R;
}

static void gtb_pt2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const Gaussian* G3, const Gaussian* G4, 
                    const GaussianAux* X13, const GaussianAux* X24,
                    const gradGaussianAux* dX13,
                    complex double* t2, gradGaussian* dt2)
{	
  complex double t2m;
  int i;

  if (G1->xi == +1 && G2->xi == +1) {

    t2m = 0.25*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
                (G1->xi*G4->xi+G2->xi*G3->xi)/2); 
    
    *t2 += 2* t2m* X13->S*X24->S* X13->R*X24->R;
  
    for (i=0; i<2; i++)
      dt2->chi[i] += 2* t2m* dX13->dS.chi[i]*X24->S* X13->R*X24->R;

    dt2->a += 2* t2m* X13->S*X24->S* dX13->dR.a*X24->R;
  
    for (i=0; i<3; i++)
      dt2->b[i] += 2* t2m* X13->S*X24->S* dX13->dR.b[i]*X24->R;
  }
}

static void gtb_nt2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const Gaussian* G3, const Gaussian* G4, 
                    const GaussianAux* X13, const GaussianAux* X24,
                    const gradGaussianAux* dX13,
                    complex double* t2, gradGaussian* dt2)
{	
  complex double t2m;
  int i;

  if (G1->xi == -1 && G2->xi == -1) {

    t2m = 0.25*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
                (G1->xi*G4->xi+G2->xi*G3->xi)/2); 
    
    *t2 += 2* t2m* X13->S*X24->S* X13->R*X24->R;
  
    for (i=0; i<2; i++)
      dt2->chi[i] += 2* t2m* dX13->dS.chi[i]*X24->S* X13->R*X24->R;

    dt2->a += 2* t2m* X13->S*X24->S* dX13->dR.a*X24->R;
  
    for (i=0; i<3; i++)
      dt2->b[i] += 2* t2m* X13->S*X24->S* dX13->dR.b[i]*X24->R;
  }
}

void calcConstraintT2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* t2)
{
  double t2one, t2two;
  OneBodyOperator op_ob_t2 = {dim: 1, opt: 1, par: NULL, me: ob_t2};
  TwoBodyOperator op_tb_t2 = {dim: 1, opt: 0, par: NULL, me: tb_t2};


  calcSlaterDetOBME(Q, X, &op_ob_t2, &t2one);
  calcSlaterDetTBME(Q, X, &op_tb_t2, &t2two);

  *t2 = t2one + t2two;
}

void calcConstraintPT2(const SlaterDet* Q, const SlaterDetAux* X, 
                       double* t2)
{
  double t2one, t2two;
  OneBodyOperator op_ob_t2 = {dim: 1, opt: 1, par: NULL, me: ob_pt2};
  TwoBodyOperator op_tb_t2 = {dim: 1, opt: 0, par: NULL, me: tb_pt2};


  calcSlaterDetOBME(Q, X, &op_ob_t2, &t2one);
  calcSlaterDetTBME(Q, X, &op_tb_t2, &t2two);

  *t2 = t2one + t2two;
}

void calcConstraintNT2(const SlaterDet* Q, const SlaterDetAux* X, 
                       double* t2)
{
  double t2one, t2two;
  OneBodyOperator op_ob_t2 = {dim: 1, opt: 1, par: NULL, me: ob_nt2};
  TwoBodyOperator op_tb_t2 = {dim: 1, opt: 0, par: NULL, me: tb_nt2};


  calcSlaterDetOBME(Q, X, &op_ob_t2, &t2one);
  calcSlaterDetTBME(Q, X, &op_tb_t2, &t2two);

  *t2 = t2one + t2two;
}



void calcgradConstraintT2(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dt2)
{
  gradOneBodyOperator gop_ob_t2 = {opt: 1, par: NULL, me: gob_t2};
  gradTwoBodyOperator gop_tb_t2 = {opt: 0, par: NULL, me: gtb_t2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_t2, dt2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_t2, dt2);
}

void calcgradConstraintPT2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dt2)
{
  gradOneBodyOperator gop_ob_t2 = {opt: 1, par: NULL, me: gob_pt2};
  gradTwoBodyOperator gop_tb_t2 = {opt: 0, par: NULL, me: gtb_pt2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_t2, dt2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_t2, dt2);
}

void calcgradConstraintNT2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* dt2)
{
  gradOneBodyOperator gop_ob_t2 = {opt: 1, par: NULL, me: gob_nt2};
  gradTwoBodyOperator gop_tb_t2 = {opt: 0, par: NULL, me: gtb_nt2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_t2, dt2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_t2, dt2);
}


void calcConstraintT2od(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X, 
			complex double* t2)
{
  complex double t2one, t2two;
  OneBodyOperator op_ob_t2 = {dim: 1, opt: 1, par: NULL, me: ob_t2};
  TwoBodyOperator op_tb_t2 = {dim: 1, opt: 0, par: NULL, me: tb_t2};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_t2, &t2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_t2, &t2two);

  *t2 = t2one + t2two;
}



void calcgradConstraintT2od(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dt2)
{
  gradOneBodyOperator gop_ob_t2 = {opt: 1, par: NULL, me: gob_t2};
  gradTwoBodyOperator gop_tb_t2 = {opt: 0, par: NULL, me: gtb_t2};


  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_t2, dt2);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_t2, dt2);
}


double outputConstraintT2(double val)
{
  return (val);
}

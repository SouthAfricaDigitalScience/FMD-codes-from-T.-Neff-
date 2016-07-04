/**

  \file ConstraintS2.c

  Constrain expectation value of S2


  (c) 2003 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Constraint.h"
#include "ConstraintS2.h"

#include "numerics/cmath.h"


Constraint ConstraintS2 = {
  val : 0.0,
  label : "S2",
  me : calcConstraintS2,
  gradme : calcgradConstraintS2,
  output : outputConstraintS2
};

Constraint ConstraintPS2 = {
  val : 0.0,
  label : "PS2",
  me : calcConstraintPS2,
  gradme : calcgradConstraintPS2,
  output : outputConstraintS2
};

Constraint ConstraintNS2 = {
  val : 0.0,
  label : "NS2",
  me : calcConstraintNS2,
  gradme : calcgradConstraintNS2,
  output : outputConstraintS2
};


Constraintod ConstraintS2od = {
  val : 0.0,
  label : "S2",
  dim : 1,
  gradneedsme : 0,
  me : calcConstraintS2od,
  gradme : calcgradConstraintS2od,
  output : outputConstraintS2
};


// matrix elements identical as in AngularMomentum.c
// S2 is calculated relativ to the origin 
// should really be calculate relativ to center of mass 

static void ob_s2(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* s2)
{
  *s2 += 0.75*X->Q;
}

static void ob_ps2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const GaussianAux* X, complex double* s2)
{
  if (G1->xi == +1)
    *s2 += 0.75*X->Q;
}

static void ob_ns2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const GaussianAux* X, complex double* s2)
{
  if (G1->xi == -1)
    *s2 += 0.75*X->Q;
}


static void tb_s2(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const Gaussian* G3, const Gaussian* G4, 
		  const GaussianAux* X13, const GaussianAux* X24, 
		  complex double* s2)
{	
  *s2 += 2*0.25*cvec3mult(X13->sig, X24->sig)* X13->T* X24->T* X13->R* X24->R;
}

static void tb_ps2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const Gaussian* G3, const Gaussian* G4, 
                   const GaussianAux* X13, const GaussianAux* X24, 
                   complex double* s2)
{	
  if (G1->xi == +1 && G2->xi == +1)
    *s2 += 2*0.25*cvec3mult(X13->sig, X24->sig)* X13->T* X24->T* X13->R* X24->R;
}

static void tb_ns2(void* par,
                   const Gaussian* G1, const Gaussian* G2, 
                   const Gaussian* G3, const Gaussian* G4, 
                   const GaussianAux* X13, const GaussianAux* X24, 
                   complex double* s2)
{	
  if (G1->xi == -1 && G2->xi == -1)
    *s2 += 2*0.25*cvec3mult(X13->sig, X24->sig)* X13->T* X24->T* X13->R* X24->R;
}


static void gob_s2(void* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX,
		   complex double* s2, gradGaussian* ds2)
{
  int i;
  complex double s2m;

  s2m = 0.75*X->S;
  
  *s2 += X->T* s2m* X->R;

  for (i=0; i<2; i++)
    ds2->chi[i] += X->T* 0.75*dX->dS.chi[i]* X->R;

  ds2->a += X->T* s2m *dX->dR.a;

  for (i=0; i<3; i++)
    ds2->b[i] += X->T* s2m* dX->dR.b[i];
}

static void gob_ps2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const GaussianAux* X, const gradGaussianAux* dX,
                    complex double* s2, gradGaussian* ds2)
{
  int i;
  complex double s2m;

  if (G1->xi == +1) {
    s2m = 0.75*X->S;
  
    *s2 += X->T* s2m* X->R;

    for (i=0; i<2; i++)
      ds2->chi[i] += X->T* 0.75*dX->dS.chi[i]* X->R;

    ds2->a += X->T* s2m *dX->dR.a;

    for (i=0; i<3; i++)
      ds2->b[i] += X->T* s2m* dX->dR.b[i];
  }
}

static void gob_ns2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const GaussianAux* X, const gradGaussianAux* dX,
                    complex double* s2, gradGaussian* ds2)
{
  int i;
  complex double s2m;

  if (G1->xi == -1) {
    s2m = 0.75*X->S;
  
    *s2 += X->T* s2m* X->R;

    for (i=0; i<2; i++)
      ds2->chi[i] += X->T* 0.75*dX->dS.chi[i]* X->R;

    ds2->a += X->T* s2m *dX->dR.a;

    for (i=0; i<3; i++)
      ds2->b[i] += X->T* s2m* dX->dR.b[i];
  }
}


static void gtb_s2(void* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24,
		   const gradGaussianAux* dX13,
		   complex double* s2, gradGaussian* ds2)
{	
  complex double s2m;
  int i;

  s2m = 0.25*cvec3mult(X13->sig, X24->sig);
  
  *s2 += 2* X13->T*X24->T* s2m *X13->R*X24->R;
  
  for (i=0; i<2; i++)
    ds2->chi[i] += 2* X13->T*X24->T*
       0.25*cvec3mult(dX13->dsig.chi[i], X24->sig)* X13->R* X24->R;

  ds2->a += 2* X13->T*X24->T* s2m* dX13->dR.a*X24->R;

  for (i=0; i<3; i++)
    ds2->b[i] += 2* X13->T*X24->T* s2m* dX13->dR.b[i]*X24->R;  
}

static void gtb_ps2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const Gaussian* G3, const Gaussian* G4, 
                    const GaussianAux* X13, const GaussianAux* X24,
                    const gradGaussianAux* dX13,
                    complex double* s2, gradGaussian* ds2)
{	
  complex double s2m;
  int i;

  if (G1->xi == +1 && G2->xi == +1) {
    s2m = 0.25*cvec3mult(X13->sig, X24->sig);
  
    *s2 += 2* X13->T*X24->T* s2m *X13->R*X24->R;
  
    for (i=0; i<2; i++)
      ds2->chi[i] += 2* X13->T*X24->T*
        0.25*cvec3mult(dX13->dsig.chi[i], X24->sig)* X13->R* X24->R;

    ds2->a += 2* X13->T*X24->T* s2m* dX13->dR.a*X24->R;

    for (i=0; i<3; i++)
      ds2->b[i] += 2* X13->T*X24->T* s2m* dX13->dR.b[i]*X24->R;  
  }
}

static void gtb_ns2(void* par,
                    const Gaussian* G1, const Gaussian* G2, 
                    const Gaussian* G3, const Gaussian* G4, 
                    const GaussianAux* X13, const GaussianAux* X24,
                    const gradGaussianAux* dX13,
                    complex double* s2, gradGaussian* ds2)
{	
  complex double s2m;
  int i;

  if (G1->xi == -1 && G2->xi == -1) {
    s2m = 0.25*cvec3mult(X13->sig, X24->sig);
  
    *s2 += 2* X13->T*X24->T* s2m *X13->R*X24->R;
  
    for (i=0; i<2; i++)
      ds2->chi[i] += 2* X13->T*X24->T*
        0.25*cvec3mult(dX13->dsig.chi[i], X24->sig)* X13->R* X24->R;

    ds2->a += 2* X13->T*X24->T* s2m* dX13->dR.a*X24->R;

    for (i=0; i<3; i++)
      ds2->b[i] += 2* X13->T*X24->T* s2m* dX13->dR.b[i]*X24->R;  
  }
}


void calcConstraintS2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* s2)
{
  double s2one, s2two;
  OneBodyOperator op_ob_s2 = {dim: 1, opt: 1, par: NULL, me: ob_s2};
  TwoBodyOperator op_tb_s2 = {dim: 1, opt: 1, par: NULL, me: tb_s2};


  calcSlaterDetOBME(Q, X, &op_ob_s2, &s2one);
  calcSlaterDetTBME(Q, X, &op_tb_s2, &s2two);

  *s2 = s2one + s2two;
}

void calcConstraintPS2(const SlaterDet* Q, const SlaterDetAux* X, 
                       double* s2)
{
  double s2one, s2two;
  OneBodyOperator op_ob_s2 = {dim: 1, opt: 1, par: NULL, me: ob_ps2};
  TwoBodyOperator op_tb_s2 = {dim: 1, opt: 1, par: NULL, me: tb_ps2};


  calcSlaterDetOBME(Q, X, &op_ob_s2, &s2one);
  calcSlaterDetTBME(Q, X, &op_tb_s2, &s2two);

  *s2 = s2one + s2two;
}

void calcConstraintNS2(const SlaterDet* Q, const SlaterDetAux* X, 
                       double* s2)
{
  double s2one, s2two;
  OneBodyOperator op_ob_s2 = {dim: 1, opt: 1, par: NULL, me: ob_ns2};
  TwoBodyOperator op_tb_s2 = {dim: 1, opt: 1, par: NULL, me: tb_ns2};


  calcSlaterDetOBME(Q, X, &op_ob_s2, &s2one);
  calcSlaterDetTBME(Q, X, &op_tb_s2, &s2two);

  *s2 = s2one + s2two;
}


void calcgradConstraintS2(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* ds2)
{
  gradOneBodyOperator gop_ob_s2 = {opt: 1, par: NULL, me: gob_s2};
  gradTwoBodyOperator gop_tb_s2 = {opt: 1, par: NULL, me: gtb_s2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_s2, ds2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_s2, ds2);
}

void calcgradConstraintPS2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* ds2)
{
  gradOneBodyOperator gop_ob_s2 = {opt: 1, par: NULL, me: gob_ps2};
  gradTwoBodyOperator gop_tb_s2 = {opt: 1, par: NULL, me: gtb_ps2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_s2, ds2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_s2, ds2);
}

void calcgradConstraintNS2(const SlaterDet* Q, const SlaterDetAux* X, 
                           const gradSlaterDetAux* dX,
                           gradSlaterDet* ds2)
{
  gradOneBodyOperator gop_ob_s2 = {opt: 1, par: NULL, me: gob_ns2};
  gradTwoBodyOperator gop_tb_s2 = {opt: 1, par: NULL, me: gtb_ns2};


  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_s2, ds2);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_s2, ds2);
}


void calcConstraintS2od(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X, 
			complex double* s2)
{
  complex double s2one, s2two;
  OneBodyOperator op_ob_s2 = {dim: 1, opt: 1, par: NULL, me: ob_s2};
  TwoBodyOperator op_tb_s2 = {dim: 1, opt: 1, par: NULL, me: tb_s2};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_s2, &s2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_s2, &s2two);

  *s2 = s2one + s2two;
}


void calcgradConstraintS2od(const SlaterDet* Q, const SlaterDet* Qp, 
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* ds2)
{
  gradOneBodyOperator gop_ob_s2 = {opt: 1, par: NULL, me: gob_s2};
  gradTwoBodyOperator gop_tb_s2 = {opt: 1, par: NULL, me: gtb_s2};


  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_s2, ds2);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_s2, ds2);
}


double outputConstraintS2(double val)
{
  return (val);
}

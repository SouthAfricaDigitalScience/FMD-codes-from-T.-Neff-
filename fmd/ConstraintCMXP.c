/**

  \file ConstraintCM.c

  constrain CM position in coordinate and momentum space


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "ConstraintCMXP.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


Constraint ConstraintX = {
  val : 0.0,
  label : "X",
  me : calcConstraintX,
  gradme : calcgradConstraintX,
  output : outputConstraintCM
};

Constraint ConstraintY = {
  val : 0.0,
  label : "Y",
  me : calcConstraintY,
  gradme : calcgradConstraintY,
  output : outputConstraintCM
};

Constraint ConstraintZ = {
  val : 0.0,
  label : "Z",
  me : calcConstraintZ,
  gradme : calcgradConstraintZ,
  output : outputConstraintCM
};

Constraint ConstraintPX = {
  val : 0.0,
  label : "PX",
  me : calcConstraintPX,
  gradme : calcgradConstraintPX,
  output : outputConstraintCM
};

Constraint ConstraintPY = {
  val : 0.0,
  label : "PY",
  me : calcConstraintPY,
  gradme : calcgradConstraintPY,
  output : outputConstraintCM
};

Constraint ConstraintPZ = {
  val : 0.0,
  label : "PZ",
  me : calcConstraintPZ,
  gradme : calcgradConstraintPZ,
  output : outputConstraintCM
};

static void ob_mx(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* x)
{
  *x += mass(G1->xi)*X->rho[0] * X->Q;
}

static void ob_my(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* x)
{
  *x += mass(G1->xi)*X->rho[1] * X->Q;
}

static void ob_mz(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* x)
{
  *x += mass(G1->xi)*X->rho[2] * X->Q;
}


static void ob_px(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double *p)
{	
  *p += X->pi[0] * X->Q;
}

static void ob_py(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double *p)
{	
  *p += X->pi[1] * X->Q;
}

static void ob_pz(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double *p)
{	
  *p += X->pi[2] * X->Q;
}


static void gob_mx(void* dummy,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* mx, gradGaussian* dmx)
{
  complex double mxme;
  int i;

  mxme = mass(G1->xi)*X->rho[0];

  *mx += mxme* X->Q;

  for (i=0; i<2; i++)
    dmx->chi[i] += mxme* dX->dQ.chi[i];

  dmx->a += (mass(G1->xi)*dX->drho.a[0]*X->Q + mxme* dX->dQ.a);
  dmx->b[0] += mass(G1->xi)*dX->drho.b* X->Q;	
  for (i=0; i<3; i++)
    dmx->b[i] += mxme*dX->dQ.b[i];
}

static void gob_my(void* dummy,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* mx, gradGaussian* dmx)
{
  complex double mxme;
  int i;

  mxme = mass(G1->xi)*X->rho[1];

  *mx += mxme* X->Q;

  for (i=0; i<2; i++)
    dmx->chi[i] += mxme* dX->dQ.chi[i];

  dmx->a += (mass(G1->xi)*dX->drho.a[1]*X->Q + mxme* dX->dQ.a);
  dmx->b[1] += mass(G1->xi)*dX->drho.b* X->Q;	
  for (i=0; i<3; i++)
    dmx->b[i] += mxme*dX->dQ.b[i];
}

static void gob_mz(void* dummy,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* mx, gradGaussian* dmx)
{
  complex double mxme;
  int i;

  mxme = mass(G1->xi)*X->rho[2];

  *mx += mxme* X->Q;

  for (i=0; i<2; i++)
    dmx->chi[i] += mxme* dX->dQ.chi[i];

  dmx->a += (mass(G1->xi)*dX->drho.a[2]*X->Q + mxme* dX->dQ.a);
  dmx->b[2] += mass(G1->xi)*dX->drho.b* X->Q;	
  for (i=0; i<3; i++)
    dmx->b[i] += mxme*dX->dQ.b[i];
}

static void gob_px(void* dummy,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* px, gradGaussian* dpx)
{
  complex double pxme;
  int i;

  pxme = X->pi[0];

  *px += pxme* X->Q;

  for (i=0; i<2; i++)
    dpx->chi[i] += pxme* dX->dQ.chi[i];

  dpx->a += (dX->dpi.a[0]*X->Q + pxme* dX->dQ.a);
  dpx->b[0] += dX->dpi.b* X->Q;	
  for (i=0; i<3; i++)
    dpx->b[i] += pxme*dX->dQ.b[i];
}

static void gob_py(void* dummy,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* px, gradGaussian* dpx)
{
  complex double pxme;
  int i;

  pxme = X->pi[1];

  *px += pxme* X->Q;

  for (i=0; i<2; i++)
    dpx->chi[i] += pxme* dX->dQ.chi[i];

  dpx->a += (dX->dpi.a[1]*X->Q + pxme* dX->dQ.a);
  dpx->b[1] += dX->dpi.b* X->Q;	
  for (i=0; i<3; i++)
    dpx->b[i] += pxme*dX->dQ.b[i];
}

static void gob_pz(void* dummy,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX, 
		   complex double* px, gradGaussian* dpx)
{
  complex double pxme;
  int i;

  pxme = X->pi[2];

  *px += pxme* X->Q;

  for (i=0; i<2; i++)
    dpx->chi[i] += pxme* dX->dQ.chi[i];

  dpx->a += (dX->dpi.a[2]*X->Q + pxme* dX->dQ.a);
  dpx->b[2] += dX->dpi.b* X->Q;	
  for (i=0; i<3; i++)
    dpx->b[i] += pxme*dX->dQ.b[i];
}


#define SQR(x) (x)*(x)

void calcConstraintX(const SlaterDet* Q, const SlaterDetAux* X,
		     double* mxcon)
{
  OneBodyOperator op_ob_mx = {dim: 1, opt: 1, par: NULL, me: ob_mx};

  calcSlaterDetOBME(Q, X, &op_ob_mx, mxcon);
}

void calcConstraintY(const SlaterDet* Q, const SlaterDetAux* X,
		     double* mxcon)
{
  OneBodyOperator op_ob_mx = {dim: 1, opt: 1, par: NULL, me: ob_my};

  calcSlaterDetOBME(Q, X, &op_ob_mx, mxcon);
}

void calcConstraintZ(const SlaterDet* Q, const SlaterDetAux* X,
		     double* mxcon)
{
  OneBodyOperator op_ob_mx = {dim: 1, opt: 1, par: NULL, me: ob_mz};

  calcSlaterDetOBME(Q, X, &op_ob_mx, mxcon);
}

void calcConstraintPX(const SlaterDet* Q, const SlaterDetAux* X,
		      double* mxcon)
{
  OneBodyOperator op_ob_mx = {dim: 1, opt: 1, par: NULL, me: ob_px};

  calcSlaterDetOBME(Q, X, &op_ob_mx, mxcon);
}

void calcConstraintPY(const SlaterDet* Q, const SlaterDetAux* X,
		      double* mxcon)
{
  OneBodyOperator op_ob_mx = {dim: 1, opt: 1, par: NULL, me: ob_py};

  calcSlaterDetOBME(Q, X, &op_ob_mx, mxcon);
}

void calcConstraintPZ(const SlaterDet* Q, const SlaterDetAux* X,
		      double* mxcon)
{
  OneBodyOperator op_ob_mx = {dim: 1, opt: 1, par: NULL, me: ob_pz};

  calcSlaterDetOBME(Q, X, &op_ob_mx, mxcon);
}


void calcgradConstraintX(const SlaterDet* Q, const SlaterDetAux* X,
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dmxcon)
{	
  gradOneBodyOperator gop_ob_mx = {opt: 1, par: NULL, me: gob_mx};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_mx, dmxcon);
}

void calcgradConstraintY(const SlaterDet* Q, const SlaterDetAux* X,
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dmxcon)
{	
  gradOneBodyOperator gop_ob_mx = {opt: 1, par: NULL, me: gob_my};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_mx, dmxcon);
}

void calcgradConstraintZ(const SlaterDet* Q, const SlaterDetAux* X,
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dmxcon)
{	
  gradOneBodyOperator gop_ob_mx = {opt: 1, par: NULL, me: gob_mz};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_mx, dmxcon);
}

void calcgradConstraintPX(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dmxcon)
{	
  gradOneBodyOperator gop_ob_mx = {opt: 1, par: NULL, me: gob_px};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_mx, dmxcon);
}

void calcgradConstraintPY(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dmxcon)
{	
  gradOneBodyOperator gop_ob_mx = {opt: 1, par: NULL, me: gob_py};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_mx, dmxcon);
}

void calcgradConstraintPZ(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dmxcon)
{	
  gradOneBodyOperator gop_ob_mx = {opt: 1, par: NULL, me: gob_pz};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_mx, dmxcon);
}


double outputConstraintCM(double val)
{
  return (val);
}

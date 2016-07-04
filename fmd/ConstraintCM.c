/**

  \file ConstraintCM.c

  constrain CM position in coordinate and momentum space


  (c) 2003 Thomas Neff

*/

#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "ConstraintCM.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


Constraint ConstraintCM = {
  val : 0.0,
  label : "CM",
  me : calcConstraintCM,
  gradme : calcgradConstraintCM,
  output : outputConstraintCM
};


typedef struct {
  double xcm[3];
  double pcm[3];
} cmconpara;


static void ob_mx(void* par,
		  const Gaussian*G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double x[3])
{
  int i;
  
  for (i=0; i<3; i++)
    x[i] += mass(G1->xi)*X->rho[i] * X->Q;
}


static void ob_p(void* par,
		 const Gaussian*G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double p[3])
{	
  int i;
  
  for (i=0; i<3; i++)
    p[i] += X->pi[i] * X->Q;
}


static void gob_cmcon(cmconpara* p,
		      const Gaussian* G1, const Gaussian* G2, 
		      const GaussianAux* X, const gradGaussianAux* dX, 
		      complex double* cmcon, gradGaussian* dcmcon)
{
  complex double cmconme;
  int i;

  cmconme = 
    (mass(G1->xi)*
     (X->rho[0]*p->xcm[0] + X->rho[1]*p->xcm[1] + X->rho[2]*p->xcm[2]) +
     (X->pi[0]*p->pcm[0] + X->pi[1]*p->pcm[1] + X->pi[2]*p->pcm[2]));

  *cmcon += 2*cmconme*X->Q;

  for (i=0; i<2; i++)
    dcmcon->chi[i] += 2*cmconme*dX->dQ.chi[i];

  dcmcon->a += 
    2*((mass(G1->xi)*
       (dX->drho.a[0]*p->xcm[0] + dX->drho.a[1]*p->xcm[1] + dX->drho.a[2]*p->xcm[2]) +
	(dX->dpi.a[0]*p->pcm[0] + dX->dpi.a[1]*p->pcm[1] + dX->dpi.a[2]*p->pcm[2]))		
       * X->Q +
       cmconme*dX->dQ.a);

  for (i=0; i<3; i++)
    dcmcon->b[i] += 
      2*((mass(G1->xi)*dX->drho.b* p->xcm[i] + dX->dpi.b* p->pcm[i])*X->Q +
	 cmconme*dX->dQ.b[i]);
}


#define SQR(x) (x)*(x)

void calcConstraintCM(const SlaterDet* Q, const SlaterDetAux* X,
		      double* cmconstraint)
{
  double xcm[3], pcm[3];
  OneBodyOperator op_ob_mx = {dim: 3, opt: 1, par: NULL, me: ob_mx};
  OneBodyOperator op_ob_p = {dim: 3, opt: 1, par: NULL, me: ob_p};

  calcSlaterDetOBME(Q, X, &op_ob_mx, xcm);
  calcSlaterDetOBME(Q, X, &op_ob_p, pcm);

  *cmconstraint = SQR(xcm[0])+SQR(xcm[1])+SQR(xcm[2])+
  		  SQR(pcm[0])+SQR(pcm[1])+SQR(pcm[2]);


  // the donlp2 routine tries to sample unphysical states 
  // from time to time, punish that try

  if (isnan(*cmconstraint))
    *cmconstraint = 100000;

}


void calcgradConstraintCM(const SlaterDet* Q, const SlaterDetAux* X,
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dcmconstraint)
{
  cmconpara p;
  OneBodyOperator op_ob_mx = {dim: 3, opt: 1, par: NULL, me: ob_mx};
  OneBodyOperator op_ob_p = {dim: 3, opt: 1, par: NULL, me: ob_p};
  gradOneBodyOperator gop_ob_cmcon = {opt: 1, par: &p, me: gob_cmcon};
  
  calcSlaterDetOBME(Q, X, &op_ob_mx, p.xcm);
  calcSlaterDetOBME(Q, X, &op_ob_p, p.pcm);
  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_cmcon, dcmconstraint);
}


double outputConstraintCM(double val)
{
  return (val);
}

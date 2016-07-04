/**

  \file Constraintl2.c

  Constrain expectation value of l2


  (c) 2004,2007 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "CenterofMass.h"
#include "gradCMSlaterDet.h"

#include "Constraint.h"
#include "Constraintl2.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


Constraint Constraintl2 = {
  val : 0.0,
  label : "l2",
  me : calcConstraintl2,
  gradme : calcgradConstraintl2,
  output : outputConstraintl2
};

typedef struct {
  double Xcm[3];
  double Vcm[3];
} CMpara;


// l2 is calculated relativ to center of mass 

static void ob_l2(const CMpara* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* l2)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double beta;
  complex double rho2cm = 0.0, pi2cm = 0.0, rhopicm = 0.0;
  complex double rhoxpicm[3];
  int i;
  
  for (i=0; i<3; i++)
    rho2cm += csqr(X->rho[i]-Xcm[i]);
  for (i=0; i<3; i++)
    pi2cm += csqr(X->pi[i]-mass(G1->xi)*Vcm[i]);
  for (i=0; i<3; i++)
    rhopicm += (X->rho[i]-Xcm[i])*(X->pi[i]-mass(G1->xi)*Vcm[i]);
  for (i=0; i<3; i++)
    rhoxpicm[i] = 
      (X->rho[(i+1)%3]-Xcm[(i+1)%3])*(X->pi[(i+2)%3]-mass(G1->xi)*Vcm[(i+2)%3]) -
      (X->rho[(i+2)%3]-Xcm[(i+2)%3])*(X->pi[(i+1)%3]-mass(G1->xi)*Vcm[(i+1)%3]);

  beta = I*X->lambda*(conj(G1->a)-G2->a);

  *l2 += (2*X->lambda*rho2cm + 2*X->alpha*pi2cm - 2*beta*rhopicm +
	 rho2cm* pi2cm - csqr(rhopicm))* X->Q;
}


static void gob_l2(const CMpara* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX,
		   complex double* l2, gradGaussian* dl2)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double l2m;
  complex double beta, dbeta;
  complex double rhocm[3], picm[3];
  complex double rhocm2, picm2, rhopicm;
  complex double rhoxpicm[3];
  gradScalar dl2m, drhocm2, dpicm2, drhopicm;
  int i;

  beta = I*X->lambda*(conj(G1->a)-G2->a);
  dbeta = I*dX->dlambda*(conj(G1->a)-G2->a) + I*X->lambda;

  for (i=0; i<3; i++)
    rhocm[i] = X->rho[i]-Xcm[i];
  rhocm2 = cvec3mult(rhocm, rhocm);
  for (i=0; i<3; i++)
    picm[i] = X->pi[i]-mass(G1->xi)*Vcm[i];
  picm2 = cvec3mult(picm, picm);
  rhopicm = cvec3mult(rhocm, picm);
  cvec3cross(rhocm, picm, rhoxpicm);

  l2m = 2*X->lambda*rhocm2 + 2*X->alpha*picm2 - 2*beta*rhopicm 
	 + rhocm2*picm2 - csqr(rhopicm);
  
  *l2 += X->T* l2m* X->S* X->R;
  
  drhocm2.a = 2*cvec3mult(rhocm, dX->drho.a);
  for (i=0; i<3; i++)
    drhocm2.b[i] = 2*rhocm[i]*dX->drho.b;

  dpicm2.a = 2*cvec3mult(picm, dX->dpi.a);
  for (i=0; i<3; i++)
    dpicm2.b[i] = 2*picm[i]*dX->dpi.b;
    
  drhopicm.a = cvec3mult(dX->drho.a, picm) + cvec3mult(rhocm, dX->dpi.a);
  for (i=0; i<3; i++)
    drhopicm.b[i] = dX->drho.b*picm[i] + rhocm[i]*dX->dpi.b;

  dl2m.a = 2*dX->dlambda*rhocm2 + 2*X->lambda*drhocm2.a +
    2*dX->dalpha*picm2 + 2*X->alpha*dpicm2.a -
    2*dbeta*rhopicm - 2*beta*drhopicm.a +
    drhocm2.a*picm2 + rhocm2*dpicm2.a -
    2*rhopicm*drhopicm.a;
  for (i=0; i<3; i++)
    dl2m.b[i] = 2*X->lambda*drhocm2.b[i] + 2*X->alpha*dpicm2.b[i] - 
      2*beta*drhopicm.b[i] +
      drhocm2.b[i]*picm2 + rhocm2*dpicm2.b[i] - 2*rhopicm*drhopicm.b[i];
    
  for (i=0; i<2; i++)
    dl2->chi[i] += X->T* l2m* dX->dS.chi[i]* X->R;

  dl2->a += X->T* (dl2m.a* X->S*X->R +
		   l2m*X->S* dX->dR.a);

  for (i=0; i<3; i++)
    dl2->b[i] += X->T* (dl2m.b[i]* X->S*X->R +
			l2m*X->S* dX->dR.b[i]);
}

static void gcmob_l2(const CMpara* par,
		    const Gaussian* G1, const Gaussian* G2,
		    const GaussianAux* X,
		    gradCM* dCMl2)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double beta;
  complex double rhocm[3], picm[3];
  complex double rhocm2, picm2, rhopicm;
  complex double rhoxpicm[3];
  gradCM dl2m, drhocm2, dpicm2, drhopicm;
  int i;

  beta = I*X->lambda*(conj(G1->a)-G2->a);

  for (i=0; i<3; i++)
    rhocm[i] = X->rho[i]-Xcm[i];
  rhocm2 = cvec3mult(rhocm, rhocm);
  for (i=0; i<3; i++)
    picm[i] = X->pi[i]-mass(G1->xi)*Vcm[i];
  picm2 = cvec3mult(picm, picm);
  rhopicm = cvec3mult(rhocm, picm);
  cvec3cross(rhocm, picm, rhoxpicm);

  for (i=0; i<3; i++)
    drhocm2.X[i] = -2*rhocm[i];
  for (i=0; i<3; i++)
    dpicm2.V[i] = -2*mass(G1->xi)*picm[i];
  for (i=0; i<3; i++)
    drhopicm.X[i] = -picm[i];
  for (i=0; i<3; i++)
    drhopicm.V[i] = -mass(G1->xi)*rhocm[i];

  for (i=0; i<3; i++)
    dl2m.X[i] = 2*X->lambda*drhocm2.X[i] - 2*beta*drhopicm.X[i] +
      drhocm2.X[i]*picm2 - 2*rhopicm*drhopicm.X[i];
  for (i=0; i<3; i++)
    dl2m.V[i] = 2*X->alpha*dpicm2.V[i] - 2*beta*drhopicm.V[i] +
      rhocm2*dpicm2.V[i] - 2*rhopicm*drhopicm.V[i];

  for (i=0; i<3; i++)
    dCMl2->X[i] += X->T* dl2m.X[i]* X->S*X->R;
  for (i=0; i<3; i++)
    dCMl2->V[i] += X->T* dl2m.V[i]* X->S*X->R;
} 


void calcConstraintl2(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* l2)
{
  CMpara para;

  calcCMPosition(Q, X, para.Xcm);
  calcCMVelocity(Q, X, para.Vcm);
  
  OneBodyOperator op_ob_l2 = {dim: 1, opt: 1, par: &para, me: ob_l2};

  calcSlaterDetOBME(Q, X, &op_ob_l2, l2);
}



void calcgradConstraintl2(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dl2)
{
  CMpara para;
  double l2;

  calcCMPosition(Q, X, para.Xcm);
  calcCMVelocity(Q, X, para.Vcm);

  gradOneBodyOperator gop_ob_l2 = {opt: 1, par: &para, me: gob_l2};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_l2, dl2);

  l2 = dl2->val;

  OneBodyOperator gcmop_ob_l2 = {dim: 6, opt: 1, par: &para, me: gcmob_l2};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_l2, dl2);

  dl2->val = l2;
}


double outputConstraintl2(double val)
{
  return (val);
}

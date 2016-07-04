/**

  \file ConstraintLS.c

  Constrain expectation value of LS


  (c) 2005 Thomas Neff

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
#include "ConstraintLS.h"

#include "numerics/cmath.h"
#include "misc/physics.h"


Constraint ConstraintLS = {
  val : 0.0,
  label : "LS",
  me : calcConstraintLS,
  gradme : calcgradConstraintLS,
  output : outputConstraintLS
};


Constraintod ConstraintLSod = {
  val : 0.0,
  label : "LS",
  dim : 1,
  gradneedsme : 0,
  me : calcConstraintLSod,
  gradme : calcgradConstraintLSod,
  output : outputConstraintLS
};


typedef struct {
  double Xcm[3];
  double Vcm[3];
} CMpara;


// matrix elements identical as in AngularMomentum.c

static void ob_ls(const CMpara* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double* ls)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double rhoxpicm[3];
  int i;
  
  for (i=0; i<3; i++)
    rhoxpicm[i] = 
      (X->rho[(i+1)%3]-Xcm[(i+1)%3])*(X->pi[(i+2)%3]-mass(G1->xi)*Vcm[(i+2)%3]) -
      (X->rho[(i+2)%3]-Xcm[(i+2)%3])*(X->pi[(i+1)%3]-mass(G1->xi)*Vcm[(i+1)%3]);

  *ls += cvec3mult(rhoxpicm, X->sig)* X->T* X->R;
}


static void tb_ls(const CMpara* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const Gaussian* G3, const Gaussian* G4, 
		  const GaussianAux* X13, const GaussianAux* X24, 
		  complex double* ls)
{	
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double rhoxpicm13[3], rhoxpicm24[3];
  int i;

  for (i=0; i<3; i++)
    rhoxpicm13[i] = 
      (X13->rho[(i+1)%3]-Xcm[(i+1)%3])*(X13->pi[(i+2)%3]-mass(G1->xi)*Vcm[(i+2)%3]) -
      (X13->rho[(i+2)%3]-Xcm[(i+2)%3])*(X13->pi[(i+1)%3]-mass(G1->xi)*Vcm[(i+1)%3]);
  for (i=0; i<3; i++)
    rhoxpicm24[i] = 
      (X24->rho[(i+1)%3]-Xcm[(i+1)%3])*(X24->pi[(i+2)%3]-mass(G2->xi)*Vcm[(i+2)%3]) -
      (X24->rho[(i+2)%3]-Xcm[(i+2)%3])*(X24->pi[(i+1)%3]-mass(G2->xi)*Vcm[(i+1)%3]);

  *ls += 0.50*(cvec3mult(rhoxpicm13, X24->sig)*X13->S+
	       cvec3mult(X13->sig, rhoxpicm24)*X24->S)*
    X13->T* X24->T* X13->R* X24->R;
}


static void gob_ls(const CMpara* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const GaussianAux* X, const gradGaussianAux* dX,
		   complex double* ls, gradGaussian* dls)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double lsm;
  complex double rhocm[3], picm[3];
  complex double rhoxpicm[3];
  gradGaussian dlsm;
  int i;

  for (i=0; i<3; i++)
    rhocm[i] = X->rho[i]-Xcm[i];

  for (i=0; i<3; i++)
    picm[i] = X->pi[i]-mass(G1->xi)*Vcm[i];

  cvec3cross(rhocm, picm, rhoxpicm);

  lsm = 0.5*cvec3mult(rhoxpicm, X->sig);
  
  *ls += X->T* 2*lsm *X->R;
  
  for (i=0; i<2; i++)
    dlsm.chi[i] = 0.5*cvec3mult(rhoxpicm, dX->dsig.chi[i]);
  dlsm.a = 0.5*(cvec3spat(dX->drho.a, picm, X->sig)+
		cvec3spat(rhocm, dX->dpi.a, X->sig));
  dlsm.b[0] = 0.5*((dX->drho.b*picm[1]-rhocm[1]*dX->dpi.b)*X->sig[2] - 
		   (dX->drho.b*picm[2]-rhocm[2]*dX->dpi.b)*X->sig[1]);
  dlsm.b[1] = 0.5*((dX->drho.b*picm[2]-rhocm[2]*dX->dpi.b)*X->sig[0] - 
		   (dX->drho.b*picm[0]-rhocm[0]*dX->dpi.b)*X->sig[2]);
  dlsm.b[2] = 0.5*((dX->drho.b*picm[0]-rhocm[0]*dX->dpi.b)*X->sig[1] - 
		   (dX->drho.b*picm[1]-rhocm[1]*dX->dpi.b)*X->sig[0]);

  for (i=0; i<2; i++)
    dls->chi[i] += X->T* 2*dlsm.chi[i]* X->R;

  dls->a += X->T* (2*dlsm.a* X->R + 2*lsm*dX->dR.a);

  for (i=0; i<3; i++)
    dls->b[i] += X->T* (2*dlsm.b[i]* X->R + 2*lsm*dX->dR.b[i]);
}


static void gtb_ls(const CMpara* par,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24,
		   const gradGaussianAux* dX13,
		   complex double* ls, gradGaussian* dls)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double rho13cm[3], rho24cm[3], pi13cm[3], pi24cm[3];
  complex double rhoxpi13cm[3], rhoxpi24cm[3];	
  complex double lsm;
  gradGaussian dlsm;
  int i;

  for (i=0; i<3; i++)
    rho13cm[i] = X13->rho[i]-Xcm[i];
  for (i=0; i<3; i++)
    pi13cm[i] = X13->pi[i]-mass(G1->xi)*Vcm[i];
  cvec3cross(rho13cm, pi13cm, rhoxpi13cm);

  for (i=0; i<3; i++)
    rho24cm[i] = X24->rho[i]-Xcm[i];
  for (i=0; i<3; i++)
    pi24cm[i] = X24->pi[i]-mass(G2->xi)*Vcm[i];
  cvec3cross(rho24cm, pi24cm, rhoxpi24cm);
  
  lsm = 0.25*(cvec3mult(rhoxpi13cm, X24->sig)*X13->S+
	      cvec3mult(X13->sig, rhoxpi24cm)*X24->S);
  
  *ls += 2* X13->T*X24->T* 2*lsm *X13->R*X24->R;
  
  for (i=0; i<2; i++)
    dlsm.chi[i] = 0.25*( 
      cvec3mult(rhoxpi13cm, X24->sig)*dX13->dS.chi[i]+
      cvec3mult(dX13->dsig.chi[i], rhoxpi24cm)*X24->S);
  dlsm.a =
    0.25*cvec3spat(dX13->drho.a, pi13cm, X24->sig)*X13->S+
    0.25*cvec3spat(rho13cm, dX13->dpi.a, X24->sig)*X13->S;

  for (i=0; i<3; i++)
    dlsm.b[i] = 0.25*
      ((dX13->drho.b*pi13cm[(i+1)%3]-rho13cm[(i+1)%3]*dX13->dpi.b)*
       X13->S* X24->sig[(i+2)%3] - 
       (dX13->drho.b*pi13cm[(i+2)%3]-rho13cm[(i+2)%3]*dX13->dpi.b)*
       X13->S* X24->sig[(i+1)%3]);	

  for (i=0; i<3; i++)
    dls->b[i] += 2* X13->T*X24->T*
      (2*dlsm.b[i]* X13->R*X24->R +
       2*lsm *dX13->dR.b[i]*X24->R);  
}


void gcmob_ls(const CMpara* par,
	      const Gaussian* G1, const Gaussian* G2,
	      const GaussianAux* X,
	      gradCM* dCMls)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double rhocm[3], picm[3];
  complex double rhoxpicm[3];
  gradCM dlsm;
  int i;

  for (i=0; i<3; i++)
    rhocm[i] = X->rho[i]-Xcm[i];

  for (i=0; i<3; i++)
    picm[i] = X->pi[i]-mass(G1->xi)*Vcm[i];

  cvec3cross(rhocm, picm, rhoxpicm);

  for (i=0; i<3; i++)
    dlsm.X[i] = -0.5*(picm[(i+1)%3]*X->sig[(i+2)%3] - 
		      picm[(i+2)%3]*X->sig[(i+1)%3]);
  for (i=0; i<3; i++)
    dlsm.V[i] = 0.5*mass(G1->xi)*(rhocm[(i+1)%3]*X->sig[(i+2)%3] - 
				  rhocm[(i+2)%3]*X->sig[(i+1)%3]);

  for (i=0; i<3; i++)
    dCMls->X[i] += X->T* 2*dlsm.X[i]* X->R;
  for (i=0; i<3; i++)
    dCMls->V[i] += X->T* 2*dlsm.V[i]* X->R;
} 


void gcmtb_ls(const CMpara* par,
	      const Gaussian* G1, const Gaussian* G2,
	      const Gaussian* G3, const Gaussian* G4,
	      const GaussianAux* X13, const GaussianAux* X24,
	      gradCM* dCMls)
{
  double* Xcm = par->Xcm;
  double* Vcm = par->Vcm;
  complex double rho13cm[3], rho24cm[3], pi13cm[3], pi24cm[3];
  gradCM dlsm;
  int i;

  for (i=0; i<3; i++)
    rho13cm[i] = X13->rho[i]-Xcm[i];
  for (i=0; i<3; i++)
    pi13cm[i] = X13->pi[i]-mass(G1->xi)*Vcm[i];

  for (i=0; i<3; i++)
    rho24cm[i] = X24->rho[i]-Xcm[i];
  for (i=0; i<3; i++)
    pi24cm[i] = X24->pi[i]-mass(G2->xi)*Vcm[i];

  for (i=0; i<3; i++)
    dlsm.X[i] = 
      -0.25*(pi13cm[(i+1)%3]*X13->S*X24->sig[(i+2)%3] -
	     pi13cm[(i+2)%3]*X13->S*X24->sig[(i+1)%3])
      -0.25*(X13->sig[(i+2)%3]*pi24cm[(i+1)%3]*X24->S -
	     X13->sig[(i+1)%3]*pi24cm[(i+2)%3]*X24->S);
  for (i=0; i<3; i++)
    dlsm.V[i] = 
      0.25*mass(G1->xi)*(rho13cm[(i+1)%3]*X13->S*X24->sig[(i+2)%3] -
			 rho13cm[(i+2)%3]*X13->S*X24->sig[(i+1)%3]) +
      0.25*mass(G2->xi)*(X13->sig[(i+2)%3]*rho24cm[(i+1)%3]*X24->S -
			 X13->sig[(i+1)%3]*rho24cm[(i+2)%3]*X24->S);

  for (i=0; i<3; i++)
    dCMls->X[i] += 2*X13->T*X24->T* 
      2*dlsm.X[i]* X13->R* X24->R;
  for (i=0; i<3; i++)
    dCMls->V[i] += 2*X13->T*X24->T* 
      2*dlsm.V[i]* X13->R* X24->R;
}


void calcConstraintLS(const SlaterDet* Q, const SlaterDetAux* X, 
		      double* ls)
{
  CMpara para;
  double lsone, lstwo;

  calcCMPosition(Q, X, para.Xcm);
  calcCMVelocity(Q, X, para.Vcm);
  
  OneBodyOperator op_ob_ls = {dim: 1, opt: 1, par: &para, me: ob_ls};
  TwoBodyOperator op_tb_ls = {dim: 1, opt: 1, par: &para, me: tb_ls};

  calcSlaterDetOBME(Q, X, &op_ob_ls, &lsone);
  calcSlaterDetTBME(Q, X, &op_tb_ls, &lstwo);

  *ls = lsone + lstwo;
}



void calcgradConstraintLS(const SlaterDet* Q, const SlaterDetAux* X, 
			  const gradSlaterDetAux* dX,
			  gradSlaterDet* dls)
{
  CMpara para;
  double ls;

  calcCMPosition(Q, X, para.Xcm);
  calcCMVelocity(Q, X, para.Vcm);

  gradOneBodyOperator gop_ob_ls = {opt: 1, par: &para, me: gob_ls};
  gradTwoBodyOperator gop_tb_ls = {opt: 1, par: &para, me: gtb_ls};

  calcgradSlaterDetOBME(Q, X, dX, &gop_ob_ls, dls);
  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_ls, dls);

  ls = dls->val;

  OneBodyOperator gcmop_ob_ls = {dim: 6, opt: 1, par: &para, me: gcmob_ls};
  TwoBodyOperator gcmop_tb_ls = {dim: 6, opt: 1, par: &para, me: gcmtb_ls};

  calcgradCMSlaterDetOBME(Q, X, dX, &gcmop_ob_ls, dls);
  calcgradCMSlaterDetTBME(Q, X, dX, &gcmop_tb_ls, dls);

  dls->val = ls;
}


void calcConstraintLSod(const SlaterDet* Q, const SlaterDet* Qp, 
			const SlaterDetAux* X, 
			complex double* ls)
{
  CMpara para = { { 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
  complex double lsone, lstwo;

  OneBodyOperator op_ob_ls = {dim: 1, opt: 1, par: &para, me: ob_ls};
  TwoBodyOperator op_tb_ls = {dim: 1, opt: 1, par: &para, me: tb_ls};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_ls, &lsone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_ls, &lstwo);

  *ls = lsone + lstwo;
}


void calcgradConstraintLSod(const SlaterDet* Q, const SlaterDet* Qp, 	
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    const double* dummy,
			    gradSlaterDet* dls)
{
  CMpara para = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };

  gradOneBodyOperator gop_ob_ls = {opt: 1, par: &para, me: gob_ls};
  gradTwoBodyOperator gop_tb_ls = {opt: 1, par: &para, me: gtb_ls};

  calcgradSlaterDetOBMEod(Q, Qp, X, dX, &gop_ob_ls, dls);
  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_ls, dls);
}


double outputConstraintLS(double val)
{
  return (val);
}

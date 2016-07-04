/**

  \file AngularMomenta.c

  calculate angular momenta


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "AngularMomenta.h"

#include "numerics/cmath.h"


static void ob_angmom(void* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const GaussianAux* X, complex double ang[3])
{
  complex double l2, s2, ls;

  l2 = (2*X->lambda*
	(conj(G1->b[0])*G2->b[0]+conj(G1->b[1])*G2->b[1]+conj(G1->b[2])*G2->b[2]) +
	X->rho2*X->pi2 - csqr(X->rhopi))*X->Q;
  s2 = 0.75*X->Q;
  ls = 0.5*cvec3mult(X->rhoxpi, X->sig)* X->T* X->R;
  
  ang[0] += l2;
  ang[1] += s2;
  ang[2] += l2 + s2 + 2*ls;
}


static void tb_angmom(void* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const Gaussian* G3, const Gaussian* G4, 
		      const GaussianAux* X13, const GaussianAux* X24, 
		      complex double ang[3])
{	
  complex double l2, s2, ls;

  l2 = cvec3mult(X13->rhoxpi, X24->rhoxpi)* X13->Q* X24->Q;
  s2 = 0.25*cvec3mult(X13->sig, X24->sig)* X13->T* X24->T* X13->R* X24->R;
  ls = 0.25*(cvec3mult(X13->rhoxpi, X24->sig)*X13->S+
	     cvec3mult(X13->sig, X24->rhoxpi)*X24->S)*
       X13->T* X24->T* X13->R* X24->R;

  ang[0] += 2*l2;
  ang[1] += 2*s2;
  ang[2] += 2*(l2 + s2 + 2*ls);
}

static void ob_l(void* par,
		 const Gaussian* G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double l[3])
{
  int i;

  for (i=0; i<3; i++)
    l[i] += X->T* X->rhoxpi[i]*X->S* X->R;
}

static void ob_s(void* par,
		 const Gaussian* G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double s[3])
{
  int i;

  for (i=0; i<3; i++)
    s[i] += X->T* 0.5*X->sig[i]* X->R;
}

static void ob_j(void* par,
		 const Gaussian* G1, const Gaussian* G2, 
		 const GaussianAux* X, complex double j[3])
{
  int i;

  for (i=0; i<3; i++)
    j[i] += X->T*(X->rhoxpi[i]*X->S + 0.5*X->sig[i])* X->R;
}


static void ob_j2(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const GaussianAux* X, complex double j2[3])
{
  complex double l2, s2, ls;
  int i;

  for (i=0; i<3; i++) {
    l2 = (X->lambda* (conj(G1->b[(i+1)%3])*G2->b[(i+1)%3]+ 
		      conj(G1->b[(i+2)%3])*G2->b[(i+2)%3]) + 
	  csqr(X->rhoxpi[i]))*X->Q;
    s2 = 0.25*X->Q;
    ls = 0.5*X->rhoxpi[i]* X->sig[i]* X->T* X->R;
  
    j2[i] += l2 + 2*ls + s2;
  }
}


static void tb_j2(void* par,
		  const Gaussian* G1, const Gaussian* G2, 
		  const Gaussian* G3, const Gaussian* G4, 
		  const GaussianAux* X13, const GaussianAux* X24, 
		  complex double j2[3])
{	
  complex double l2, s2, ls;
  int i;

  for (i=0; i<3; i++) {
    l2 = X13->rhoxpi[i]* X24->rhoxpi[i]* X13->Q* X24->Q;
    s2 = 0.25*X13->sig[i]* X24->sig[i]* X13->T* X24->T* X13->R* X24->R;
    ls = 0.25*(X13->rhoxpi[i]* X24->sig[i]* X13->S+
	       X13->sig[i]* X24->rhoxpi[i]* X24->S)*
      X13->T* X24->T* X13->R* X24->R;

    j2[i] += 2*(l2 + s2 + 2*ls);
  }
}


void calcAngularMomenta(const SlaterDet* Q, const SlaterDetAux* X, 
			double* l2, double* s2, double* j2)
{
  double angone[3], angtwo[3];
  OneBodyOperator op_ob_angmom = {dim: 3, opt: 1, par: NULL, me: ob_angmom};
  TwoBodyOperator op_tb_angmom = {dim: 3, opt: 1, par: NULL, me: tb_angmom};


  calcSlaterDetOBME(Q, X, &op_ob_angmom, angone);
  calcSlaterDetTBME(Q, X, &op_tb_angmom, angtwo);

  *l2 = angone[0] + angtwo[0];
  *s2 = angone[1] + angtwo[1];
  *j2 = angone[2] + angtwo[2];
}


void calcAngularMomentaod(const SlaterDet* Q, const SlaterDet* Qp,
			  const SlaterDetAux* X, 
			  complex double* l2, complex double* s2, 
			  complex double* j2)
{
  complex double angone[3], angtwo[3];
  OneBodyOperator op_ob_angmom = {dim: 3, opt: 1, par: NULL, me: ob_angmom};
  TwoBodyOperator op_tb_angmom = {dim: 3, opt: 1, par: NULL, me: tb_angmom};


  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_angmom, angone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_angmom, angtwo);

  *l2 = angone[0] + angtwo[0];
  *s2 = angone[1] + angtwo[1];
  *j2 = angone[2] + angtwo[2];
}


void calcl2HF(const SlaterDet* Q, const SlaterDetAux* X,
	      void* mes)
{
  OneBodyOperator op_ob_angmom = {dim: 3, opt: 1, par: NULL, me: ob_angmom};
  int A=Q->A;
  complex double angmommes[A*A][3];
  complex double *l2mes = mes;

  calcSlaterDetOBHFMEs(Q, X, &op_ob_angmom, angmommes);

  int k,m;
  for (m=0; m<A; m++)
    for (k=0; k<A; k++)
      l2mes[k+m*A] = angmommes[k+m*A][0];
}


void calcj2HF(const SlaterDet* Q, const SlaterDetAux* X,
	      void* mes)
{
  OneBodyOperator op_ob_angmom = {dim: 3, opt: 1, par: NULL, me: ob_angmom};
  int A=Q->A;
  complex double angmommes[A*A][3];
  complex double *j2mes = mes;

  calcSlaterDetOBHFMEs(Q, X, &op_ob_angmom, angmommes);

  int k,m;
  for (m=0; m<A; m++)
    for (k=0; k<A; k++)
      j2mes[k+m*A] = angmommes[k+m*A][2];
}


void calcJ2(const SlaterDet* Q, const SlaterDetAux* X, 
	    double j2[3])
{
  double j2one[3], j2two[3];
  int i;
  OneBodyOperator op_ob_j2 = {dim: 3, opt: 1, par: NULL, me: ob_j2};
  TwoBodyOperator op_tb_j2 = {dim: 3, opt: 1, par: NULL, me: tb_j2};

  calcSlaterDetOBME(Q, X, &op_ob_j2, j2one);
  calcSlaterDetTBME(Q, X, &op_tb_j2, j2two);

  for (i=0; i<3; i++)
    j2[i] = j2one[i] + j2two[i];
}


void calcJ2od(const SlaterDet* Q, const SlaterDet* Qp,
	      const SlaterDetAux* X, 
	      complex double j2[3])
{
  complex double j2one[3], j2two[3];
  int i;
  OneBodyOperator op_ob_j2 = {dim: 3, opt: 1, par: NULL, me: ob_j2};
  TwoBodyOperator op_tb_j2 = {dim: 3, opt: 1, par: NULL, me: tb_j2};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_j2, j2one);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_j2, j2two);

  for (i=0; i<3; i++)
    j2[i] = j2one[i] + j2two[i];
}


void calcL(const SlaterDet* Q, const SlaterDetAux* X,
	   double l[3])
{
  OneBodyOperator op_ob_l = {dim: 3, opt: 1, par: NULL, me: ob_l};

  calcSlaterDetOBME(Q, X, &op_ob_l, l);
}

void calcS(const SlaterDet* Q, const SlaterDetAux* X,
	   double s[3])
{
  OneBodyOperator op_ob_s = {dim: 3, opt: 1, par: NULL, me: ob_s};

  calcSlaterDetOBME(Q, X, &op_ob_s, s);
}

void calcJ(const SlaterDet* Q, const SlaterDetAux* X,
	   double j[3])
{
  OneBodyOperator op_ob_j = {dim: 3, opt: 1, par: NULL, me: ob_j};

  calcSlaterDetOBME(Q, X, &op_ob_j, j);
}


void calcLod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X,
	     complex double l[3])
{
  OneBodyOperator op_ob_l = {dim: 3, opt: 1, par: NULL, me: ob_l};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_l, l);
}

void calcSod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X,
	     complex double s[3])
{
  OneBodyOperator op_ob_s = {dim: 3, opt: 1, par: NULL, me: ob_s};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_s, s);
}

void calcJod(const SlaterDet* Q, const SlaterDet* Qp,
	     const SlaterDetAux* X,
	     complex double j[3])
{
  OneBodyOperator op_ob_j = {dim: 3, opt: 1, par: NULL, me: ob_j};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_j, j);
}


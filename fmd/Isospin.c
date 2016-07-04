/**

  \file Isospin.c

  calculate T^2


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Isospin.h"

#include "numerics/cmath.h"


static void ob_t3(void* par,
		  const Gaussian* G1, const Gaussian* G2,
		  const GaussianAux* X, complex double* t3)
{
  *t3 += (G1->xi == +1 ? 0.5 : -0.5)* X->Q;
}


static void ob_isospin(void* par,
		       const Gaussian* G1, const Gaussian* G2, 
		       const GaussianAux* X, complex double* t2)
{
  *t2 += 0.75*X->Q;
}


static void tb_isospin(void* par,
		       const Gaussian* G1, const Gaussian* G2, 
		       const Gaussian* G3, const Gaussian* G4, 
		       const GaussianAux* X13, const GaussianAux* X24, 
		       complex double* t2)
{	
  *t2 += 0.5*(((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
	      (G1->xi*G4->xi+G2->xi*G3->xi)/2)*
         X13->S* X24->S* X13->R* X24->R;
}


void calcIsospin(const SlaterDet* Q, const SlaterDetAux* X, 
		 double* t2)
{
  double isoone, isotwo;
  OneBodyOperator op_ob_isospin = {dim: 1, opt: 1, par: NULL, me: ob_isospin};
  TwoBodyOperator op_tb_isospin = {dim: 1, opt: 0, par: NULL, me: tb_isospin};


  calcSlaterDetOBME(Q, X, &op_ob_isospin, &isoone);
  calcSlaterDetTBME(Q, X, &op_tb_isospin, &isotwo);

  *t2 = isoone+isotwo;
}


void calcIsospinod(const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double* t2)
{
  complex double isoone, isotwo;
  OneBodyOperator op_ob_isospin = {dim: 1, opt: 1, par: NULL, me: ob_isospin};
  TwoBodyOperator op_tb_isospin = {dim: 1, opt: 0, par: NULL, me: tb_isospin};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_isospin, &isoone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_isospin, &isotwo);

  *t2 = isoone+isotwo;
}


void calct3HF(const SlaterDet* Q, const SlaterDetAux* X,
	      void* mes)
{
  OneBodyOperator op_ob_t3 = {dim: 1, opt: 1, par: NULL, me: ob_t3};

  calcSlaterDetOBHFMEs(Q, X, &op_ob_t3, mes);
}

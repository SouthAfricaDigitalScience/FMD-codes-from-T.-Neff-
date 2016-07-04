/**

  \file Parity.c

  calculate parity


  (c) 2003 Thomas Neff

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Parity.h"


// these are used as workspace

static SlaterDet Qpp;
static SlaterDetAux Xpp;


void ob_parity(void* par,
	       const Gaussian* G1, const Gaussian* G2,
	       const GaussianAux* X, complex double* pi)
{
  Gaussian G2p = *G2;
  GaussianAux Xp;

  invertGaussian(&G2p);
  calcGaussianAux(G1, &G2p, &Xp);

  *pi += Xp.Q;
}


void calcParity(const SlaterDet* Q, const SlaterDetAux* X, 
		double* pi)
{
  // do we have to initialize the workspace first ?
  if (Qpp.A == 0) {
    initSlaterDet(Q, &Qpp);
    initSlaterDetAux(&Qpp, &Xpp);
  }

  copySlaterDet(Q, &Qpp);
  calcSlaterDetAuxod(Q, &Qpp, &Xpp);

  double norm = Xpp.ovlap;

  invertSlaterDet(&Qpp);
  calcSlaterDetAuxod(Q, &Qpp, &Xpp);

  *pi = Xpp.ovlap/norm;
}


void calcParityod(const SlaterDet* Q, const SlaterDet* Qp,
		  const SlaterDetAux* X, 
		  complex double* pi)
{
  // do we have to initialize the workspace first ?
  if (Qpp.A == 0) {
    initSlaterDet(Qp, &Qpp);
    initSlaterDetAux(&Qpp, &Xpp);
  }

  copySlaterDet(Qp, &Qpp);
  invertSlaterDet(&Qpp);

  calcSlaterDetAuxod(Q, &Qpp, &Xpp);

  *pi = Xpp.ovlap;
}


void calcparHF(const SlaterDet* Q, const SlaterDetAux* X, void* mes)
{
  OneBodyOperator op_ob_par = {dim: 1, opt: 1, par: NULL, me: ob_parity};

  calcSlaterDetOBHFMEs(Q, X, &op_ob_par, mes);
}


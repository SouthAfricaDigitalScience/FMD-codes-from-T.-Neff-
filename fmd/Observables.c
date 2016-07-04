/**

  \file Observables.c

  calculate observables for Slater determinants


  (c) 2003, 2005 Thomas Neff

*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "SlaterDet.h"
#include "Interaction.h"
#include "KineticEnergy.h"
#include "CenterofMass.h"
#include "Potential.h"
#include "Radii.h"
#include "Quadrupole.h"
#include "AngularMomenta.h"
#include "Isospin.h"
#include "Parity.h"

#include "misc/physics.h"

#include "Observables.h"



void calcObservables(const Interaction* Int,
		     const SlaterDet* Q,
		     Observables* obs)
{
  SlaterDetAux X;
  initSlaterDetAux(Q, &X);
  calcSlaterDetAux(Q, &X);

  calcT(Q, &X, &obs->t);
  calcTCM(Q, &X, &obs->tcm);
  calcPotential(Int, Q, &X, obs->v);

  if (Int->cm) obs->t -= obs->tcm;
  obs->h = obs->t + obs->v[0];

  calcRadii2(Q, &X, &obs->r2m, &obs->r2p, &obs->r2n);
  calcAngularMomenta(Q, &X, &obs->l2, &obs->s2, &obs->j2);
  calcParity(Q, &X, &obs->pi);
  calcIsospin(Q, &X, &obs->t2);

  freeSlaterDetAux(&X);
}


void calcObservablesParity(const Interaction* Int, int parity,
			   const SlaterDet* Q,
			   Observables* obs)
{
  Observablesod obsd, obsp;

  SlaterDet Qp;
  cloneSlaterDet(Q, &Qp);
  invertSlaterDet(&Qp);

  SlaterDetAux X;
  initSlaterDetAux(Q, &X);
  calcSlaterDetAuxod(Q, Q, &X);
  calcObservablesod(Int, Q, Q, &X, &obsd); 

  calcSlaterDetAuxod(Q, &Qp, &X);
  calcObservablesod(Int, Q, &Qp, &X, &obsp);

  double n = obsd.n + parity* obsp.n;

  obs->h = (obsd.h + parity* obsp.h)/n;
  obs->tcm = (obsd.tcm + parity* obsp.tcm)/n;
  obs->r2m = (obsd.r2m + parity* obsp.r2m)/n;
  obs->r2p = (obsd.r2p + parity* obsp.r2p)/n;
  obs->r2n = (obsd.r2n + parity* obsp.r2n)/n;
  obs->l2 = (obsd.l2 + parity* obsp.l2)/n;
  obs->s2 = (obsd.s2 + parity* obsp.s2)/n;
  obs->j2 = (obsd.j2 + parity* obsp.j2)/n;
  obs->pi = (obsd.pi + parity* obsp.pi)/n;
  obs->t2 = (obsd.t2 + parity* obsp.t2)/n;
  obs->t = (obsd.t + parity* obsp.t)/n;
 
  int i;
  for (i=0; i<Int->n; i++)
    obs->v[i] = (obsd.v[i] + parity* obsp.v[i])/n;

  freeSlaterDet(&Qp);
  freeSlaterDetAux(&X);
}


void fprintObservables(FILE *fp, 
		       const Interaction* Int, const SlaterDet* Q,
		       const Observables* obs)
{
  fprintf(fp, "\nH = %9.3f MeV,   T = %9.3f MeV,   V = %9.3f MeV\n",
	  hbc*obs->h, hbc*obs->t, hbc*obs->v[0]);

  fprintf(fp, "\n       Tcm =  %9.3f MeV\n\n", hbc*obs->tcm);

  int i;
  for (i=1; i<Int->n; i++)
    fprintf(fp, "%10s = %10.3f MeV\n", Int->label[i], hbc*obs->v[i]);

  fprintf(fp, "\nrm = %4.3f fm,   rp = %4.3f fm,   rn = %4.3f fm\n",
	  sqrt(obs->r2m), sqrt(obs->r2p), sqrt(obs->r2n));
  fprintf(fp, "rcharge = %4.3f fm\n", sqrt(r2charge(obs->r2p, Q->N, Q->Z)));

  fprintf(fp, "\nL2 = %6.3f,   S2 = %6.3f,   J2 = %6.3f\n",
	  obs->l2, obs->s2, obs->j2);
  fprintf(fp, "pi = %6.3f\n", obs->pi);
  fprintf(fp, "\nT2 = %6.3f\n", obs->t2);
}



void calcObservablesod(const Interaction* Int,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       const SlaterDetAux* X, 
		       Observablesod* obs)
{
  obs->n = X->ovlap;

  calcTod(Q, Qp, X, &obs->t);
  calcTCMod(Q, Qp, X, &obs->tcm);
  calcPotentialod(Int, Q, Qp, X, obs->v);

  if (Int->cm) obs->t -= obs->tcm;
  obs->h = obs->t + obs->v[0];

  calcRadii2od(Q, Qp, X, &obs->r2m, &obs->r2p, &obs->r2n);
  calcAngularMomentaod(Q, Qp, X, &obs->l2, &obs->s2, &obs->j2);
  calcParityod(Q, Qp, X, &obs->pi);
  calcIsospinod(Q, Qp, X, &obs->t2);
}


void scaleObservablesod(const Interaction* Int, Observablesod* obs)
{
  int i;
  complex double sumV = 0.0;

  for (i=1; i<Int->n; i++) {
    obs->v[i] *= Int->scale[i];
    sumV += obs->v[i];
  }
  obs->v[0] = sumV;
  obs->h = obs->t + obs->v[0];
}

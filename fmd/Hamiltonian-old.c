/**

  \file Hamiltonian.c

  calculate matrix element of Hamiltonian


  (c) 2003 Thomas Neff

*/


#include <complex.h>

#include "SlaterDet.h"
#include "Interaction.h"
#include "KineticEnergy.h"
#include "CenterofMass.h"
#include "Potential.h"
#include "Hamiltonian.h"


void calcHamiltonian(const Interaction *Int,
		     const SlaterDet* Q, const SlaterDetAux* X,
		     double* h)
{
  double v[Int->n];
  double t, tcm = 0.0;

  calcT(Q, X, &t); 
  if (Int->cm) 
    calcTCM(Q, X, &tcm);
  calcPotential(Int, Q, X, v);

  *h = t - tcm + v[0];
}


void calcHamiltonianod(const Interaction *Int,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       const SlaterDetAux* X,
		       complex double* h)
{
  complex double v[Int->n];
  complex double t, tcm = 0.0;

  calcTod(Q, Qp, X, &t);
  if (Int->cm)
    calcTCMod(Q, Qp, X, &tcm);
  calcPotentialod(Int, Q, Qp, X, v);
  
  *h = t - tcm + v[0];
}


// Tcm ?
void calcHamiltonianHF(const Interaction* Int,
		       const SlaterDet* Q, const SlaterDetAux* X,
		       void* mes)
{
  int A=Q->A;
  complex double t[A*A];
  complex double tcm[A*A];
  complex double v[A*A];

  calcTHF(Q, X, t);
  if (Int->cm)
    calcTCMHF(Q, X, tcm);
  calcPotentialHF(Int, Q, X, v);

  complex double *h = mes;
  int k,m;

  for (m=0; m<A; m++)
    for (k=0; k<A; k++) {
      h[k+m*A] = t[k+m*A] + v[k+m*A];
      if (Int->cm)
	h[k+m*A] -= tcm[k+m*A];
    }
}

  

/**

  \file gradHamiltonian.c

  calculate gradients of Hamiltonian matrix elements


  (c) 2003 Thomas Neff

*/


#include "gradKineticEnergy.h"
#include "gradCenterofMass.h"
#include "gradPotential.h"

#include "gradHamiltonian.h"


void calcgradHamiltonian(const Interaction *P,
			 const SlaterDet* Q, const SlaterDetAux* X, 
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dh)
{
  dh->ngauss = Q->ngauss;
  dh->val = 0.0;
  zerogradSlaterDet(dh);

  calcgradT(Q, X, dX, dh);
  if (P->cm)
    calcgradTCM(Q, X, dX, dh);
  calcgradPotential(P, Q, X, dX, dh);
}


void calcgradHamiltonianod(const Interaction *P,
			   const SlaterDet* Q, const SlaterDet* Qp,
			   const SlaterDetAux* X, 
			   const gradSlaterDetAux* dX,
			   gradSlaterDet* dh)
{
  dh->ngauss = Q->ngauss;
  dh->val = 0.0;
  zerogradSlaterDet(dh);

  calcgradTod(Q, Qp, X, dX, dh);
  if (P->cm)
    calcgradTCMod(Q, Qp, X, dX, dh);
  calcgradPotentialod(P, Q, Qp, X, dX, dh);
}

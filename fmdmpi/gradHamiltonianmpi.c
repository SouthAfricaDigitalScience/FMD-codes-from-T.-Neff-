/**

  \file gradHamiltonianmpi.c

  calculate gradients of Hamiltonian matrix elements


  (c) 2003 Thomas Neff

*/

#include <complex.h>
#include <mpi.h>

#include "fmd/SlaterDet.h"
#include "fmd/Interaction.h"
#include "fmd/gradGaussian.h"
#include "fmd/gradKineticEnergy.h"
#include "fmd/gradCenterofMass.h"

#include "Communication.h"
#include "gradHamiltonianmpi.h"


void calcgradHamiltonianmpi(const Interaction *P,
			    const SlaterDet* Q, const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* dh)
{
  dh->val = 0.0;
  zerogradSlaterDet(dh);

  calcgradT(Q, X, dX, dh);
  if (P->cm)
    calcgradTCM(Q, X, dX, dh);
  calcgradPotentialmpi(P, Q, X, dX, dh);
}


void calcgradHamiltonianodmpi(const Interaction *P,
			      const SlaterDet* Q, const SlaterDet* Qp,
			      const SlaterDetAux* X, 
			      const gradSlaterDetAux* dX,
			      gradSlaterDet* dh)
{
  dh->val = 0.0;
  zerogradSlaterDet(dh);

  calcgradTod(Q, Qp, X, dX, dh);
  if (P->cm)
    calcgradTCMod(Q, Qp, X, dX, dh);
  calcgradPotentialodmpi(P, Q, Qp, X, dX, dh);
}


#define SQR(x) ((x)*(x))

void calcgradPotentialmpi(const Interaction *P,
                          const SlaterDet* Q, const SlaterDetAux* X, 
                          const gradSlaterDetAux* dX,
                          gradSlaterDet* dv)
{
  int task = TASKGRADPOTENTIAL;
  int i;
  int A=Q->A;
  int todo=A*A;
  int ngauss=Q->ngauss;
  int k[2]={0,0}; 
  int done=0; int rank;
  gradSlaterDet dvrowcol;
  MPI_Status status;
  
  allocategradSlaterDet(&dvrowcol, A);

  BroadcastTask(&task);
  BroadcastSlaterDet(Q);
  // BroadcastSlaterDetAux(Q, X);
  // BroadcastgradSlaterDetAux(Q, dX);

  // supply slaves with work
  for (rank=1; rank<mpisize; rank++) {

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=0;
    }
  }

  while (1) {
    MPI_Recv(&dvrowcol.val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, 
	     TAGGRADPOTENTIALVAL,
             MPI_COMM_WORLD, &status);
    MPI_Recv(dvrowcol.gradval, ngauss*sizeof(gradGaussian), 
             MPI_BYTE, MPI_ANY_SOURCE, TAGGRADPOTENTIALGRAD, 
             MPI_COMM_WORLD, &status); 
    rank = status.MPI_SOURCE;

    dv->val += dvrowcol.val;
    for (i=0; i<ngauss; i++)
      addmulttogradGaussian(&dv->gradval[i], &dvrowcol.gradval[i], 1.0);

    done++;

    // all work distributed ?
    if (k[0] == A)
      break;

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);

    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=0;
    }
  }

  // collect remaining results
  while (done<todo) {
    MPI_Recv(&dvrowcol.val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, 
	     TAGGRADPOTENTIALVAL,
             MPI_COMM_WORLD, &status);
    MPI_Recv(dvrowcol.gradval, ngauss*sizeof(gradGaussian), 
             MPI_BYTE, MPI_ANY_SOURCE, TAGGRADPOTENTIALGRAD, 
             MPI_COMM_WORLD, &status); 

    dv->val += dvrowcol.val;
    for (i=0; i<ngauss; i++)
      addmulttogradGaussian(&dv->gradval[i], &dvrowcol.gradval[i], 1.0);
    done++;
  }

  // tell slaves we have finished
  k[0] = -1;
  for (rank=1; rank<mpisize; rank++) {
    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
  }

  freegradSlaterDet(&dvrowcol);
}


void calcgradPotentialodmpi(const Interaction *P,
			    const SlaterDet* Q, const SlaterDet* Qp,
			    const SlaterDetAux* X, 
			    const gradSlaterDetAux* dX,
			    gradSlaterDet* dv)
{
  int task = TASKGRADPOTENTIALOD;
  int i;
  int A=Q->A;
  int todo=A*A;
  int ngauss=Q->ngauss;
  int k[2]={0,0}; 
  int done=0; int rank;
  gradSlaterDet dvrowcol;
  MPI_Status status;
  
  allocategradSlaterDet(&dvrowcol, A);

  BroadcastTask(&task);
  BroadcastSlaterDet(Q);
  BroadcastSlaterDet(Qp);
  // BroadcastSlaterDetAux(Q, X);
  // BroadcastgradSlaterDetAux(Q, dX);

  // supply slaves with work
  for (rank=1; rank<mpisize; rank++) {

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=0;
    }
  }

  while (1) {
    MPI_Recv(&dvrowcol.val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, 
	     TAGGRADPOTENTIALODVAL,
             MPI_COMM_WORLD, &status);
    MPI_Recv(dvrowcol.gradval, ngauss*sizeof(gradGaussian), 
             MPI_BYTE, MPI_ANY_SOURCE, TAGGRADPOTENTIALODGRAD, 
             MPI_COMM_WORLD, &status); 
    rank = status.MPI_SOURCE;

    dv->val += dvrowcol.val;
    for (i=0; i<ngauss; i++)
      addmulttogradGaussian(&dv->gradval[i], &dvrowcol.gradval[i], 1.0);

    done++;

    // all work distributed ?
    if (k[0] == A)
      break;

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);

    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=0;
    }
  }

  // collect remaining results
  while (done<todo) {
    MPI_Recv(&dvrowcol.val, 1, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, 
	     TAGGRADPOTENTIALODVAL,
             MPI_COMM_WORLD, &status);
    MPI_Recv(dvrowcol.gradval, ngauss*sizeof(gradGaussian), 
             MPI_BYTE, MPI_ANY_SOURCE, TAGGRADPOTENTIALODGRAD, 
             MPI_COMM_WORLD, &status); 

    dv->val += dvrowcol.val;
    for (i=0; i<ngauss; i++)
      addmulttogradGaussian(&dv->gradval[i], &dvrowcol.gradval[i], 1.0);
    done++;
  }

  // tell slaves we have finished
  k[0] = -1;
  for (rank=1; rank<mpisize; rank++) {
    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
  }

  freegradSlaterDet(&dvrowcol);
}

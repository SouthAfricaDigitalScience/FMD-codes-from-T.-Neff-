/**

  \file Hamiltonian.c

  calculate matrix element of Hamiltonian


  (c) 2003 Thomas Neff

*/


#include <complex.h>
#include <mpi.h>

#include "fmd/SlaterDet.h"
#include "fmd/Interaction.h"
#include "fmd/KineticEnergy.h"
#include "fmd/CenterofMass.h"

#include "Communication.h"
#include "Hamiltonianmpi.h"

#include "numerics/cmath.h"
#include "numerics/coulomb.h"


void calcHamiltonianmpi(const Interaction *P,
			const SlaterDet* Q, const SlaterDetAux* X,
			double* h)
{
  double v[P->n];
  double t, tcm = 0.0;

  calcT(Q, X, &t); 
  if (P->cm) 
    calcTCM(Q, X, &tcm);
  calcPotentialmpi(P, Q, X, v);

  *h = t - tcm + v[0];
}


void calcHamiltonianodmpi(const Interaction *P,
			  const SlaterDet* Q, const SlaterDet* Qp,
			  const SlaterDetAux* X,
			  complex double* h)
{
  complex double v[P->n];
  complex double t, tcm = 0.0;

  calcTod(Q, Qp, X, &t); 
  if (P->cm) 
    calcTCMod(Q, Qp, X, &tcm);
  calcPotentialodmpi(P, Q, Qp, X, v);

  *h = t - tcm + v[0];
}


#define SQR(x) ((x)*(x))

void calcPotentialmpi(const Interaction *P,
                      const SlaterDet* Q, const SlaterDetAux* X, double v[])
{
  int i;
  int A=Q->A;
  int todo=(A*(A-1))/2;

  int k[2]={0,1}; 
  int done=0; int rank;
  int dim=P->n;
  double vrowcol[dim];
  MPI_Status status;

  for (i=0; i<dim; i++)
    v[i] = 0.0;

  int task = TASKPOTENTIAL;
  BroadcastTask(&task);

  BroadcastSlaterDet(Q);	
  // BroadcastSlaterDetAux(Q, X);

  // supply slaves with work
  for (rank=1; rank<mpisize; rank++) {
    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=k[0]+1;
    }
  }

  while (1) {
    MPI_Recv(vrowcol, dim, MPI_DOUBLE, MPI_ANY_SOURCE, TAGPOTENTIAL,
             MPI_COMM_WORLD, &status);
    rank = status.MPI_SOURCE;

    for (i=0; i<dim; i++)
      v[i] += vrowcol[i];

    done++;

    // all work distributed ?
    if (k[0] == A-1)
      break;

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);

    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=k[0]+1;
    }
  }

  // collect remaining results
  while (done<todo) {
    MPI_Recv(vrowcol, dim, MPI_DOUBLE, MPI_ANY_SOURCE, TAGPOTENTIAL,
             MPI_COMM_WORLD, &status);

    for (i=0; i<dim; i++)
      v[i] += vrowcol[i];
  
    done++;
  }

  // tell slaves we have finished
  k[0] = -1;
  for (rank=1; rank<mpisize; rank++) {
    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
  }
}


void calcPotentialodmpi(const Interaction *P,
			const SlaterDet* Q, const SlaterDet* Qp,
			const SlaterDetAux* X, complex double v[])
{
  int task = TASKPOTENTIALOD;
  int i;
  int A=Q->A;
  int todo=(A*(A-1))/2;
  int k[2]={0,1}; 
  int done=0; int rank;
  int dim=P->n;
  complex double vrowcol[dim];
  MPI_Status status;

  for (i=0; i<dim; i++)
    v[i] = 0.0;

  BroadcastTask(&task);
  BroadcastSlaterDet(Q);
  BroadcastSlaterDet(Qp);
  // BroadcastSlaterDetAux(Q, X);

  // supply slaves with work
  for (rank=1; rank<mpisize; rank++) {

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=k[0]+1;
    }
  }

  while (1) {
    MPI_Recv(vrowcol, dim, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, TAGPOTENTIALOD,
             MPI_COMM_WORLD, &status);
    rank = status.MPI_SOURCE;

    for (i=0; i<dim; i++)
      v[i] += vrowcol[i];

    done++;

    // all work distributed ?
    if (k[0] == A-1)
      break;

    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);

    k[1]++;
    if (k[1] == A) {
      k[0]++; k[1]=k[0]+1;
    }
  }

  // collect remaining results
  while (done<todo) {
    MPI_Recv(vrowcol, dim, MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, TAGPOTENTIALOD,
             MPI_COMM_WORLD, &status);

    for (i=0; i<dim; i++)
      v[i] += vrowcol[i];
  
    done++;
  }

  // tell slaves we have finished
  k[0] = -1;
  for (rank=1; rank<mpisize; rank++) {
    MPI_Send(k, 2, MPI_INT, rank, TAGROWCOL, MPI_COMM_WORLD);
  }
}

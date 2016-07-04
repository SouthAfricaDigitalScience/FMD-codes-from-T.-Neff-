/**

  \file Communication.c

  Communicate using MPI


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/Gaussian.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"

#include "Communication.h"


#define SQR(x) ((x)*(x))

int mpirank;
int mpisize;


void BroadcastA(int* A)
{
  MPI_Bcast(A, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


void BroadcastTask(int* task)
{
  MPI_Bcast(task, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


// slaves overwride sldet
void BroadcastSlaterDet(SlaterDet* Q)
{
  if (mpirank == 0) {
    // send
    MPI_Bcast(Q, sizeof(SlaterDet), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Q->idx, Q->A, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(Q->ng, Q->A, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(Q->G, Q->ngauss*sizeof(Gaussian), 
	      MPI_BYTE, 0, MPI_COMM_WORLD);
  } else {
    int* idx = Q->idx; 
    int* ng = Q->ng;
    Gaussian* G = Q->G;

    MPI_Bcast(Q, sizeof(SlaterDet), MPI_BYTE, 0, MPI_COMM_WORLD);
    Q->idx = idx;
    Q->ng = ng;
    Q->G = G;
    MPI_Bcast(Q->idx, Q->A, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(Q->ng, Q->A, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(Q->G, Q->ngauss*sizeof(Gaussian), 
	      MPI_BYTE, 0, MPI_COMM_WORLD);
    
  }
}


void BroadcastSlaterDetAux(const SlaterDet* Q, SlaterDetAux* X)
{
    MPI_Bcast(X->Gaux, SQR(Q->ngauss)*sizeof(GaussianAux), 
	      MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(X->n, SQR(Q->A)*sizeof(complex double),
	      MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(X->o, SQR(Q->A)*sizeof(complex double),
	      MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&X->ovlap, 1, 
	      MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}


void BroadcastgradSlaterDetAux(const SlaterDet* Q, gradSlaterDetAux* dX)
{
  MPI_Bcast(&dX->ngauss, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dX->dGaux, SQR(dX->ngauss)*sizeof(gradGaussianAux),
	    MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(dX->dno, Q->A*dX->ngauss*sizeof(gradGaussian),
	    MPI_BYTE, 0, MPI_COMM_WORLD);
}


/// slaves allocate space for Interaction
/// labels are not broadcastet
void BroadcastInteraction(Interaction* Int)
{
  if (mpirank == 0) {
    // send
    MPI_Bcast(Int, sizeof(Interaction), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(Int->c, Int->ncomponents*sizeof(InteractionComponent),
	      MPI_BYTE, 0, MPI_COMM_WORLD);
  } else {
    // receive
    MPI_Bcast(Int, sizeof(Interaction), MPI_BYTE, 0, MPI_COMM_WORLD);
    Int->c = malloc(Int->ncomponents*sizeof(InteractionComponent));
    MPI_Bcast(Int->c, Int->ncomponents*sizeof(InteractionComponent),
	      MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  
}


/// generic broadcast for Parameters
/// will only broadcast scalar components 
void BroadcastParameters(void* par, int size)
{
  if (mpirank == 0) {
    // send
    MPI_Bcast(par, size, MPI_BYTE, 0, MPI_COMM_WORLD);
  } else {
    // receive
    MPI_Bcast(par, size, MPI_BYTE, 0, MPI_COMM_WORLD);
  }
}

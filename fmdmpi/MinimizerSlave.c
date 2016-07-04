/**

  \file MinimizerSlave.c

  MPI slave for calculation of Potential and gradPotential matrix elements


  (c) 2003,2004 Thomas Neff

*/


#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/Potential.h"
#include "fmd/gradPotential.h"

#include "Communication.h"
#include "MinimizerSlave.h"


void MinimizerSlave(void)
{
  MPI_Status status;

  Interaction Int;
  int A;
  SlaterDet Q;
  SlaterDetAux X;
  gradSlaterDetAux dX;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastInteraction(&Int);
  BroadcastA(&A);

  double v[Int.n];
  gradSlaterDet dv;
  allocategradSlaterDet(&dv, A);

  allocateSlaterDet(&Q, A);
  allocateSlaterDetAux(&X, A);
  allocategradSlaterDetAux(&dX, A);

  while (1) {

    BroadcastTask(&task);
    if (task == TASKPOTENTIAL) {
      int k[2];

      BroadcastSlaterDet(&Q);
      // BroadcastSlaterDetAux(&Q, &X);
      calcSlaterDetAux(&Q, &X);

      while (1) {
	MPI_Recv(k, 2, MPI_INT, 0, TAGROWCOL, MPI_COMM_WORLD, &status);
	if (k[0]==-1)
	  break;
	calcPotentialrowcol(&Int, &Q, &X, v, k[0], k[1]);
	MPI_Send(v, Int.n, MPI_DOUBLE, 0, TAGPOTENTIAL, MPI_COMM_WORLD);
      }

    } else if (task == TASKGRADPOTENTIAL) {
      int k[2];

      BroadcastSlaterDet(&Q);
      // BroadcastSlaterDetAux(&Q, &X);
      // BroadcastgradSlaterDetAux(&Q, &dX);
      calcSlaterDetAux(&Q, &X);
      calcgradSlaterDetAux(&Q, &X, &dX);


      while (1) {
	MPI_Recv(k, 2, MPI_INT, 0, TAGROWCOL, MPI_COMM_WORLD, &status);
	if (k[0]==-1)
	  break;
	calcgradPotentialrowcol(&Int, &Q, &X, &dX, &dv, k[0], k[1]);
	MPI_Send(&dv.val, 1, MPI_DOUBLE_COMPLEX, 0, TAGGRADPOTENTIALVAL, MPI_COMM_WORLD);
	MPI_Send(dv.gradval, Q.ngauss*sizeof(gradGaussian),
		 MPI_BYTE, 0, TAGGRADPOTENTIALGRAD, MPI_COMM_WORLD);
      }

    } else if (task == TASKFIN) 
      return;
  }
}


void MinimizerSlaveod(void)
{
  MPI_Status status;

  Interaction Int;
  int A;
  SlaterDet Q;
  SlaterDet Qp;
  SlaterDetAux X;
  gradSlaterDetAux dX;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastInteraction(&Int);
  BroadcastA(&A);

  complex double v[Int.n];
  gradSlaterDet dv;
  allocategradSlaterDet(&dv, A);

  allocateSlaterDet(&Q, A);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDetAux(&X, A);
  allocategradSlaterDetAux(&dX, A);

  while (1) {
    int task;

    BroadcastTask(&task);
    if (task == TASKPOTENTIALOD) {
      int k[2];

      BroadcastSlaterDet(&Q);
      BroadcastSlaterDet(&Qp);
      // BroadcastSlaterDetAux(&Q, &X);
      calcSlaterDetAuxod(&Q, &Qp, &X);
      
      while (1) {
	MPI_Recv(k, 2, MPI_INT, 0, TAGROWCOL, MPI_COMM_WORLD, &status);
	if (k[0]==-1)
	  break;
	calcPotentialodrowcol(&Int, &Q, &Qp, &X, v, k[0], k[1]);
	MPI_Send(v, Int.n, MPI_DOUBLE_COMPLEX, 0, TAGPOTENTIALOD, 
		 MPI_COMM_WORLD);
      }

    } else if (task == TASKGRADPOTENTIALOD) {
      int k[2];

      BroadcastSlaterDet(&Q);
      BroadcastSlaterDet(&Qp);
      // BroadcastSlaterDetAux(&Q, &X);
      // BroadcastgradSlaterDetAux(&Q, &dX);
      calcSlaterDetAuxod(&Q, &Qp, &X);
      calcgradSlaterDetAuxod(&Q, &Qp, &X, &dX);

      while (1) {
	MPI_Recv(k, 2, MPI_INT, 0, TAGROWCOL, MPI_COMM_WORLD, &status);
	if (k[0]==-1)
	  break;
	calcgradPotentialodrowcol(&Int, &Q, &Qp, &X, &dX, &dv, k[0], k[1]);
	MPI_Send(&dv.val, 1, MPI_DOUBLE_COMPLEX, 0, TAGGRADPOTENTIALODVAL,
		 MPI_COMM_WORLD);
	MPI_Send(dv.gradval, Q.ngauss*sizeof(gradGaussian),
		 MPI_BYTE, 0, TAGGRADPOTENTIALODGRAD, MPI_COMM_WORLD);
      }

    } else if (task == TASKFIN) 
      return;
  }
}

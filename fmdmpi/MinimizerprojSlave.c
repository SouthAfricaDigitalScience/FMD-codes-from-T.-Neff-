/**

  \file MinimizerprojSlave.c

  MPI slave for calculation of projected matrix elements
  and gradients


  (c) 2006-2011 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Hamiltonian.h"
#include "fmd/ConstraintJ2.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/gradHamiltonian.h"

#include "Communication.h"
#include "MinimizerprojSlave.h"



void MinimizerprojSlave(void)
{
  MPI_Status status;

  Interaction Int;
  int A;
  SlaterDet Q, Qp, Qpp;
  SlaterDetAux X;
  gradSlaterDetAux dX;
  
  complex double h[2], n[2], j2[2];
  gradSlaterDet dh[2], dn[2];

  int task;
  int cmproj=0;

  BroadcastTask(&task);
  if (task != TASKSTART && task != TASKSTARTCM)
    return;

  if (task == TASKSTARTCM)
    cmproj = 1;

  BroadcastInteraction(&Int);
  BroadcastA(&A);

  allocateSlaterDet(&Q, A);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDet(&Qpp, A);
  allocateSlaterDetAux(&X, A);
  allocategradSlaterDetAux(&dX, A);
  allocategradSlaterDet(&dh[0], A);
  allocategradSlaterDet(&dh[1], A);
  allocategradSlaterDet(&dn[0], A);
  allocategradSlaterDet(&dn[1], A);

  while (1) {

    BroadcastTask(&task);
    if (task != TASKHAMILTONIANOD && 
	task != TASKGRADHAMILTONIANOD)
      return;

    BroadcastSlaterDet(&Q);
    BroadcastSlaterDet(&Qp);

    double projpar[6];
    double *angle, *R;

    while (1) {

      if (cmproj) {
	MPI_Recv(projpar, 6, MPI_DOUBLE, 0, TAGPROJECT6, MPI_COMM_WORLD, &status);
	if (projpar[0] < 0.0)
	  break;

	angle = &projpar[0];
	R = &projpar[3];

	copySlaterDet(&Qp, &Qpp);
	moveSlaterDet(&Qpp, R);
	rotateSlaterDet(&Qpp, angle[0], angle[1], angle[2]);
      } else {
	MPI_Recv(projpar, 3, MPI_DOUBLE, 0, TAGPROJECT3, MPI_COMM_WORLD, &status);
	if (projpar[0] < 0.0)
	  break;

	angle = &projpar[0];

	copySlaterDet(&Qp, &Qpp);
	rotateSlaterDet(&Qpp, angle[0], angle[1], angle[2]);
      }

      if (task == TASKHAMILTONIANOD) {
	
	calcSlaterDetAuxod(&Q, &Qpp, &X);
	n[0] = X.ovlap;
	calcHamiltonianod(&Int, &Q, &Qpp, &X, &h[0]);
	calcConstraintJ2od(&Q, &Qpp, &X, &j2[0]);
	
	invertSlaterDet(&Qpp);

	calcSlaterDetAuxod(&Q, &Qpp, &X); 
	n[1] = X.ovlap;
	calcHamiltonianod(&Int, &Q, &Qpp, &X, &h[1]);
	calcConstraintJ2od(&Q, &Qpp, &X, &j2[1]);

	MPI_Send(h, 2, MPI_DOUBLE_COMPLEX, 0, TAGHAMILTONIANOD, MPI_COMM_WORLD);
	MPI_Send(n, 2, MPI_DOUBLE_COMPLEX, 0, TAGOVLAPOD, MPI_COMM_WORLD);
	MPI_Send(j2, 2, MPI_DOUBLE_COMPLEX, 0, TAGANGULARMOMENTUMOD, MPI_COMM_WORLD);
      }

      if (task == TASKGRADHAMILTONIANOD) {

	calcSlaterDetAuxod(&Q, &Qpp, &X);
	calcgradSlaterDetAuxod(&Q, &Qpp, &X, &dX);
	calcgradOvlapod(&Q, &Qpp, &X, &dX, &dn[0]);
	calcgradHamiltonianod(&Int, &Q, &Qpp, &X, &dX, &dh[0]);

	invertSlaterDet(&Qpp);

	calcSlaterDetAuxod(&Q, &Qpp, &X); 
	calcgradSlaterDetAuxod(&Q, &Qpp, &X, &dX);
	calcgradOvlapod(&Q, &Qpp, &X, &dX, &dn[1]);
	calcgradHamiltonianod(&Int, &Q, &Qpp, &X, &dX, &dh[1]);

	MPI_Send(&dh[0].val, 1, MPI_DOUBLE_COMPLEX, 0, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD);
	MPI_Send(dh[0].gradval, Q.ngauss*sizeof(gradGaussian),
		 MPI_BYTE, 0, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD);
	MPI_Send(&dh[1].val, 1, MPI_DOUBLE_COMPLEX, 0, TAGGRADHAMILTONIANODVAL, MPI_COMM_WORLD);
	MPI_Send(dh[1].gradval, Q.ngauss*sizeof(gradGaussian),
		 MPI_BYTE, 0, TAGGRADHAMILTONIANODGRAD, MPI_COMM_WORLD);
	MPI_Send(&dn[0].val, 1, MPI_DOUBLE_COMPLEX, 0, TAGGRADOVLAPODVAL, MPI_COMM_WORLD);
	MPI_Send(dn[0].gradval, Q.ngauss*sizeof(gradGaussian),
		 MPI_BYTE, 0, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD);
	MPI_Send(&dn[1].val, 1, MPI_DOUBLE_COMPLEX, 0, TAGGRADOVLAPODVAL, MPI_COMM_WORLD);
	MPI_Send(dn[1].gradval, Q.ngauss*sizeof(gradGaussian),
		 MPI_BYTE, 0, TAGGRADOVLAPODGRAD, MPI_COMM_WORLD);
      }

    }	

  }

}

/**

  \file ProjectionSlave.c

  MPI slave for calculation of projected Matrix Elements


  (c) 2003,2004,2005,2006 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Observables.h"

#include "Communication.h"
#include "ProjectionSlave.h"


// calculates Overlaps and Observables

void ProjectionSlave(void)
{
  MPI_Status status;

  Interaction Int;
  int A;
  SlaterDet Q, Qp, Qpp;
  SlaterDetAux X;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastInteraction(&Int);
  BroadcastA(&A);

  allocateSlaterDet(&Q, A);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDet(&Qpp, A);
  allocateSlaterDetAux(&X, A);

  while (1) {

    BroadcastTask(&task);
    if (task != TASKPROJECTOVLAPOD && 
	task != TASKPROJECTOBSERVABLESOD)
      return;

    double projpar[6];
    double *angle, *R;

    BroadcastSlaterDet(&Q);
    BroadcastSlaterDet(&Qp);

    while (1) {

      MPI_Recv(projpar, 6, MPI_DOUBLE, 0, TAGPROJECT6, MPI_COMM_WORLD, &status);
      if (projpar[0] < 0.0)
	break;
	
      angle = &projpar[0];
      R = &projpar[3];

      copySlaterDet(&Qp, &Qpp);

      moveSlaterDet(&Qpp, R);
      rotateSlaterDet(&Qpp, angle[0], angle[1], angle[2]);

      if (task == TASKPROJECTOVLAPOD) {
	complex double ovl[2];

	calcSlaterDetAuxod(&Q, &Qpp, &X);
	ovl[0] = X.ovlap;

	invertSlaterDet(&Qpp);

	calcSlaterDetAuxod(&Q, &Qpp, &X); 
	ovl[1] = X.ovlap;

	MPI_Send(ovl, 2*sizeof(complex double), MPI_BYTE, 
		 0, TAGMEOD, MPI_COMM_WORLD);
      }

      if (task == TASKPROJECTOBSERVABLESOD) {
	Observablesod obs[2];

	calcSlaterDetAuxod(&Q, &Qpp, &X); 
	calcObservablesod(&Int, &Q, &Qpp, &X, &obs[0]);

	invertSlaterDet(&Qpp);

	calcSlaterDetAuxod(&Q, &Qpp, &X); 
	calcObservablesod(&Int, &Q, &Qpp, &X, &obs[1]);

	MPI_Send(obs, 2*sizeof(Observablesod), MPI_BYTE, 
		 0, TAGMEOD, MPI_COMM_WORLD);
      }

    }	

  }

}

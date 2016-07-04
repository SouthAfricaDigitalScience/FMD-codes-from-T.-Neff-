/**

  \file FormfactorSlave.c

  MPI slave for calculation of Formfactor matrix elements


  (c) 2008 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/Formfactors.h"

#include "Communication.h"
#include "FormfactorSlave.h"



void FormfactorSlave(void)
{
  MPI_Status status;

  FormfactorPara FfP;
  int A;
  SlaterDet Q, Qp, Qpp;
  SlaterDetAux X;
  complex double* ffme;
  int dimff;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastParameters(&FfP, sizeof(FormfactorPara));
  dimff = 2*FfP.npoints;
  ffme = malloc(2*dimff*sizeof(complex double));

  BroadcastA(&A);

  allocateSlaterDet(&Q, A);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDet(&Qpp, A);
  allocateSlaterDetAux(&X, A);

  while (1) {

    BroadcastTask(&task);
    if (task != TASKPROJECTMONOPOLEFORMFACTOROD)
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

      calcSlaterDetAuxod(&Q, &Qpp, &X); 
      calcMonopoleFormfactorod(&FfP, &Q, &Qpp, &X, &ffme[0]);

      invertSlaterDet(&Qpp);

      calcSlaterDetAuxod(&Q, &Qpp, &X); 
      calcMonopoleFormfactorod(&FfP, &Q, &Qpp, &X, &ffme[dimff]);
      
      MPI_Send(ffme, 2*dimff*sizeof(complex double), MPI_BYTE, 
	       0, TAGMEOD, MPI_COMM_WORLD);

    }	

  }

}

/**

  \file DiagonalDensityHOSlave.c

  MPI slave for calculation of diagonal density matrix in HO basis


  (c) 2008 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/ProjectedDensityMatrixHO.h"

#include "Communication.h"
#include "DiagonalDensityHOSlave.h"



void DiagonalDensityHOSlave(void)
{
  MPI_Status status;

  DensityMatrixHOPar DMpar;
  int A;
  SlaterDet Q, Qp, Qpp;
  SlaterDetAux X;
  complex double* dmme;
  int dimdm;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastParameters(&DMpar, sizeof(DensityMatrixHOPar));
  initHOBasis(DMpar.nmax);

  dimdm = (DMpar.nmax+1)*(DMpar.nmax+2);
  dmme = malloc(2*dimdm*sizeof(complex double));

  BroadcastA(&A);

  allocateSlaterDet(&Q, A);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDet(&Qpp, A);
  allocateSlaterDetAux(&X, A);

  while (1) {

    BroadcastTask(&task);
    if (task != TASKPROJECTEDDIAGONALDENSITYMATRIXHOOD)
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
      calcDiagonalDensityMatrixHOod(&DMpar, &Q, &Qpp, &X, &dmme[0]);

      invertSlaterDet(&Qpp);

      calcSlaterDetAuxod(&Q, &Qpp, &X); 
      calcDiagonalDensityMatrixHOod(&DMpar, &Q, &Qpp, &X, &dmme[dimdm]);
      
      MPI_Send(dmme, 2*dimdm*sizeof(complex double), MPI_BYTE, 
	       0, TAGMEOD, MPI_COMM_WORLD);

    }	

  }

}

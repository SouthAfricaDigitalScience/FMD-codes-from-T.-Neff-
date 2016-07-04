/**

  \file TwoNucleonOvlapsSlave.c

  MPI slave for calculation of projected Matrix Elements


  (c) 2012 Thomas Neff

*/

#include <stdio.h>
#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/TwoNucleonOvlaps.h"

#include "Communication.h"
#include "TwoNucleonOvlapsSlave.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define SQR(x) ((x)*(x))


void BroadcastTwoNucleonOvlapsPara(TwoNucleonOvlapsPara* TNOpar)
{
  // send or receive
  MPI_Bcast(TNOpar, sizeof(TwoNucleonOvlapsPara), MPI_BYTE, 0, MPI_COMM_WORLD);
}    


void TwoNucleonOvlapsSlave(void)
{
  MPI_Status status;

  TwoNucleonOvlapsPara TNOpar;
  int A;
  SlaterDet Q, Qp, Qpp;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastTwoNucleonOvlapsPara(&TNOpar);
  int dim = MAX(DIMSPECT,DIMSPECY)* SQR(TNOpar.npoints);
  BroadcastA(&A);

  complex double specamp[2][dim];

  allocateSlaterDet(&Q, A-2);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDet(&Qpp, A);

  while (1) {

    BroadcastTask(&task);
    if (task != TASKPROJECTTWONUCLEONOVLAPSTOD &&
	task != TASKPROJECTTWONUCLEONOVLAPSYOD)
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

      if (task == TASKPROJECTTWONUCLEONOVLAPSTOD)
	calcTwoNucleonOvlapsTod(&TNOpar, &Q, &Qpp, NULL, specamp[0]);
      else if (task == TASKPROJECTTWONUCLEONOVLAPSYOD)
	calcTwoNucleonOvlapsYod(&TNOpar, &Q, &Qpp, NULL, specamp[0]);

      invertSlaterDet(&Qpp);

      if (task == TASKPROJECTTWONUCLEONOVLAPSTOD)
	calcTwoNucleonOvlapsTod(&TNOpar, &Q, &Qpp, NULL, specamp[1]);
      else if (task == TASKPROJECTTWONUCLEONOVLAPSYOD)
	calcTwoNucleonOvlapsYod(&TNOpar, &Q, &Qpp, NULL, specamp[1]);

      MPI_Send(specamp, 2*dim*sizeof(complex double), MPI_BYTE, 
	       0, TAGMEOD, MPI_COMM_WORLD);
    }	

  }

}

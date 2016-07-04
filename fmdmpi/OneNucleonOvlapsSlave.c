/**

  \file OneNucleonOvlapsSlave.c

  MPI slave for calculation of projected Matrix Elements


  (c) 2012 Thomas Neff

*/

#include <stdio.h>
#include <complex.h>
#include <mpi.h>

#include "fmd/Interaction.h"
#include "fmd/SlaterDet.h"
#include "fmd/OneNucleonOvlaps.h"

#include "Communication.h"
#include "OneNucleonOvlapsSlave.h"

#define SQR(x) ((x)*(x))


void BroadcastOneNucleonOvlapsPara(OneNucleonOvlapsPara* ONOpar)
{
  // send or receive
  MPI_Bcast(ONOpar, sizeof(OneNucleonOvlapsPara), MPI_BYTE, 0, MPI_COMM_WORLD);
}    


void OneNucleonOvlapsSlave(void)
{
  MPI_Status status;

  OneNucleonOvlapsPara ONOpar;
  int A;
  SlaterDet Q, Qp, Qpp;

  int task;

  BroadcastTask(&task);
  if (task != TASKSTART)
    return;

  BroadcastOneNucleonOvlapsPara(&ONOpar);
  int dim = DIMSPEC* ONOpar.npoints;
  BroadcastA(&A);

  complex double specamp[2][dim];

  allocateSlaterDet(&Q, A-1);
  allocateSlaterDet(&Qp, A);
  allocateSlaterDet(&Qpp, A);

  while (1) {

    BroadcastTask(&task);
    if (task != TASKPROJECTONENUCLEONOVLAPSOD)
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

      calcOneNucleonOvlapsod(&ONOpar, &Q, &Qpp, NULL, specamp[0]);

      invertSlaterDet(&Qpp);

      calcOneNucleonOvlapsod(&ONOpar, &Q, &Qpp, NULL, specamp[1]);

      MPI_Send(specamp, 2*dim*sizeof(complex double), MPI_BYTE, 
	       0, TAGMEOD, MPI_COMM_WORLD);
    }	

  }

}

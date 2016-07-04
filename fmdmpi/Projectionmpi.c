/**

  \file Projectionmpi.c

  angular momentum projection of matrix elements


  (c) 2003, 2004, 2005 Thomas Neff

*/

// 11 December 2014: KRH: Commented out sections doing Normalisation
// Put back the line   weight[processor] = 0.5*weightcm*weightang; instead, to deal with weights


#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <mpi.h>

#include "fmd/SlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Hamiltonian.h"
#include "fmd/Observables.h"

#include "numerics/wignerd.h"
#include "numerics/clebsch.h"

#include "misc/physics.h"
#include "misc/utils.h"

#include "Communication.h"
#include "Projectionmpi.h"


#define SQR(x) (x)*(x)


void calcprojectedMBMEmpi(const Projection* P, const ManyBodyOperator* Op,
			  const SlaterDet* Q, const SlaterDet* Qp,
			  Symmetry S, Symmetry Sp,
			  void* mbme)
{
  // only implemented for some operators

  if (!strncmp(Op->name, "Ovlap", 5)) {
    int task = TASKPROJECTOVLAPOD;
    BroadcastTask(&task);
  } else if (!strncmp(Op->name, "Observables-", 12)) {
    int task = TASKPROJECTOBSERVABLESOD;
    BroadcastTask(&task);
  } else if (!strncmp(Op->name, "EMonopoleFormfactor-", 20)) {
    int task = TASKPROJECTMONOPOLEFORMFACTOROD;
    BroadcastTask(&task);
  } else if (!strncmp(Op->name, "DiagonalDensityMatrixHO-", 24)) {
    int task = TASKPROJECTEDDIAGONALDENSITYMATRIXHOOD;
    BroadcastTask(&task);
  } else {
    fprintf(stderr, "No parallelized implementation for this Operator yet !\n");
    exit(127);
  }

  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  complex double (**val)[(rank+1)*size] = mbme;

  SlaterDetAux  X;

  int jmax = P->jmax;
  int odd = P->odd;

  // norms of SlaterDets
  //initSlaterDetAux(Q, &X);
  // calcSlaterDetAuxod(Q, Q, &X);
  // double norm = sqrt(creal(X.ovlap));
  
  //initSlaterDetAux(Qp, &X);
  //calcSlaterDetAuxod(Qp, Qp, &X);
  // double normp = sqrt(creal(X.ovlap));  

  // set up cm integration
  cmintegrationpara cmpara;
  initcmintegration(P, Q, Qp, &cmpara);
  int ncm = cmpara.n;

  // set up angular momentum integration
  angintegrationpara angpara;
  initangintegration(P, Q, Qp, S, Sp, &angpara);
  int nang = angpara.n;

  int l, r;
  int p, j, m, k;

  // set matrix elements to zero
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      for (k=-j; k<=j; k=k+2)
	for (m=-j; m<=j; m=m+2)
	  for (l=0; l<dim; l++)
	    for (r=0; r<=rank; r++)
	      val[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = 0.0;
    }	

  BroadcastSlaterDet(Q);
  BroadcastSlaterDet(Qp);

  // loop over all the orientations
  complex double w;
  int pi;

  int todo = nang*ncm;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double pos[mpisize][3];
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[6];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;

  double weightang, weightcm;
  complex double sval[2][(rank+1)*size];

  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getcmintegrationpoint(next/nang, &cmpara, pos[processor], &weightcm);
      getangintegrationpoint(next%nang, &angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weightang); 

      //  weight[processor] = 1.0/(2*norm*normp)*weightcm*weightang;
        weight[processor] = 0.5*weightcm*weightang;
      
      for (i=0; i<3; i++) projpar[i] = angle[processor][i];
      for (i=0; i<3; i++) projpar[i+3] = pos[processor][i];

      MPI_Send(projpar, 6, MPI_DOUBLE, processor, TAGPROJECT6, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;
	
    } else {					// collect results from slaves

      MPI_Recv(sval, 2*(rank+1)*size, MPI_DOUBLE_COMPLEX, 
	       MPI_ANY_SOURCE, TAGMEOD, MPI_COMM_WORLD, &status);
      processor = status.MPI_SOURCE;

      for (pi=0; pi<=1; pi++)
	for (p=0; p<=1; p++)
	  for (j=odd; j<jmax; j=j+2)
	    for (k=-j; k<=j; k=k+2)
	      for (m=-j; m<=j; m=m+2) {
		if ((Op-rank !=0 || SymmetryAllowed(S, p, j, m)) &&
		    SymmetryAllowed(Sp, p, j, k)) {
		  w = weight[processor] * (p && pi%2 ? -1 : 1)*
		    (j+1)/(8*SQR(M_PI))*Djmkstar(j,m,k,angle[processor][0],angle[processor][1],angle[processor][2]);
		  for (l=0; l<dim; l++)
		    for (r=0; r<=rank; r++)
		      val[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] += w*sval[pi][r+l*(rank+1)];
		}
	      }	

      working[processor] = 0;
      running--;
      done++;
    }

  }

  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 6, MPI_DOUBLE, processor, TAGPROJECT6, MPI_COMM_WORLD);
  }

}	


void calcprojectedMBMEsmpi(const Projection* P, const ManyBodyOperators* Ops,
			   const SlaterDet* Q, const SlaterDet* Qp,
			   Symmetry S, Symmetry Sp,
			   void* mbme)
{
  // only implemented for some operators

  if (!strncmp(Ops->name, "OneNucleonOvlaps-", 17)) {
    int task = TASKPROJECTONENUCLEONOVLAPSOD;
    BroadcastTask(&task);
  } else if (!strncmp(Ops->name, "TwoNucleonOvlapsT-", 18)) {
    int task = TASKPROJECTTWONUCLEONOVLAPSTOD;
    BroadcastTask(&task);
  } else if (!strncmp(Ops->name, "TwoNucleonOvlapsY-", 18)) {
    int task = TASKPROJECTTWONUCLEONOVLAPSYOD;
    BroadcastTask(&task);
  } else {
    fprintf(stderr, "No parallelized implementation for this Operator yet !\n");
    exit(127);
  }

  int size=Ops->size;
  int dim=Ops->dim;
  complex double ***val = mbme;

  SlaterDetAux  X;

  int jmax = P->jmax;
  int odd = P->odd;

  // norms of SlaterDets
  // initSlaterDetAux(Q, &X);
  // calcSlaterDetAuxod(Q, Q, &X);
  // double norm = sqrt(creal(X.ovlap));
  
  // initSlaterDetAux(Qp, &X);
  //  calcSlaterDetAuxod(Qp, Qp, &X);
  //  double normp = sqrt(creal(X.ovlap));  

  // set up cm integration
  cmintegrationpara cmpara;
  initcmintegration(P, Q, Qp, &cmpara);
  int ncm = cmpara.n;

  // set up angular momentum integration
  angintegrationpara angpara;
  initangintegration(P, Q, Qp, S, Sp, &angpara);
  int nang = angpara.n;

  int l, r;
  int o, p, j, m, k;

  // helpful for indexing matrix elements

  int ranko[Ops->n];
  for (o=0; o<Ops->n; o++)
    ranko[o] = Ops->Op[o].rank;

  int sizeo[Ops->n];
  for (o=0; o<Ops->n; o++)
    sizeo[o] = size*(ranko[o]+1);

  int io[Ops->n]; io[0] = 0;
  for (o=0; o<Ops->n-1; o++)
    io[o+1] = io[o]+sizeo[o];

  int no=0;
  for (o=0; o<Ops->n; o++)
    no += sizeo[o];


  // set matrix elements to zero
  for (o=0; o<Ops->n; o++)
    for (p=0; p<=1; p++)
      for (j=odd; j<jmax; j=j+2) {
	for (k=-j; k<=j; k=k+2)
	  for (m=-j; m<=j; m=m+2)
	    for (l=0; l<dim; l++)
	      for (r=0; r<=ranko[o]; r++)
		val[o][idxpij(jmax,p,j)][r+l*(ranko[o]+1)+idxjmk(j,m,k)*sizeo[o]] = 0.0;
    }	

  BroadcastSlaterDet(Q);
  BroadcastSlaterDet(Qp);

  // loop over all the orientations
  complex double w;
  int pi;

  int todo = nang*ncm;
  int running=0;
  int done=0;
  int processor;
  MPI_Status status;

  // processor 0 is master
  double pos[mpisize][3];
  double angle[mpisize][3];
  double weight[mpisize];
  double projpar[6];
  
  // which slaves are working, at beginning none
  int i;
  int working[mpisize];
  for (i=1; i<mpisize; i++)
    working[i] = 0;

  double weightang, weightcm;
  complex double sval[2][no];

  int next=0;
  while (done<todo) {

    if (next<todo && running < mpisize-1) {	// supply slaves with work

      // find empty slave
      for (i=1; i<mpisize; i++)
	if (!working[i]) {
	  processor = i;
	  break;
	}

      getcmintegrationpoint(next/nang, &cmpara, pos[processor], &weightcm);
      getangintegrationpoint(next%nang, &angpara, &angle[processor][0], &angle[processor][1], &angle[processor][2], &weightang); 

    //  weight[processor] = 1.0/(2*norm*normp)*weightcm*weightang;
         weight[processor] = 0.5*weightcm*weightang;
      
      for (i=0; i<3; i++) projpar[i] = angle[processor][i];
      for (i=0; i<3; i++) projpar[i+3] = pos[processor][i];

      MPI_Send(projpar, 6, MPI_DOUBLE, processor, TAGPROJECT6, MPI_COMM_WORLD);

      working[processor] = 1;
      running++;
      next++;
	
    } else {					// collect results from slaves

      MPI_Recv(sval, 2*no, MPI_DOUBLE_COMPLEX, 
	       MPI_ANY_SOURCE, TAGMEOD, MPI_COMM_WORLD, &status);
      processor = status.MPI_SOURCE;

      for (pi=0; pi<=1; pi++)
	for (o=0; o<Ops->n; o++)
	  for (p=0; p<=1; p++)
	    for (j=odd; j<jmax; j=j+2)
	      for (k=-j; k<=j; k=k+2)
		for (m=-j; m<=j; m=m+2) {
		  if ((ranko[o] !=0 || SymmetryAllowed(S, p, j, m)) &&
		      SymmetryAllowed(Sp, p, j, k)) {
		    w = weight[processor] * (p && pi%2 ? -1 : 1)*
		      (j+1)/(8*SQR(M_PI))*Djmkstar(j,m,k,angle[processor][0],angle[processor][1],angle[processor][2]);
		    for (l=0; l<dim; l++)
		      for (r=0; r<=ranko[o]; r++)
			val[o][idxpij(jmax,p,j)][r+l*(ranko[o]+1)+idxjmk(j,m,k)*sizeo[o]] += w*sval[pi][r+l*(ranko[o]+1)+io[o]];
		  }
		}	

      working[processor] = 0;
      running--;
      done++;
    }

  }

  // tell slaves we have finished
  projpar[0] = -1.0;
  for (processor=1; processor<mpisize; processor++) {
    MPI_Send(projpar, 6, MPI_DOUBLE, processor, TAGPROJECT6, MPI_COMM_WORLD);
  }

}	

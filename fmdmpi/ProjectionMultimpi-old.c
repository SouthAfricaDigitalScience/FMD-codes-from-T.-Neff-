/**

  \file ProjectionMultimpi.c

  angular momentum projection of matrix elements


  (c) 2003, 2004, 2005 Thomas Neff

*/


#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <mpi.h>

#include "fmd/SlaterDet.h"
#include "fmd/MultiSlaterDet.h"
#include "fmd/Observables.h"

#include "numerics/wignerd.h"
#include "numerics/clebsch.h"

#include "misc/physics.h"
#include "misc/utils.h"

#include "Communication.h"
#include "Projectionmpi.h"


#define SQR(x) (x)*(x)
#define MAX(x,y) ((x)<(y) ? (y) : (x))


// parallelization only in the inner-loop for angular momentum integration
// this is not very efficient, especially if MultiSlaterDets have axial symmetry

// todo: parallelize including outer loops 

static inline double dmin(double a, double b)
{
  return (a < b ? a : b);
}

static inline double dmax(double a, double b)
{
  return (a > b ? a : b);
}

static double dminarray(double a[], int n)
{
  double m = a[0];
  int i;
  for (i=1; i<n; i++)
    m = dmin(m, a[i]);

  return m;
}

static double dmaxarray(double a[], int n)
{
  double m = a[0];
  int i;
  for (i=1; i<n; i++)
    m = dmax(m, a[i]);

  return m;
}



void calcprojectedMultiMBMEmpi(const Projection* P, const ManyBodyOperator* Op,
			       const MultiSlaterDet* MBA, 
			       const MultiSlaterDet* MBB,
			       void** mbme)
{
  // only implemented for Observablesod by now

  int task;
  if (!strncmp(Op->name, "Observables-", 12)) {
    task = TASKPROJECTOBSERVABLESOD;
  } else {
    fprintf(stderr, "No parallelized implementation for this Operator yet !\n");
    exit(127);
  }

  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  complex double (***val)[(rank+1)*size] = mbme;

  int jmax = P->jmax;
  int odd = P->odd;

  // loop over the individual SlaterDets in the MultiSlaterDets
  int nA = MBA->n; int nB = MBB->n;
  int iA, iB;

  int NA = MBA->N; int NB = MBB->N;
  int IA, IB;

  SlaterDet Q, Qp, Qpp;
  allocateSlaterDet(&Q, MBA->A);
  allocateSlaterDet(&Qp, MBB->A);
  allocateSlaterDet(&Qpp, MBB->A);

  SlaterDetAux X;
  allocateSlaterDetAux(&X, MAX(MBA->A, MBB->A));

  Symmetry S, Sp;
  int axialsym;

  int l, r;
  int p, j, m, k;


  // set matrix elements to zero
  for (IB=0; IB<NB; IB++)
    for (IA=0; IA<NA; IA++)
      for (p=0; p<=1; p++)
	for (j=odd; j<jmax; j=j+2) {
	  for (k=-j; k<=j; k=k+2)
	    for (m=-j; m<=j; m=m+2)
	      for (l=0; l<dim; l++)
		for (r=0; r<=rank; r++)
		  val[IA+IB*NA][idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = 0.0;
	}			


  // we need a common denominaotor for the NA*NB matrix elements
  // all states have axial symmetry ? use axial as indicator
  // also spherical will be tretated like axial

  S=0;
  axialsym=1;
  for (IA=0; IA<NA; IA++)
    axialsym &= hasAxialSymmetry(MBA->symmetry(MBA, IA));
  if (axialsym)
    setSymmetry(&S, axial);

  Sp=0;
  axialsym=1;
  for (IB=0; IB<NB; IB++)
    axialsym &= hasAxialSymmetry(MBB->symmetry(MBB, IB));
  if (axialsym)
    setSymmetry(&Sp, axial);

  fprintf(stderr, "Symmetry - MBA: %s, MBB: %s\n", SymmetrytoStr(S), SymmetrytoStr(Sp)); 

  double kappaA[nA], kappaB[nB];
  double acmA[nA], acmB[nB];
  
  for (iA=0; iA<nA; iA++) {
    MBA->get(MBA, iA, &Q);
    kappaA[iA] = _estimateangkappa(&Q);
    acmA[iA] = _estimateacm(&Q);
  }

  for (iB=0; iB<nB; iB++) {
    MBB->get(MBB, iB, &Qp);
    kappaB[iB] = _estimateangkappa(&Qp);
    acmB[iB] = _estimateacm(&Qp);
  }

  fprintf(stderr, "kappacrit: %6.2f\n", _getangkappacrit());
  fprintf(stderr, "kappa - MBA: [%6.2f - %6.2f], MBB: [%6.2f - %6.2f]\n",
          dminarray(kappaA, nA), dmaxarray(kappaA, nA),
          dminarray(kappaB, nB), dmaxarray(kappaB, nB));

  int c=0; 
  int cmax=nA*nB; 

  for (iB=0; iB<nB; iB++) {
    MBB->get(MBB, iB, &Qp);

    for (iA=0; iA<nA; iA++) {
      MBA->get(MBA, iA, &Q);

      // progress indicator
      c++;
      if (c%100==0) fprintf(stderr, "%d%%", (100*c)/cmax);
      if (c%10==0) fprintf(stderr, ".");

      // norms of SlaterDets
      calcSlaterDetAuxod(&Q, &Q, &X);
      double norm = sqrt(creal(X.ovlap));

      calcSlaterDetAuxod(&Qp, &Qp, &X);
      double normp = sqrt(creal(X.ovlap));  

      // set up cm integration
      cmintegrationpara cmpara;
      double cmalpha = 0.5/(acmA[iA]+acmB[iB]);
      _initcmintegration(P, cmalpha, &cmpara);
      int ncm = cmpara.n;

      // set up ang integration
      angintegrationpara angpara;
      double angkappa = dmin(kappaA[iA], kappaB[iB]);
      _initangintegration(P, angkappa, S, Sp, &angpara);
      int nang = angpara.n;

      BroadcastTask(&task);
      BroadcastSlaterDet(&Q);
      BroadcastSlaterDet(&Qp);

      // loop over all the orientations
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

	 weight[processor] = 1.0/(2*norm*normp)*weightcm*weightang;
      
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


	  complex double w, wA, wB;
	  Symmetry SA, SB;
	  for (IB=0; IB<NB; IB++) {
	    SB = MBB->symmetry(MBB, IB);
	    wB = MBB->weight(MBB, IB, iB);
	    for (IA=0; IA<NA; IA++) {
	      SA = MBA->symmetry(MBA, IA);
	      wA = MBA->weight(MBA, IA, iA);
	      for (pi=0; pi<=1; pi++)
		for (p=0; p<=1; p++)
		  for (j=odd; j<jmax; j=j+2)
		    for (k=-j; k<=j; k=k+2)
		      for (m=-j; m<=j; m=m+2) {
			if ((Op->rank != 0 || SymmetryAllowed(SA, p, j, m)) &&
			    SymmetryAllowed(SB, p, j, k)) {
			  w = conj(wA)*wB*	 
			    weight[processor] * (p && pi%2 ? -1 : 1)*
			    (j+1)/(8*SQR(M_PI))*Djmkstar(j,m,k,angle[processor][0],angle[processor][1],angle[processor][2]);
			  for (l=0; l<dim; l++)
			    for (r=0; r<=rank; r++)
			      val[IA+IB*NA][idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] += w*sval[pi][r+l*(rank+1)];
		    }	
		  }
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

      freeAngintegration(&angpara);
      freecmintegration(&cmpara);

    }
  }

}	

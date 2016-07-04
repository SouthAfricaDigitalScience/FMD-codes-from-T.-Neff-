/**

  \file calcoccupationnumbershoprojmes.c

  calculate harmonic oscillator occupation number matrix elements


  (c) 2006-2008 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/HOBasis.h"
#include "fmd/ProjectedDensityMatrixHO.h"

#include "numerics/zcw.h"
#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/Projectionmpi.h"
#include "fmdmpi/DiagonalDensityHOSlave.h"
#endif


#define MAXSTATES 100


void cleanup(int ret)
{
#ifdef MPI
  int task=TASKFIN;
  BroadcastTask(&task);

  MPI_Finalize();
#endif

  exit(ret);
}


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    DiagonalDensityHOSlave();

    MPI_Finalize();
  } else {
#endif

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPAR NUCSFILE"
	    "\n   -n NMAX	oscillator cut off"
	    "\n   -o OMEGA	Oscillator constant [MeV]\n",
	    argv[0]);
    exit(-1);
  }

  int nmax = HONMAX;
  double omega = 0.0;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "n:o:")) != -1)
    switch (c) {
    case 'n':
      nmax = atoi(optarg);
      break;
    case 'o':
      omega = atof(optarg)/hbc;
      break;
    }

  char* projpar = argv[optind];
  char* nucsfile = argv[optind+1];

  char* mbfile[MAXSTATES];
  int n;

  if (readstringsfromfile(nucsfile, &n, mbfile))
    cleanup(-1);

  SlaterDet Q[n];
  Symmetry S[n];

  int i;
  for (i=0; i<n; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      cleanup(-1);
  }


  // odd numer of nucleons ?
  int odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  if (omega == 0.0) {
    fprintf(stderr, "You have to provide a oscillator parameter\n");
    exit(-1);
  }

  // density matrix

  DensityMatrixHOPar DMpar = {
    nmax : nmax,
    omega : omega,
    dim : dimHOBasis(nmax)
  };

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastParameters(&DMpar, sizeof(DensityMatrixHOPar));
  BroadcastA(&Q[0].A);
#endif

#ifndef MPI
  initHOBasis(nmax);
#endif
  initOpDiagonalDensityMatrixHO(&DMpar);

  void* dmme = initprojectedMBME(&P, &OpDiagonalDensityMatrixHO);
  
  // read or calculate matrix elements

  int a,b;
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {

      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P,
				    &OpDiagonalDensityMatrixHO, S[a], S[b],
				    dmme)) {
#ifdef MPI
	calcprojectedMBMEmpi(&P, &OpDiagonalDensityMatrixHO, &Q[a], &Q[b],
			     S[a], S[b], dmme);
#else
	calcprojectedMBME(&P, &OpDiagonalDensityMatrixHO, &Q[a], &Q[b],
			  S[a], S[b], dmme);
#endif
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpDiagonalDensityMatrixHO, S[a], S[b], 
				 dmme);
      }
    }	


  cleanup(0);

#ifdef MPI
  }
#endif

}

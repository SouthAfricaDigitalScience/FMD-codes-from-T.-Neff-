/**

  \file recalcenergymultiproj.c

  calculate energies with new interaction for given many-body states


  (c) 2010 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/Projection.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"


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

  // fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    ProjectionSlave();

    MPI_Finalize();
  } else {
#endif

  /* enough arguments ? */

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [OPTIONS] INTERACTION mcstate"
	    "\n   -A             show all eigenstates\n", argv[0]);
    exit(-1);
  }

  int all=0;
  int hermit=0;

  char c;

  /* manage command-line options */

  while ((c = getopt(argc, argv, "A")) != -1)
    switch (c) {
    case 'A':
      all=1;
      break;
    }

  char* interactionfile = argv[optind];
  char* mcstatefile = argv[optind+1];
  char** mbfile;

  // open multiconfigfile
  Projection P;
  SlaterDet* Q;
  Symmetry* S;
  Eigenstates E;
  int n;

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    cleanup(-1);
  Int.cm = 1;

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastInteraction(&Int);
  BroadcastA(&Q[0].A);
#endif

  initOpObservables(&Int);

  readMulticonfigfile(mcstatefile, &mbfile, &P, &Q, &S, &E, &n);

  int a,b; 

  Observablesod** obsme = malloc(n*n*sizeof(Observablesod**));
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++)
      obsme[a+b*n] = initprojectedMBME(&P, &OpObservables);

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++)
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, &OpObservables, 
				    S[a], S[b], obsme[a+b*n])) {
#ifdef MPI
	calcprojectedMBMEmpi(&P, &OpObservables, 
			     &Q[a], &Q[b], S[a], S[b], obsme[a+b*n]);
#else
	calcprojectedMBME(&P, &OpObservables, 
			  &Q[a], &Q[b], S[a], S[b], obsme[a+b*n]);
#endif
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], 
				 &P, &OpObservables, S[a], S[b], obsme[a+b*n]);
      }

  // scale matrix elements
  if (Int.mescaling) {
    fprintf(stderr, "... scaling matrix elements\n");
    for (b=0; b<n; b++)
      for (a=0; a<n; a++)
        scaleprojectedObservablesMBME(&P, &Int, obsme[a+b*n]);
  }

  if (hermit)
    hermitizeprojectedMBME(&P, &OpObservables, obsme, n);

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  Observablesod** obs = initprojectedVector(&P, &OpObservables, n);
  fprintf(stderr, "... calcexpectprojectedMBME\n");
  calcexpectprojectedMBME(&P, &OpObservables, obsme, S, &E, obs);

  // output

  showprojectedObservables(stdout, &P, &Int, &Q[0], obs, &E, NULL, "");


  cleanup(0);

#ifdef MPI
  }
#endif

}



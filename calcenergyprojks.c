/**

  \file calcenergyprojks.c

  project SlaterDet on good angular momentum


  (c) 2003,2004,2005 Thomas Neff

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
#include "fmd/Symmetry.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/Projectionmpi.h"
#include "fmdmpi/ProjectionSlave.h"
#endif


void cleanup(int ret)
{
  // stop all mpi clients

#ifdef MPI
  int task=TASKFIN;
  BroadcastTask(&task);

  MPI_Finalize();
#endif
 
  // terminate program with return code
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
    ProjectionSlave();

    MPI_Finalize();
  } else {
#endif

  // enough arguments ?

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA INTERACTION [SYMMETRY:]MBSTATE"
	    "\n   -h                hermitize matrix elements"
	    "\n   -K K              max K"
	    "\n   -A                show really all eigenstates\n",
	    filepart(argv[0]));
    cleanup(-1);
  }

  int all=0;
  int hermit=0;
  int Kmax=999, K;
  int odd;
  double threshkmix=0.00001;
  double minnorm=0.005;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "hK:At:")) != -1)
    switch (c) {
    case 'h':
      hermit=1;
      break;
    case 'K':
      Kmax = atoi(optarg);
      break;
    case 'A':
      all=1;
      break;
    case 't':
      threshkmix = atof(optarg);
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* mbfile = argv[optind+2];

  // get Symmetry
  Symmetry S;
  extractSymmetryfromString(&mbfile, &S);

  // only SlaterDets can be used as many-body states up to now
  SlaterDet Q;
  if (readSlaterDetfromFile(&Q, mbfile))
    cleanup(-1);

  // get rid of path in filename
  char* mbname = filepart(mbfile);

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    cleanup(-1);
  Int.cm = 1;

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);
  
  BroadcastInteraction(&Int);
  BroadcastA(&Q.A);
#endif

  // odd numer of nucleons ?
  odd = Q.A % 2;

  // defaults for Kmax
  if (Kmax == 999)
    Kmax = odd ? 7 : 8; 

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  initOpObservables(&Int);

  // initialize space for matrix elements
  Observablesod **obsme = initprojectedMBME(&P, &OpObservables);

  // try to read ME's from file, if not possible calculate and save

  if (readprojectedMBMEfromFile(mbfile, mbfile, &P, &OpObservables, S, S, obsme)) {
#ifdef MPI
    calcprojectedMBMEmpi(&P, &OpObservables, &Q, &Q, S, S, obsme);
#else
    calcprojectedMBME(&P, &OpObservables, &Q, &Q, S, S, obsme);
#endif
    writeprojectedMBMEtoFile(mbfile, mbfile, &P, &OpObservables, S, S, obsme);
  }

  // scale matrix elements
  if (Int.mescaling) {
    fprintf(stderr, "... scaling matrix elements\n");
    scaleprojectedObservablesMBME(&P, &Int, obsme);
  }

  // hermitize MEs
  if (hermit)
    hermitizeprojectedMBME(&P, &OpObservables, &obsme, 1);

  // write output
  FILE *outfp; 
  char outfilename[255];

  snprintf(outfilename, 255, "%s%s.%s-K", mbfile, 
	   S==0 ? "" : strjoin(":", SymmetrytoStr(S)),  
	   ProjectiontoStr(&P));

  if (!(outfp = fopen(outfilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfilename);
    cleanup(-1);
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);
  fprintf(outfp, "# used Symmetry - %s\n", SymmetrytoStr(S));

  // first full K-mixing

  fprintf(outfp, "\n# full K-mixing\n");

  // solve the eigenvalue problem
  Eigenstates E;
  initEigenstates(&P, &E, 1);
  calcEigenstates(&P, &Int, &obsme, &E, threshkmix);

  // calculate expectation values
  Observablesod** obs = initprojectedVector(&P, &OpObservables, 1);
  calcexpectprojectedMBME(&P, &OpObservables, &obsme, &S, &E, obs);
    
  // sort Eigenstates
  sortEigenstates(&P, &Int, obs, &E, minnorm, all);

  // output Observables
  showprojectedObservables(outfp, &P, &Int, &Q, obs, &E, NULL, "");

  // now for each K-projection

  for (K=-Kmax; K<=Kmax; K=K+2) {

    char pre[12];

    if (P.odd) sprintf(pre, "[K = %+d/2]", K);
    else sprintf(pre, "[K = %+d]", K/2);

    calcEigenstatesK(&P, &Int, &obsme, &E, K, threshkmix);

    // calculate expectation values
    calcexpectprojectedMBME(&P, &OpObservables, &obsme, &S, &E, obs);
    
    // sort Eigenstates
    sortEigenstates(&P, &Int, obs, &E, minnorm, all);

    // output Observables
    showprojectedObservables(outfp, &P, &Int, &Q, obs, &E, NULL, pre);

  }

  fprintf(outfp, "\n# end of file\n");

  fclose(outfp);

  cleanup(0);

#ifdef MPI
  }
#endif
}



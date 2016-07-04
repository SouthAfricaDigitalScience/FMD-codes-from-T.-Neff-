/**

  \file calcenergyproj.c

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

  // fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    ProjectionSlave();

    MPI_Finalize();
  } else {
#endif

  // enough arguments ?

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA INTERACTION [SYMMETRY:]MBSTATE"
	    "\n   -h                hermitize matrix elements"
	    "\n   -K K              use only K projection"
            "\n   -n NORM           set minimal norm for K-mixing eigenstates"
	    "\n   -t THRESH         set threshold for K-mixing SVD"
	    "\n   -A                show really all eigenstates"
            "\n   -l                write energy level file"
	    "\n   -s                write Eigenstates into file\n",
	    filepart(argv[0]));
    cleanup(-1);
  }

  int all=0;
  int levels=0;
  int hermit=0;
  int Ksel=0, K;
  int odd;
  int savestates=0;
  double threshkmix=0.01;
  double minnormkmix=0.001;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "hK:Alst:n:")) != -1)
    switch (c) {
    case 'h':
      hermit=1;
      break;
    case 'K':
      Ksel=1;
      K = atoi(optarg);
      break;
    case 'A':
      all=1;
      break;
    case 'l':
      levels=1;
      break;
    case 's':
      savestates=1;
      break;
    case 'n':
      minnormkmix = atof(optarg);
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

  // odd number of nucleons ?
  odd = Q.A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  _setangkappacrit(25.0);

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

  // solve the eigenvalue problem
  Eigenstates E;
  initEigenstates(&P, &E, 1);
  if (Ksel) 
    calcEigenstatesK(&P, &Int, &obsme, &E, K, threshkmix);
  else
    calcEigenstates(&P, &Int, &obsme, &E, threshkmix);

  // calculate expectation values
  Observablesod** obs = initprojectedVector(&P, &OpObservables, 1);
  calcexpectprojectedMBME(&P, &OpObservables, &obsme, &S, &E, obs);
    
  // sort Eigenstates
  sortEigenstates(&P, &Int, obs, &E, minnormkmix, all);

  // write output
  FILE *outfp; 
  char outfilename[255];

  snprintf(outfilename, 255, "%s%s.%s", mbfile, 
	   S==0 ? "" : strjoin(":", SymmetrytoStr(S)),  
	   ProjectiontoStr(&P));

  if (!(outfp = fopen(outfilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfilename);
    cleanup(-1);
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);
  fprintf(outfp, "# used Symmetry - %s\n", SymmetrytoStr(S));
  fprintf(outfp, "# K-mixing,        threshold: %g, minnorm: %g\n", 
          threshkmix, minnormkmix);

  // output Observables

  showprojectedObservables(outfp, &P, &Int, &Q, obs, &E, NULL, "");

  fprintf(outfp, "\n# end of file\n");

  fclose(outfp);

  if (savestates)
    writeMulticonfigfile(outfilename, &P, &S, &mbname, &E, 1);

  // levels output ?
  if (levels) {

    snprintf(outfilename, 255, "%s.%s.levels", mbfile, ProjectiontoStr(&P));
    if (!(outfp = fopen(outfilename, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", outfilename);
      cleanup(-1);
    }

    fprintinfo(outfp);
    fprintProjectinfo(outfp, &P);
    fprintf(outfp, "# K-mixing,        threshold: %g, minnorm: %g\n", 
            threshkmix, minnormkmix);

    // output levels

    showSpectrum(outfp, &P, obs, &E);
  
    fclose(outfp);
  }

  cleanup(0);

#ifdef MPI
  }
#endif
}



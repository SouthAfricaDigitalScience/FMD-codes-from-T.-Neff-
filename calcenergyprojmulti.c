/**

  \file calcenergyprojmulti.c

  project Many-body States on good angular momentum


  (c) 2003-2009 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/MultiSlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/ProjectionMultimpi.h"
#include "fmdmpi/ProjectionSlave.h"
#endif


#define SQR(x) ((x)*(x))


void cleanup(int ret)
{
  // indicate mpi slaves to finish

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
    // do the work
    ProjectionSlave();

    MPI_Finalize();
  } else {
#endif

  // enough arguments ?

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA INTERACTION [IDX:]MBSTATE"
	    "\n   -h                hermitize matrix elements"
	    "\n   -K K              use only K projection"
	    "\n   -A                show really all eigenstates"
	    "\n   -s                write Eigenstates into file\n",
	    filepart(argv[0]));
    cleanup(-1);
  }

  int all=0;
  int hermit=0;
  int Ksel=0, K;
  int odd;
  int savestates=0;
  double threshkmix=0.00001;
  double minnorm=0.005;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "hK:Ast:")) != -1)
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
    case 's':
      savestates=1;
      break;
    case 't':
      threshkmix = atof(optarg);
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* mbfile = argv[optind+2];

  // many-body states up to now
  MultiSlaterDet Q;
  Indices In;

  extractIndicesfromString(&mbfile, &In);

  if (readMultiSlaterDetfromFile(&Q, &In, mbfile))
    cleanup(-1);

  if (In.n > 1) {
    fprintf(stderr, "specify a single index\n");
    cleanup(-1);
  }

  int idx = In.idx[0];

  int N = Q.N;

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

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  _setangkappacrit(25.0);

  initOpObservables(&Int);

  // initialize space for matrix elements
  Observablesod ***obsme = initprojectedMultiMBME(&P, &OpObservables,
						  &Q, &Q);

  // try to read ME's from file, if not possible calculate and save

  if (readprojectedMultiMBMEfromFile(mbfile, mbfile, &Q, &Q, 
				     &P, &OpObservables, obsme)) {
#ifdef MPI
    calcprojectedMultiMBMEmpi(&P, &OpObservables, &Q, &Q, obsme);
#else
    calcprojectedMultiMBME(&P, &OpObservables, &Q, &Q, obsme);
#endif
    writeprojectedMultiMBMEtoFile(mbfile, mbfile, &Q, &Q, 
				  &P, &OpObservables, obsme);
  }

  // hermitize MEs
  if (hermit)
    hermitizeprojectedMBME(&P, &OpObservables, &obsme[idx+idx*N], 1);

  fprintf(stderr, "... solving Eigenvalue problems\n");

  // solve the eigenvalue problem
  Eigenstates E;
  initEigenstates(&P, &E, 1);

  if (Ksel) 
    calcEigenstatesK(&P, &Int, &obsme[idx+idx*N], &E, K, threshkmix);
  else
    calcEigenstates(&P, &Int, &obsme[idx+idx*N], &E, threshkmix);

  fprintf(stderr, "... calculating expectation values\n");
  // calculate expectation values
  Symmetry S = Q.symmetry(&Q, idx);
  Observablesod** obs = initprojectedVector(&P, &OpObservables, 1);
  calcexpectprojectedMBME(&P, &OpObservables, &obsme[idx+idx*N], &S, &E, obs);
    
  // sort Eigenstates
  sortEigenstates(&P, &Int, obs, &E, minnorm, all); 
  // write output
  FILE *outfp; 
  char outfilename[255];

  snprintf(outfilename, 255, "%s:%d.%s", mbfile, idx,  
	   ProjectiontoStr(&P));

  if (!(outfp = fopen(outfilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfilename);
    cleanup(-1);
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  // output Observables

  // only needed for output
  SlaterDet Q0;
  allocateSlaterDet(&Q0, Q.A);
  Q.get(&Q, idx, &Q0);
  
  showprojectedObservables(outfp, &P, &Int, &Q0, obs, &E, NULL, "");

  fprintf(outfp, "\n# end of file\n");

  fclose(outfp);

  S = Q.symmetry(&Q, idx);

  //  if (savestates)
  //  writeMultiEigenstatestoFile(outfilename, &P, mbname, &In, &E);

  if (savestates)
    writeMultiMulticonfigfile(outfilename, &P, &mbfile, &In, 1, &E);

  cleanup(0);

#ifdef MPI
  }
#endif
}



/**

  \file calcenergymultiproj.c

  multiconfiguration calculations with projected Slater determinants


  (c) 2003, 2004, 2005 Thomas Neff

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

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPAR INTERACTION NUCSFILE"
	    "\n   -h                hermitize matrix elements"
	    "\n   -A                show really all eigenstates"
	    "\n   -s                write Eigenstates into file"
            "\n   -l                write Energy Level file"
            "\n   -n NORM           set minimal norm for K-mixing eigenstates"
	    "\n   -t THRESH         set threshold for K-mixing SVD"
            "\n   -N NORM           set minimal norm for Multiconfig eigenstates"
	    "\n   -T THRESH         set threshold for Multiconfig SVD\n", 
	    argv[0]);
    return -1;
  }

  int all=0;
  int levels=0;
  int hermit=0;
  int odd;
  double threshkmix=0.01;
  double minnormkmix=0.001;
  double threshmulti=0.0000001;
  double minnormmulti=0.001;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "hAln:t:N:T:")) != -1)
    switch (c) {
    case 'h':
      hermit=1;
      break;
    case 'A':
      all=1;
      break;
    case 'l':
      levels=1;
      break;
    case 'n':
      minnormkmix=atof(optarg);
      break;
    case 't':
      threshkmix = atof(optarg);
      break;
    case 'N':
      minnormmulti=atof(optarg);
      break;
    case 'T':
      threshmulti = atof(optarg);
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* nucsfile = argv[optind+2];

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

  // odd numer of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);
    
  initOpObservables(&Int);

  int a,b; 

  // initialize space for matrix elements
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


  // read or calculate the Eigenstates
      
  int allp=0;
  Eigenstates Ep[n];
  Observablesod** obsp = initprojectedVector(&P, &OpObservables, 1);
  
  for (i=0; i<n; i++) {
    if (readEigenstatesfromFile(mbfile[i], &P, &Ep[i], 1)) {
      fprintf(stderr, "... calculating Eigenstates for %s\n", mbfile[i]);
      initEigenstates(&P, &Ep[i], 1);
      calcEigenstates(&P, &Int, &obsme[i+i*n], &Ep[i], threshkmix);
      calcexpectprojectedMBME(&P, &OpObservables, &obsme[i+i*n], &S[i], &Ep[i], obsp);
      sortEigenstates(&P, &Int, obsp, &Ep[i], minnormkmix, allp);
    }
  }

  // solve the n SlaterDet eigenvalue problem
  Eigenstates multiE;
  Amplitudes multiA;
  initEigenstates(&P, &multiE, n);
  initAmplitudes(&P, &multiA, n);
  fprintf(stderr, "... calcMultiEigenstates\n");
  calcMultiEigenstates(&P, &Int, obsme, &Ep, &multiE, &multiA, threshmulti);

  // calculate expectation values
  Observablesod** obs = initprojectedVector(&P, &OpObservables, n);
  fprintf(stderr, "... calcexpectprojectedMBME\n");
  calcexpectprojectedMBME(&P, &OpObservables, obsme, S, &multiE, obs);

  // sort Eigenstates
  fprintf(stderr, "... sortEigenstates\n");
  sortEigenstates(&P, &Int, obs, &multiE, minnormmulti, all);

  // output
  FILE *outfp;
  char outfilename[255];

  snprintf(outfilename, 255, "%s.multi-%s", nucsfile, ProjectiontoStr(&P));
  if (!(outfp = fopen(outfilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfilename);
    cleanup(-1);
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);
  fprintf(outfp, "# K-mixing,        threshold: %g, minnorm: %g\n", 
          threshkmix, minnormkmix);
  fprintf(outfp, "# Diagonalization, threshold: %g, minnorm: %g\n",
          threshmulti, minnormmulti);

  // output Observables

  showprojectedObservables(outfp, &P, &Int, &Q[0], obs, &multiE, &multiA, "");

  fclose(outfp);

  writeMulticonfigfile(outfilename, &P, S, mbfile, &multiE, n);

  // levels output ?
  if (levels) {

    snprintf(outfilename, 255, "%s.multi-%s.levels", nucsfile, ProjectiontoStr(&P));
    if (!(outfp = fopen(outfilename, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", outfilename);
      cleanup(-1);
    }

    fprintinfo(outfp);
    fprintProjectinfo(outfp, &P);
    fprintf(outfp, "# K-mixing,        threshold: %g, minnorm: %g\n", 
            threshkmix, minnormkmix);
    fprintf(outfp, "# Diagonalization, threshold: %g, minnorm: %g\n",
            threshmulti, minnormmulti);

    // output levels

    showSpectrum(outfp, &P, obs, &multiE);
  
    fclose(outfp);
  }

  cleanup(0);

#ifdef MPI
  }
#endif

}

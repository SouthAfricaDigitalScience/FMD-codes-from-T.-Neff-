/**

  \file calcenergymultiprojmulti.c

  multiconfiguration calculations with projected Slater determinants


  (c) 2003-2009 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/MultiSlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/Symmetry.h"
#include "fmd/Projection.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/ProjectionMultimpi.h"
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
  int odd;
  double threshkmix=0.01;
  double minnormkmix=0.001;
  double threshmulti=0.0000001;
  double minnormmulti=0.001;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "Aln:t:N:T:")) != -1)
    switch (c) {
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
      threshkmix=atof(optarg);
      break;
    case 'N':
      minnormmulti=atof(optarg);
      break;
    case 'T':
      threshmulti=atof(optarg);
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* nucsfile = argv[optind+2];

  char* mbfile[MAXSTATES];
  int n;

  if (readstringsfromfile(nucsfile, &n, mbfile))
    return -1;

  if (n>MAXSTATES) {
    fprintf(stderr, "MAXSTATES too small, aborting\n");
    return -1;
  }

  MultiSlaterDet Q[n]; 
  Indices In[n];
  
  int i;
  for (i=0; i<n; i++) {
    extractIndicesfromString(&mbfile[i], &In[i]);
    if (readMultiSlaterDetfromFile(&Q[i], &In[i], mbfile[i]))
      return -1;
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    return -1;
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
    
  // change kappacrit for switching from Gauss-Legendre to Gauss-Exponential in beta integration
  _setangkappacrit(25.0);

  initOpObservables(&Int);

  int a,b; 

  // initialize space for matrix elements
  Observablesod*** obsme[n*n];
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++)
      obsme[a+b*n] = initprojectedMultiMBME(&P, &OpObservables, &Q[a], &Q[b]);

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++)
      if (readprojectedMultiMBMEfromFile(mbfile[a], mbfile[b], &Q[a], &Q[b],
					 &P, &OpObservables, obsme[a+b*n])) {
#ifdef MPI
	calcprojectedMultiMBMEmpi(&P, &OpObservables, &Q[a], &Q[b], obsme[a+b*n]);
#else
	calcprojectedMultiMBME(&P, &OpObservables, &Q[a], &Q[b], obsme[a+b*n]);
#endif
	writeprojectedMultiMBMEtoFile(mbfile[a], mbfile[b], &Q[a], &Q[b], 
				      &P, &OpObservables, obsme[a+b*n]);
      }

  // scale matrix elements
  if (Int.mescaling) {
    fprintf(stderr, "scaling not yet implemented, aborting\n");
    return -1;
  }


  // dimension of many-body basis
  // beware: matrix elements are calculated in (possibly) larger basis

  int dim=0;
  for (i=0; i<n; i++)
    dim += In[i].n;


  // read or calculate the Eigenstates
      
  int allp=0;
  Eigenstates Ep[dim];
  Observablesod** obsp = initprojectedVector(&P, &OpObservables, 1);
  Symmetry Sp;
  
  int idx=-1;
  int ii;
  for (i=0; i<n; i++)
    for (ii=0; ii<In[i].n; ii++) {
      idx++;
      int N=In[i].N;
      int idxii=In[i].idx[ii];
      Sp = Q[i].symmetry(&Q[i], idxii);
      // if (readMultiEigenstatesfromFile(mbfile[i], idxii, &P, &Ep[idx], 1)) {
      {
	fprintf(stderr, "... calculating Eigenstates for %s:%d\n", 
		mbfile[i], idxii);
	initEigenstates(&P, &Ep[idx], 1);
	calcEigenstates(&P, &Int, &obsme[i+i*n][idxii+idxii*N], &Ep[idx], threshkmix);
	calcexpectprojectedMBME(&P, &OpObservables, &obsme[i+i*n][idxii+idxii*N], &Sp, &Ep[idx], obsp);
	sortEigenstates(&P, &Int, obsp, &Ep[idx], minnormkmix, allp);
      }
    }

  // solve the n SlaterDet eigenvalue problem
  Eigenstates multiE;
  Amplitudes multiA;
  initEigenstates(&P, &multiE, dim);
  initAmplitudes(&P, &multiA, dim);
  calcMultiEigenstatesMulti(&P, &Int, obsme, 
			    Ep, In, n, 
			    &multiE, &multiA, threshmulti);

  // calculate expectation values
  Observablesod** obs = initprojectedVector(&P, &OpObservables, dim);
  calcexpectprojectedMultiMBME(&P, &OpObservables, obsme, 
			       Q, In, n, &multiE, obs);

  // sort Eigenstates
  sortEigenstates(&P, &Int, obs, &multiE, minnormmulti, all);

  // output
  FILE *outfp;
  char outfilename[255];

  snprintf(outfilename, 255, "%s.multi-%s", nucsfile, ProjectiontoStr(&P));
  if (!(outfp = fopen(outfilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfilename);
    return -1;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);
  fprintf(outfp, "# K-mixing,        threshold: %g, minnorm: %g\n", 
          threshkmix, minnormkmix);
  fprintf(outfp, "# Diagonalization, threshold: %g, minnorm: %g\n",
          threshmulti, minnormmulti);

  // output Observables

  // only needed for output
  SlaterDet Q0;
  allocateSlaterDet(&Q0, Q[0].A);
  Q[0].get(&Q[0], 0, &Q0);
  showprojectedObservables(outfp, &P, &Int, &Q0, obs, &multiE, &multiA, "");

  fclose(outfp);

  writeMultiMulticonfigfile(outfilename, &P, &mbfile, In, n, &multiE);

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


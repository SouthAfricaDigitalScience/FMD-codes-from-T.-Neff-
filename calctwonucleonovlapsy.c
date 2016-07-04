/**

  \file calctwonucleonovlapsy.c

  calculate two-nucleon ovlaps in Y-coordinates


  (c) 2012 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/TwoNucleonOvlaps.h"
#include "fmd/Projection.h"

#include "numerics/zcw.h"
#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/Projectionmpi.h"
#include "fmdmpi/TwoNucleonOvlapsSlave.h"
#endif


#define ABS(x) ((x) < 0 ? -(x) : (x))


int triangle(int jf, int k, int ji)
{
  if (ABS(ji-k) <= jf && jf <= ji+k)
    return 1;
  else
    return 0;
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
    TwoNucleonOvlapsSlave();
  } else {
#endif

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] J,Pi,IDX,MBSTATE J,Pi,IDX,MBSTATE"
	    "\n     -q QMAX         calculate to max momentum"
	    "\n     -n NPOINTS      number of points"
	    "\n     -p NALPHA-NBETA number of angles"
	    "\n     -r              recoil\n",
	    argv[0]);
    exit(-1);
  }

  // spectroscopic amplitude parameters

  int npoints = 51;
  double qmax = 5.0;
  int nalpha = 5;
  int ncosb = 4;
  int recoil = 0;

  // project to this states
  int jA, pA, aA;
  int jB, pB, aB;

  /* manage command-line options */

  char* projparas;

  char c;
  while ((c = getopt(argc, argv, "q:n:p:r")) != -1)
    switch (c) {
    case 'q':
      qmax = atof(optarg);
      break;
    case 'n':
      npoints = atoi(optarg);
      break;
    case 'p':
      projparas = optarg;
      char* ang = strtok(projparas, "-");
      nalpha = atoi(ang);
      ang = strtok(NULL, "-");
      ncosb = atoi(ang);
      break;
    case 'r':
      recoil = 1;
      break;
    }

  char* Aparas = argv[optind];
  char* Bparas = argv[optind+1];

  char* p = strtok(Aparas, ",");
  jA = atoi(p);
  p = strtok(NULL, ",");
  pA = (*p == '+' ? 0 : 1);
  p = strtok(NULL, ",");
  aA = atoi(p);
  p = strtok(NULL, ",");
  char* mcAstatefile = p;

  p = strtok(Bparas, ",");
  jB = atoi(p);
  p = strtok(NULL, ",");
  pB = (*p == '+' ? 0 : 1);
  p = strtok(NULL, ",");
  aB = atoi(p);
  p = strtok(NULL, ",");
  char* mcBstatefile = p;

  // open multiconfigfiles
  Projection PA;
  Projection PB;

  SlaterDet* QA;
  Symmetry* SA;
  Eigenstates EA;
  int nA;

  SlaterDet* QB;
  Symmetry* SB;
  Eigenstates EB;
  int nB;

  char** mbAfile;
  if (readMulticonfigfile(mcAstatefile, &mbAfile, &PA, 
			  &QA, &SA, &EA, &nA)) {
    fprintf(stderr, "couldn't open %s for reading\n", mcAstatefile);
    exit(-1);
  }

  char** mbBfile;
  if (readMulticonfigfile(mcBstatefile, &mbBfile, &PB, 
			  &QB, &SB, &EB, &nB)) {
    fprintf(stderr, "couldn't open %s for reading\n", mcBstatefile);
    exit(-1);
  }

  // PA and PB should be identical besides of odd !
  Projection P = PB;


  TwoNucleonOvlapsPara SAP = {
    qmax : qmax,
    npoints : npoints,
    nalpha : nalpha,
    ncosb : ncosb,
    recoil : recoil
  };

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastTwoNucleonOvlapsPara(&SAP);
  BroadcastA(&QB->A);
#endif

  int i;
  int a,b;

  // which spectroscopic amplitudes are possible ?
  // THIS IS NOT GENERIC YET, ONLY OK FOR SCALAR OVLAPS

  int multi[NSPECY] = {0, 0, 0, 0, 0};

  fprintf(stderr, "jA: %d, pA: %d\n", jA, pA);
  fprintf(stderr, "jB: %d, pB: %d\n", jB, pB);

  // only scalar ovlaps implemented for now
  if (pB == pA && jB == jA) {
    for (i=0; i<NSPECY; i++)
      multi[i] =1;
  } 

  initOpTwoNucleonOvlapsY(&SAP);

  void* specamplitudeme[NSPECY][nA*nB];

  // read or calculate matrix elements

  int havetocalc;

  for (b=0; b<nB; b++)
    for (a=0; a<nA; a++) {

      for (i=0; i<NSPECY; i++) 
	specamplitudeme[i][a+b*nA] = 
	  initprojectedMBME(&P, &OpTwoNucleonOvlapY[i]);

      havetocalc = 0;
      for (i=0; i<NSPECY; i++)
	if (multi[i])
	  if (readprojectedMBMEfromFile(mbAfile[a], mbBfile[b], &P,
					&OpTwoNucleonOvlapY[i], 
					SA[a], SB[b],
					specamplitudeme[i][a+b*nA]))
	    havetocalc = 1;

      if (havetocalc) {
	void* me[NSPECY];
	for (i=0; i<NSPECY; i++)
	  me[i] = specamplitudeme[i][a+b*nA];

#ifdef MPI
	calcprojectedMBMEsmpi(&P, &OpTwoNucleonOvlapsY, &QA[a], &QB[b],
			      SA[a], SB[b], me);
#else
	calcprojectedMBMEs(&P, &OpTwoNucleonOvlapsY, &QA[a], &QB[b],
			   SA[a], SB[b], me);
#endif
	for (i=0; i<NSPECY; i++)
	  writeprojectedMBMEtoFile(mbAfile[a], mbBfile[b], &P, 
				   &OpTwoNucleonOvlapY[i], SA[a], SB[b], 
				   specamplitudeme[i][a+b*nA]);
      }
    }


  void* specamplitude[NSPECY];

  for (i=0; i<NSPECY; i++)   
    if (multi[i]) {
      specamplitude[i] = 
	initprojectedtransitionVector(&P, &OpTwoNucleonOvlapY[i], nA, nB);
      calctransitionprojectedMBME(&P, &OpTwoNucleonOvlapY[i], 
				  specamplitudeme[i], 
				  SA, SB, &EA, &EB, specamplitude[i]);
    }

  // output

  char outfile[255];

  snprintf(outfile, 255, "%s-%s.%d-%s-%s.%d.twonucleonovlapsy", 
	   stripstr(mcAstatefile, ".states"), AngmomtoStr(jA, pA), aA,
	   stripstr(mcBstatefile, ".states"), AngmomtoStr(jB, pB), aB);

  // write spectroscopic factors to data files

  char* speclabel[NSPECY] = { "0s12s12", "0p12p12", "0p32p32", "0d32d32", "0d52d52" };

  for (i=0; i<NSPECY; i++)
    if (multi[i]) {
      
      char datafile[255];
      FILE *datafp;

      snprintf(datafile, 255, "%s-%s--%05.2f-%d-%d-%d%s.data", 
	       outfile, speclabel[i],
	       qmax, npoints, nalpha, ncosb, recoil ? "-recoil" : "");
      if (!(datafp = fopen(datafile, "w"))) {
	fprintf(stderr, "couldn't open %s for writing\n", datafile);
	exit(-1);
      }

      writeTwoNucleonOvlapsY(datafp, 
			     &P, &SAP,
			     jA, pA, aA,
			     jB, pB, aB,
			     specamplitude[i],
			     &EA, &EB);
      fclose(datafp);
    }

  cleanup:
#ifdef MPI
  task=TASKFIN;
  BroadcastTask(&task);
  }

  MPI_Finalize();
#endif

  return 0;
}

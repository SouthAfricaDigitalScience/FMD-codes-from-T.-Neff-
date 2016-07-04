/**

  \file calconenucleonovlaps.c

  calculate spectroscopic amplitudes


  (c) 2005-2007 Thomas Neff, Benjamin Hellwig

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/OneNucleonOvlaps.h"
#include "fmd/Projection.h"

#include "numerics/zcw.h"
#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/Projectionmpi.h"
#include "fmdmpi/OneNucleonOvlapsSlave.h"
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
    OneNucleonOvlapsSlave();
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


  OneNucleonOvlapsPara SAP = {
    qmax : qmax,
    npoints : npoints,
    nalpha : nalpha,
    ncosb : ncosb,
    recoil : recoil
  };

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastOneNucleonOvlapsPara(&SAP);
  BroadcastA(&QB->A);
#endif


  int i;
  int a,b;

  // which spectroscopic amplitudes are possible ?

  int multi[NSPEC] = {0, 0, 0, 0, 0};

  fprintf(stderr, "jA: %d, pA: %d\n", jA, pA);
  fprintf(stderr, "jB: %d, pB: %d\n", jB, pB);

  if (pB == pA && triangle(jA, 1, jB))
    multi[0] = 1;

  if (pB != pA && triangle(jA, 1, jB))
    multi[1] = 1;

  if (pB != pA && triangle(jA, 3, jB))
    multi[2] = 1;

  if (pB == pA && triangle(jA, 3, jB))
    multi[3] = 1;

  if (pB == pA && triangle(jA, 5, jB))
    multi[4] = 1;


  initOpOneNucleonOvlaps(&SAP);

  void* specamplitudeme[NSPEC][nA*nB];

  // read or calculate matrix elements

  int havetocalc;

  for (b=0; b<nB; b++)
    for (a=0; a<nA; a++) {

      for (i=0; i<NSPEC; i++) 
	specamplitudeme[i][a+b*nA] = 
	  initprojectedMBME(&P, &OpOneNucleonOvlap[i]);

      havetocalc = 0;
      for (i=0; i<NSPEC; i++)
	if (multi[i])
	  if (readprojectedMBMEfromFile(mbAfile[a], mbBfile[b], &P,
					&OpOneNucleonOvlap[i], 
					SA[a], SB[b],
					specamplitudeme[i][a+b*nA]))
	    havetocalc = 1;

      if (havetocalc) {
	void* me[NSPEC];
	for (i=0; i<NSPEC; i++)
	  me[i] = specamplitudeme[i][a+b*nA];

#ifdef MPI
	calcprojectedMBMEsmpi(&P, &OpOneNucleonOvlaps, &QA[a], &QB[b],
			      SA[a], SB[b], me);
#else
	calcprojectedMBMEs(&P, &OpOneNucleonOvlaps, &QA[a], &QB[b],
			   SA[a], SB[b], me);
#endif
	for (i=0; i<NSPEC; i++)
	  writeprojectedMBMEtoFile(mbAfile[a], mbBfile[b], &P, 
				   &OpOneNucleonOvlap[i], SA[a], SB[b], 
				   specamplitudeme[i][a+b*nA]);
      }
    }


  void* specamplitude[NSPEC];

  for (i=0; i<NSPEC; i++)   
    if (multi[i]) {
      specamplitude[i] = 
	initprojectedtransitionVector(&P, &OpOneNucleonOvlap[i], nA, nB);
      calctransitionprojectedMBME(&P, &OpOneNucleonOvlap[i], 
				  specamplitudeme[i], 
				  SA, SB, &EA, &EB, specamplitude[i]);
    }

  // output

  char outfile[255];

  snprintf(outfile, 255, "%s-%s.%d-%s-%s.%d.onenucleonovlaps", 
	   stripstr(mcAstatefile, ".states"), AngmomtoStr(jA, pA), aA,
	   stripstr(mcBstatefile, ".states"), AngmomtoStr(jB, pB), aB);

  // write spectroscopic factors to data files

  for (i=0; i<NSPEC; i++)
    if (multi[i]) {
      int j,l;
      switch(i) {
      case 0:
	j=1; l=0;
	break;
      case 1:
	j=1; l=2;
	break;
      case 2:
	j=3; l=2;
	break;
      case 3:
	j=3; l=4;
	break;
      case 4:
	j=5; l=4;
	break;
      }
      
      char datafile[255];
      FILE *datafp;

      snprintf(datafile, 255, "%s-%d-%d--%05.2f-%d-%d-%d%s.data", 
	       outfile, j, l,
	       qmax, npoints, nalpha, ncosb, recoil ? "-recoil" : "");
      if (!(datafp = fopen(datafile, "w"))) {
	fprintf(stderr, "couldn't open %s for writing\n", datafile);
	exit(-1);
      }

      writeOneNucleonOvlaps(datafp, 
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

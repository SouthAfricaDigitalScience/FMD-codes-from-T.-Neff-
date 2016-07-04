/**

  \file calcoccupationnumbershoproj.c

  calculate harmonic oscillator occupation numbers


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



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] mcstate"
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

  char* mcstatefile = argv[optind];
  char** mbfile;

  // open multiconfigfile
  Projection P;
  SlaterDet* Q;
  Symmetry* S;
  Eigenstates E;
  int n;

  if (readMulticonfigfile(mcstatefile, &mbfile, &P, &Q, &S, &E, &n)) {
    fprintf(stderr, "couldn't open %s for reading\n", mcstatefile);
    exit(-1);
  }

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

  int i;
  int a,b;

  initHOBasis(nmax);
  initOpDiagonalDensityMatrixHO(&DMpar);

  void* dmme[n*n];

  // read or calculate matrix elements

  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      dmme[a+b*n] = initprojectedMBME(&P, &OpDiagonalDensityMatrixHO);

	  if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P,
					&OpDiagonalDensityMatrixHO, S[a], S[b],
					dmme[a+b*n])) {
	    calcprojectedMBME(&P, &OpDiagonalDensityMatrixHO, &Q[a], &Q[b],
			      S[a], S[b], dmme[a+b*n]);
	    writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				     &OpDiagonalDensityMatrixHO, S[a], S[b], 
				     dmme[a+b*n]);
	  }
    }	

  fprintf(stderr, "calculate occupation numbers\n");

  void* dmexp = initprojectedVector(&P, &OpDiagonalDensityMatrixHO, n);
  calcexpectprojectedMBME(&P, &OpDiagonalDensityMatrixHO, dmme, 
			  S, &E, dmexp);

  // output

  char outfile[255];
  FILE* outfp;

  snprintf(outfile, 255, "%s.occuho-%d-%05.2f", 
	   stripstr(mcstatefile, ".states"), nmax, hbc*omega);
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    exit(-1);
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedDiagonalDensityMatricesHO(outfp, &P, &DMpar, dmexp, &E);

  fclose(outfp);

  exit(0);
}

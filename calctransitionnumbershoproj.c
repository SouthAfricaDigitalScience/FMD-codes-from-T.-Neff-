/**

  \file calctransitionnumbershoproj.c

  calculate transition density in harmonic oscillator basis


  (c) 2006 Thomas Neff

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
	    "\n   -j J          2j of initial state"
	    "\n   -p +1|-1      parity of initial state"
	    "\n   -a INDEX      index of initial state"
	    "\n   -b INDEX      index of final state"	    
	    "\n   -n NMAX	oscillator cut off"
	    "\n   -o OMEGA	Oscillator constant [MeV]\n",
	    argv[0]);
    exit(-1);
  }

  int nmax = HONMAX;
  double omega = 0.0;

  // project to this state
  int j=-1, pi=-1, alpha=0, beta=0;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "j:p:a:b:n:o:")) != -1)
    switch (c) {
    case 'j':
      j = atoi(optarg);
      break;
    case 'p':
      pi = atoi(optarg) == -1 ? 1 : 0;
      break;
    case 'a':
      alpha = atoi(optarg);
      break;
    case 'b':
      beta = atoi(optarg);
      break;
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

  // J or Pi unset
  if (j == -1) j = (Q[0].A % 2 ? 1 : 0);
  if (pi == -1) pi = (Q[0].A % 2 ? 1 : 0); 

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

  fprintf(stderr, "calculate transition densities\n");

  void* dmtrans = initprojectedtransitionVectornull(&P, &OpDiagonalDensityMatrixHO, n, n);
  calctransitionprojectedMBMEipj(&P, &OpDiagonalDensityMatrixHO, dmme, 
				 S, S, &E, &E,
				 j, pi, beta, j, pi, alpha,
				 dmtrans);

  // output

  char outfile[255];
  FILE* outfp;

  snprintf(outfile, 255, "%s-%s.%d-%s.%d.transho-%d-%05.2f", 
	   stripstr(mcstatefile, ".states"), 
	   AngmomtoStr(j, pi), alpha,
	   AngmomtoStr(j, pi), beta,
	   nmax, hbc*omega);

  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    exit(-1);
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedDiagonalTransitionDensityMatrixHO(outfp, &P, &DMpar,
						  j, pi, beta, alpha,
						  dmtrans, &E);

  fclose(outfp);

  return 0;
}

/**

  \file calcprojectedoccupationnumbershome.c

  calculate harmonic oscillator occupation numbers


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

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA [SYMMETRY:]MBSTATE [SYMMETRY:]MBSTATE"
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

  char* projpar = argv[optind];
  char** mbfile = &argv[optind+1];


  SlaterDet Q[2];
  Symmetry S[2];

  int i;
  for (i=0; i<2; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      exit(-1);
  }

  // odd numer of nucleons ?
  int odd = Q[0].A % 2;


  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

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

  initHOBasis(nmax);
  initOpDiagonalDensityMatrixHO(&DMpar);

  void* dmme;

  // read or calculate matrix elements

  dmme = initprojectedMBME(&P, &OpDiagonalDensityMatrixHO);

  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P,
				&OpDiagonalDensityMatrixHO, S[0], S[1],
				dmme)) {
    calcprojectedMBME(&P, &OpDiagonalDensityMatrixHO, &Q[0], &Q[1],
		      S[0], S[1], dmme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
			     &OpDiagonalDensityMatrixHO, S[0], S[1], 
			     dmme);
  }

  return 0;
}

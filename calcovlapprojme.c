/**

  \file calcovlapprojme.c

  calculate projected many-body overlap matrixelement


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"
#include "fmd/Ovlap.h"

#include "misc/utils.h"
#include "misc/physics.h"



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s PROJPARA MBSTATE1 MBSTATE2\n",
	    filepart(argv[0]));
    goto cleanup;
  }

  int odd;

  char* projpar = argv[optind];
  char** mbfile = &argv[optind+1];

  SlaterDet Q[2];
  Symmetry S[2];

  int i;
  for (i=0; i<2; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      goto cleanup;
  }

  // odd numer of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  complex double** ovlme = initprojectedMBME(&P, &OpOvlap);
  
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, &OpOvlap,
				S[0], S[1], ovlme)) {
    calcprojectedMBME(&P, &OpOvlap, &Q[0], &Q[1], S[0], S[1], ovlme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P,
			     &OpOvlap, S[0], S[1], ovlme);
  }
  else 
    fprintf(stderr, "... matrix elements exist already\n");

 cleanup:
  return 0;
}



/**

  \file calcradiiallme.c

  calculate radii matrix elements for projected states


  (c) 2009 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/RadiiAll.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA [SYMMETRY:]MBSTATE [SYMMETRY:]MBSTATE\n", argv[0]);
    exit(-1);
  }

  /* manage command-line options */

  char* projpar = argv[optind];
  char** mbfile = &argv[optind+1];

  SlaterDet Q[2];
  Symmetry S[2];
  int odd;

  int i;
  for (i=0; i<2; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      exit(-1);
  }

  // odd numer of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  void* radiime = initprojectedMBME(&P, &OpRadiiAll); 

  // read or calculate matrix elements

  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, 
                                &OpRadiiAll, S[0], S[1], 
                                radiime)) {
    calcprojectedMBME(&P, &OpRadiiAll, &Q[0], &Q[1], 
                      S[0], S[1], radiime);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
                             &OpRadiiAll, S[0], S[1], 
                             radiime);
  } 

  return 0;
}



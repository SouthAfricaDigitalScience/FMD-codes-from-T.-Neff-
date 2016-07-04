/**

  \file calcenergyprojme.c

  calculate projected many-body matrixelement


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


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 5) {
    fprintf(stderr, "\nusage: %s PROJPARA INTERACTION [SYMMETRY:]MBSTATE [SYMMETRY:]MBSTATE\n",
	    filepart(argv[0]));
    exit(-1);
  }

  int odd;

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char** mbfile = &argv[optind+2];

  SlaterDet Q[2];
  Symmetry S[2];

  int i;
  for (i=0; i<2; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      exit(-1);
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    exit(-1);
  Int.cm = 1;

  // odd numer of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  initOpObservables(&Int);


  Observablesod** obsme = initprojectedMBME(&P, &OpObservables);
  
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, &OpObservables, 
				S[0], S[1], obsme)) {
    calcprojectedMBME(&P, &OpObservables, &Q[0], &Q[1], S[0], S[1], obsme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
			     &OpObservables, S[0], S[1], obsme);
  }
  else 
    fprintf(stderr, "... matrix elements have already been calculated\n");

}

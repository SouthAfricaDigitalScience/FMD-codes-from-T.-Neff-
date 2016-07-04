/**

  \file calcenergyprojmemulti.c

  calculate projected many-body matrixelements


  (c) 2006-2009 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/MultiSlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/Projection.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 5) {
    fprintf(stderr, "\nusage: %s PROJPARA INTERACTION MBSTATE MBSTATE\n",
	    filepart(argv[0]));
    exit(-1);
  }

  int odd;

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char** mbfile = &argv[optind+2];

  MultiSlaterDet MSD[2];
  Indices In[2];

  int i;
  for (i=0; i<2; i++) {
    extractIndicesfromString(&mbfile, &In[i]);
    if (readMultiSlaterDetfromFile(&MSD[i], &In[i], mbfile[i]))
      exit(-1);
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    exit(-1);
  Int.cm = 1;

  // odd numer of nucleons ?
  odd = MSD[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  _setangkappacrit(25.0);

  initOpObservables(&Int);

  // initialize space for matrix elements
  Observablesod*** obsme = initprojectedMultiMBME(&P, &OpObservables,
						  &MSD[0], &MSD[1]);
  
  if (readprojectedMultiMBMEfromFile(mbfile[0], mbfile[1], &MSD[0], &MSD[1],
				     &P, &OpObservables, obsme)) {
    calcprojectedMultiMBME(&P, &OpObservables, &MSD[0], &MSD[1], obsme);
    writeprojectedMultiMBMEtoFile(mbfile[0], mbfile[1], &MSD[0], &MSD[1], 
				  &P, &OpObservables, obsme);
  }
  else 
    fprintf(stderr, "... matrix elements have already been calculated\n");

  return 0;
}

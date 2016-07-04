/**

  \file calcenergyprojmemulti.c

  calculate projected many-body matrixelements


  (c) 2003, 2004, 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/MultiSlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Projection.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s PROJPARA MBSTATE MBSTATE\n",
	    filepart(argv[0]));
    exit(-1);
  }

  int odd;

  char* projpar = argv[optind];
  char** mbfile = &argv[optind+1];

  MultiSlaterDet MSD[2];
  Indices In[2];

  // strip possible indices from filename 
  int i;
  for (i=0; i<2; i++) {
    extractIndicesfromString(&mbfile[i], &In[i]);
    if (readMultiSlaterDetfromFile(&MSD[i], &In[i], mbfile[i]))
      exit(-1);
  }

  // odd numer of nucleons ?
  odd = MSD[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  _setangkappacrit(25.0);

  // initialize space for matrix elements
  complex double*** ovlme = initprojectedMultiMBME(&P, &OpOvlap,
						  &MSD[0], &MSD[1]);
  
  if (readprojectedMultiMBMEfromFile(mbfile[0], mbfile[1], &MSD[0], &MSD[1],
				     &P, &OpOvlap, ovlme)) {
    calcprojectedMultiMBME(&P, &OpOvlap, &MSD[0], &MSD[1], ovlme);
    writeprojectedMultiMBMEtoFile(mbfile[0], mbfile[1], &MSD[0], &MSD[1], 
				  &P, &OpOvlap, ovlme);
  }
  else 
    fprintf(stderr, "... matrix elements have already been calculated\n");

  return 0;
}

/**

  \file calctransprojme.c

  calculate projected many-body matrixelement of em transition operators


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/ElectroMagneticMultipole.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s PROJPARA [SYMMETRY:]MBSTATE [SYMMETRY:]MBSTATE\n",
	    argv[0]);
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

  void* emome   = initprojectedMBME(&P, &OpEMonopole); 
  void* edipme  = initprojectedMBME(&P, &OpEDipole); 
  void* mdipme  = initprojectedMBME(&P, &OpMDipole);
  void* equadme = initprojectedMBME(&P, &OpEQuadrupole);

  // read or calculate matrix elements
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, 
				&OpEMonopole, S[0], S[1], emome)) {
    calcprojectedMBME(&P, &OpEMonopole, &Q[0], &Q[1], 
		      S[0], S[1], emome);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
				 &OpEMonopole, S[0], S[1], emome);
  }
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, 
				&OpEDipole, S[0], S[1], edipme)) {
    calcprojectedMBME(&P, &OpEDipole, &Q[0], &Q[1], 
		      S[0], S[1], edipme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
			     &OpEDipole, S[0], S[1], edipme);
  }
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, 
				&OpMDipole, S[0], S[1], mdipme)) {
    calcprojectedMBME(&P, &OpMDipole, &Q[0], &Q[1], 
		      S[0], S[1], mdipme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
			     &OpMDipole, S[0], S[1], mdipme);
  }
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, 
				&OpEQuadrupole, S[0], S[1], equadme)) {
    calcprojectedMBME(&P, &OpEQuadrupole, &Q[0], &Q[1], 
		      S[0], S[1], equadme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
			     &OpEQuadrupole, S[0], S[1], equadme);
  }

 cleanup:

  return 0;
}


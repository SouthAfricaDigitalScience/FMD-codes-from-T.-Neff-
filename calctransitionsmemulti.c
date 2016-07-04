/**

  \file calctransitionsmemulti.c

  calculate projected many-body matrixelement of em transition operators


  (c) 2006 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/MultiSlaterDet.h"
#include "fmd/ElectroMagneticMultipole.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s PROJPARA MBSTATE MBSTATE\n",
	    argv[0]);
    exit(-1);
  }

  int odd;

  char* projpar = argv[optind];
  char** mbfile = &argv[optind+1];

  MultiSlaterDet Q[2];
  Indices In[2];

  int i;
  for (i=0; i<2; i++) {
    extractIndicesfromString(&mbfile, &In[i]);
    if (readMultiSlaterDetfromFile(&Q[i], &In[i], mbfile[i]))
      exit(-1);
  }

  // odd numer of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  void* emome   = initprojectedMultiMBME(&P, &OpEMonopole, &Q[0], &Q[1]); 
  void* edipme  = initprojectedMultiMBME(&P, &OpEDipole, &Q[0], &Q[1]); 
  void* mdipme  = initprojectedMultiMBME(&P, &OpMDipole, &Q[0], &Q[1]);
  void* equadme = initprojectedMultiMBME(&P, &OpEQuadrupole, &Q[0], &Q[1]);

  // read or calculate matrix elements

  if (readprojectedMultiMBMEfromFile(mbfile[0], mbfile[1], &Q[0], &Q[1],
				     &P, &OpEMonopole, emome)) {
    calcprojectedMultiMBME(&P, &OpEMonopole, &Q[0], &Q[1], emome);
    writeprojectedMultiMBMEtoFile(mbfile[0], mbfile[1], &Q[0], &Q[1], 
				  &P, &OpEMonopole, emome);
  }
  if (readprojectedMultiMBMEfromFile(mbfile[0], mbfile[1], &Q[0], &Q[1],
				     &P, &OpEDipole, edipme)) {
    calcprojectedMultiMBME(&P, &OpEDipole, &Q[0], &Q[1], edipme);
    writeprojectedMultiMBMEtoFile(mbfile[0], mbfile[1], &Q[0], &Q[1], 
				  &P, &OpEDipole, edipme);
  }
  if (readprojectedMultiMBMEfromFile(mbfile[0], mbfile[1], &Q[0], &Q[1],
				     &P, &OpMDipole, mdipme)) {
    calcprojectedMultiMBME(&P, &OpMDipole, &Q[0], &Q[1], mdipme);
    writeprojectedMultiMBMEtoFile(mbfile[0], mbfile[1], &Q[0], &Q[1], 
				  &P, &OpMDipole, mdipme);
  }
  if (readprojectedMultiMBMEfromFile(mbfile[0], mbfile[1], &Q[0], &Q[1],
				     &P, &OpEQuadrupole, equadme)) {
    calcprojectedMultiMBME(&P, &OpEQuadrupole, &Q[0], &Q[1], equadme);
    writeprojectedMultiMBMEtoFile(mbfile[0], mbfile[1], &Q[0], &Q[1], 
				  &P, &OpEQuadrupole, equadme);  
  }  

}


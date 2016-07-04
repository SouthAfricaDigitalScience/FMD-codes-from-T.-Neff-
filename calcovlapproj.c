/**

  \file calcovlappro.c

  calculate projected many-body overlap


  (c) 2007 Thomas Neff

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

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s MBSTATEF MBSTATEI\n",
	    filepart(argv[0]));
    exit(-1);
  }

  char* mcstatefilef = argv[optind];
  char* mcstatefilei = argv[optind+1];
  char** mbfilef; char** mbfilei;

  // openm multiconfigfiles 
  Projection P;
  SlaterDet *Qf, *Qi;
  Symmetry *Sf, *Si;
  Eigenstates Ef, Ei;
  int nf, ni;

  readMulticonfigfile(mcstatefilef, &mbfilef, &P, &Qf, &Sf, &Ef, &nf);
  readMulticonfigfile(mcstatefilei, &mbfilei, &P, &Qi, &Si, &Ei, &ni);

  // read or calculate overlap matrix elements

  complex double** ovlme[nf*ni]; 

  int a,b;
  for (b=0; b<ni; b++)
    for (a=0; a<nf; a++)
      ovlme[a + b*nf] = initprojectedMBME(&P, &OpOvlap);
  
  for (b=0; b<ni; b++)
    for (a=0; a<nf; a++)

      if (readprojectedMBMEfromFile(mbfilef[a], mbfilei[b], &P, &OpOvlap, Sf[a], Si[b], 
                                    ovlme[a + b*nf])) {
        calcprojectedMBME(&P, &OpOvlap, &Qf[a], &Qi[b], Sf[a], Si[b], ovlme[a + b*nf]);
        writeprojectedMBMEtoFile(mbfilef[a], mbfilei[b], &P, &OpOvlap, Sf[a], Si[b], 
                                 ovlme[a + b*nf]);
      }
  

  // calculate many-body overlaps
  complex double**** ovltrans = initprojectedtransitionVector(&P, &OpOvlap, nf, ni);
  calctransitionprojectedMBME(&P, &OpOvlap, ovlme, Sf, Si, &Ef, &Ei, ovltrans);

  // output

  writeprojectedtransitionOvlap(stdout, &P, ovltrans, &Ef, &Ei);

  return 0;
}



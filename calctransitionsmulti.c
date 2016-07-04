/**

  \file calctransitions.c

  calculate electro-magnetic moments and transitions


  (c) 2003, 2005 Thomas Neff

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

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] mcstate"
	    "\n   -A             show all eigenstates\n", argv[0]);
    exit(-1);
  }

  int all=0;

  char c;

  /* manage command-line options */

  while ((c = getopt(argc, argv, "A")) != -1)
    switch (c) {
    case 'A':
      all=1;
      break;
    }

  char* mcstatefile = argv[optind];
  char** mbfile;

  // open multiconfigfile
  Projection P;
  MultiSlaterDet* Q;
  Indices* In;
  Eigenstates E;
  int n;

  readMultiMulticonfigfile(mcstatefile, &mbfile, &P, &Q, &In, &n, &E);

  fprintf(stderr, "n: %d\n", n);
  for (int i=0; i<n; i++)
    fprintf(stderr, "In[%d].n: %d\n", i, In[i].n);

  int a,b; 

  void* emome[n*n]; 
  void* edipme[n*n]; 
  void* mdipme[n*n]; 
  void* equadme[n*n];

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      emome[a+b*n] = initprojectedMultiMBME(&P, &OpEMonopole, &Q[a], &Q[b]);
      edipme[a+b*n] = initprojectedMultiMBME(&P, &OpEDipole, &Q[a], &Q[b]);
      mdipme[a+b*n] = initprojectedMultiMBME(&P, &OpMDipole, &Q[a], &Q[b]);
      equadme[a+b*n] = initprojectedMultiMBME(&P, &OpEQuadrupole, &Q[a], &Q[b]);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMultiMBMEfromFile(mbfile[a], mbfile[b], &Q[a], &Q[b],
					 &P, &OpEMonopole, emome[a+b*n])) {
	calcprojectedMultiMBME(&P, &OpEMonopole, &Q[a], &Q[b], emome[a+b*n]);
	writeprojectedMultiMBMEtoFile(mbfile[a], mbfile[b], &Q[a], &Q[b], 
				      &P, &OpEMonopole, emome[a+b*n]);
      }
      if (readprojectedMultiMBMEfromFile(mbfile[a], mbfile[b], &Q[a], &Q[b],
					 &P, &OpEDipole, edipme[a+b*n])) {
	calcprojectedMultiMBME(&P, &OpEDipole, &Q[a], &Q[b], edipme[a+b*n]);
	writeprojectedMultiMBMEtoFile(mbfile[a], mbfile[b], &Q[a], &Q[b], 
				      &P, &OpEDipole, edipme[a+b*n]);
      }
      if (readprojectedMultiMBMEfromFile(mbfile[a], mbfile[b], &Q[a], &Q[b],
					 &P, &OpMDipole, mdipme[a+b*n])) {
	calcprojectedMultiMBME(&P, &OpMDipole, &Q[a], &Q[b], mdipme[a+b*n]);
	writeprojectedMultiMBMEtoFile(mbfile[a], mbfile[b], &Q[a], &Q[b], 
				      &P, &OpMDipole, mdipme[a+b*n]);
      }
      if (readprojectedMultiMBMEfromFile(mbfile[a], mbfile[b], &Q[a], &Q[b],
					 &P, &OpEQuadrupole, equadme[a+b*n])) {
	calcprojectedMultiMBME(&P, &OpEQuadrupole, &Q[a], &Q[b], equadme[a+b*n]);
	writeprojectedMultiMBMEtoFile(mbfile[a], mbfile[b], &Q[a], &Q[b], 
				      &P, &OpEQuadrupole, equadme[a+b*n]);
      }
    }


  fprintf(stderr, "calculate exectation values\n");

  // dimension of many-body basis
  // beware: matrix elements are calculated in larger basis

  int dim=0;
  for (int i=0; i<n; i++)
    dim += In[i].n;

  // calculate expectation values
  void* emoexp = initprojectedVector(&P, &OpEMonopole, dim);
  calcexpectprojectedMultiMBME(&P, &OpEMonopole, emome, Q, In, n, &E, emoexp);

  void* edipexp = initprojectedVector(&P, &OpEDipole, dim);
  calcexpectprojectedMultiMBME(&P, &OpEDipole, edipme, Q, In, n, &E, edipexp);

  void* mdipexp = initprojectedVector(&P, &OpMDipole, dim);
  calcexpectprojectedMultiMBME(&P, &OpMDipole, mdipme, Q, In, n, &E, mdipexp);

  void* equadexp = initprojectedVector(&P, &OpEQuadrupole, dim);
  calcexpectprojectedMultiMBME(&P, &OpEQuadrupole, equadme, Q, In, n, &E, equadexp);

  fprintf(stderr, "calculate transition strengths\n");

  // calculate transition strengths
  void* emotrans = initprojectedtransitionVector(&P, &OpEMonopole, dim, dim);
  calctransitionprojectedMultiMBME(&P, &OpEMonopole, emome, Q, Q, In, In, n, n, &E, &E, emotrans);

  void* ediptrans = initprojectedtransitionVector(&P, &OpEDipole, dim, dim);
  calctransitionprojectedMultiMBME(&P, &OpEDipole, edipme, Q, Q, In, In, n, n, &E, &E, ediptrans);

  void* mdiptrans = initprojectedtransitionVector(&P, &OpMDipole, dim, dim);
  calctransitionprojectedMultiMBME(&P, &OpMDipole, mdipme, Q, Q, In, In, n, n, &E, &E, mdiptrans);

  void* equadtrans = initprojectedtransitionVector(&P, &OpEQuadrupole, dim, dim);
  calctransitionprojectedMultiMBME(&P, &OpEQuadrupole, equadme, Q, Q, In, In, n, n, &E, &E, equadtrans);


  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.trans", stripstr(mcstatefile, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedEMMultipoles(outfp,
			     &P, emoexp, edipexp, mdipexp, equadexp, &E);

  writeprojectedtransitionEMMultipoles(outfp,
				       &P, emotrans, ediptrans, mdiptrans, equadtrans, &E);

  fclose(outfp);

 cleanup:

  return 0;
}



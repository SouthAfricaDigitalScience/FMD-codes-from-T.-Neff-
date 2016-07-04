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

#include "fmd/SlaterDet.h"
#include "fmd/ElectroMagneticMultipole.h"
#include "fmd/Projection.h"

#include "numerics/zcw.h"
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
  int hermit=0;

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
  SlaterDet* Q;
  Symmetry* S;
  Eigenstates E;
  int n;

  readMulticonfigfile(mcstatefile, &mbfile, &P, &Q, &S, &E, &n);

  int a,b; 

  void* emome[n*n]; 
  void* edipme[n*n]; 
  void* mdipme[n*n]; 
  void* equadme[n*n];

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      emome[a+b*n] = initprojectedMBME(&P, &OpEMonopole);
      edipme[a+b*n] = initprojectedMBME(&P, &OpEDipole);
      mdipme[a+b*n] = initprojectedMBME(&P, &OpMDipole);
      equadme[a+b*n] = initprojectedMBME(&P, &OpEQuadrupole);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpEMonopole, S[a], S[b], emome[a+b*n])) {
	calcprojectedMBME(&P, &OpEMonopole, &Q[a], &Q[b], 
			  S[a], S[b], emome[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpEMonopole, S[a], S[b], emome[a+b*n]);
      }
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpEDipole, S[a], S[b], edipme[a+b*n])) {
	calcprojectedMBME(&P, &OpEDipole, &Q[a], &Q[b], 
			  S[a], S[b], edipme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpEDipole, S[a], S[b], edipme[a+b*n]);
      }
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpMDipole, S[a], S[b], mdipme[a+b*n])) {
	calcprojectedMBME(&P, &OpMDipole, &Q[a], &Q[b], 
			  S[a], S[b], mdipme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpMDipole, S[a], S[b], mdipme[a+b*n]);
      }
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpEQuadrupole, S[a], S[b], equadme[a+b*n])) {
	calcprojectedMBME(&P, &OpEQuadrupole, &Q[a], &Q[b], 
			  S[a], S[b], equadme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpEQuadrupole, S[a], S[b], equadme[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpEMonopole, emome, n);
    hermitizeprojectedMBME(&P, &OpEDipole, mdipme, n);
    hermitizeprojectedMBME(&P, &OpMDipole, mdipme, n);
    hermitizeprojectedMBME(&P, &OpEQuadrupole, equadme, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* emoexp = initprojectedVector(&P, &OpEMonopole, n);
  calcexpectprojectedMBME(&P, &OpEMonopole, emome, S, &E, emoexp);

  void* edipexp = initprojectedVector(&P, &OpEDipole, n);
  calcexpectprojectedMBME(&P, &OpEDipole, edipme, S, &E, edipexp);

  void* mdipexp = initprojectedVector(&P, &OpMDipole, n);
  calcexpectprojectedMBME(&P, &OpMDipole, mdipme, S, &E, mdipexp);

  void* equadexp = initprojectedVector(&P, &OpEQuadrupole, n);
  calcexpectprojectedMBME(&P, &OpEQuadrupole, equadme, S, &E, equadexp);

  fprintf(stderr, "calculate transition strengths\n");

  // calculate transition strengths
  void* emotrans = initprojectedtransitionVector(&P, &OpEMonopole, n, n);
  calctransitionprojectedMBME(&P, &OpEMonopole, emome, S, S, &E, &E, emotrans);

  void* ediptrans = initprojectedtransitionVector(&P, &OpEDipole, n, n);
  calctransitionprojectedMBME(&P, &OpEDipole, edipme, S, S, &E, &E, ediptrans);

  void* mdiptrans = initprojectedtransitionVector(&P, &OpMDipole, n, n);
  calctransitionprojectedMBME(&P, &OpMDipole, mdipme, S, S, &E, &E, mdiptrans);

  void* equadtrans = initprojectedtransitionVector(&P, &OpEQuadrupole, n, n);
  calctransitionprojectedMBME(&P, &OpEQuadrupole, equadme, S, S, &E, &E, equadtrans);


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



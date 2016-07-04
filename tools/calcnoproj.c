/**

  \file calcnoproj.c

  calculate oscillator quanta for projected states


  (c) 2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/NOsci.h"
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

  void* nome[n*n]; 

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      nome[a+b*n] = initprojectedMBME(&P, &OpNOsci);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpNOsci, S[a], S[b], nome[a+b*n])) {
	calcprojectedMBME(&P, &OpNOsci, &Q[a], &Q[b], 
			  S[a], S[b], nome[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpNOsci, S[a], S[b], nome[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpNOsci, nome, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* noexp = initprojectedVector(&P, &OpNOsci, n);
  calcexpectprojectedMBME(&P, &OpNOsci, nome, S, &E, noexp);

  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.nosci", stripstr(mcstatefile, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedNOsci(outfp, &P, noexp, &E);

  fclose(outfp);

 cleanup:

  return 0;
}



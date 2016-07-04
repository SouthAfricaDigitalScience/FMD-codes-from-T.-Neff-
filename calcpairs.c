/**

  \file calcpairs.c

  calculate two-body pairs in S,T-channels


  (c) 2010 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/TwoBodyDensity.h"
#include "fmd/Projection.h"

#include "misc/utils.h"


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

  void* pairsme[n*n]; 

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      pairsme[a+b*n] = initprojectedMBME(&P, &OpPairs);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpPairs, S[a], S[b], pairsme[a+b*n])) {
	calcprojectedMBME(&P, &OpPairs, &Q[a], &Q[b], 
			  S[a], S[b], pairsme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpPairs, S[a], S[b], pairsme[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpPairs, pairsme, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* pairsexp = initprojectedVector(&P, &OpPairs, n);
  calcexpectprojectedMBME(&P, &OpPairs, pairsme, S, &E, pairsexp);

  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.pairs", stripstr(mcstatefile, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedPairs(outfp, &P, pairsexp, &E);

  fclose(outfp);

 cleanup:

  return 0;
}



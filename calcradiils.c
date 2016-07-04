/**

  \file calcradiils.c

  calculate spin-orbit contribution to charge radii


  (c) 2011 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/RadiiLS.h"
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

  void* radiilsme[n*n]; 

  int a,b;
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      radiilsme[a+b*n] = initprojectedMBME(&P, &OpRadiiLS);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpRadiiLS, S[a], S[b], 
				    radiilsme[a+b*n])) {
	calcprojectedMBME(&P, &OpRadiiLS, &Q[a], &Q[b], 
			  S[a], S[b], radiilsme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpRadiiLS, S[a], S[b], 
				 radiilsme[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpRadiiLS, radiilsme, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* radiilsexp = initprojectedVector(&P, &OpRadiiLS, n);
  calcexpectprojectedMBME(&P, &OpRadiiLS, radiilsme, S, &E, radiilsexp);


  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.radiils", stripstr(mcstatefile, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedRadiiLS(outfp, Q, &P, radiilsexp, &E);

  fclose(outfp);

 cleanup:

  return 0;
}



/**

  \file calcsdradii.c

  calculate radii for projected states
  do NOT correct for center-of-mass


  (c) 2009 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/SDRadii.h"
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

  void* radiime[n*n]; 

  int a,b;
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      radiime[a+b*n] = initprojectedMBME(&P, &OpSDRadii);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpSDRadii, S[a], S[b], 
				    radiime[a+b*n])) {
	calcprojectedMBME(&P, &OpSDRadii, &Q[a], &Q[b], 
			  S[a], S[b], radiime[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpSDRadii, S[a], S[b], 
				 radiime[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpSDRadii, radiime, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* radiiexp = initprojectedVector(&P, &OpSDRadii, n);
  calcexpectprojectedMBME(&P, &OpSDRadii, radiime, S, &E, radiiexp);


  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.radii", stripstr(mcstatefile, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedSDRadii(outfp, Q, &P, radiiexp, &E);

  fclose(outfp);

 cleanup:

  return 0;
}



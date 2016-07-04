/**

  \file calctimereversal.c

  calculate time reversal eigenvalue


  (c) 2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/TimeReversal.h"
#include "fmd/Projection.h"


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

  void* trme[n*n]; 

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      trme[a+b*n] = initprojectedMBME(&P, &OpTimeReversal);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpTimeReversal, S[a], S[b], trme[a+b*n])) {
	calcprojectedMBME(&P, &OpTimeReversal, &Q[a], &Q[b], 
			  S[a], S[b], trme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpTimeReversal, S[a], S[b], trme[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpTimeReversal, trme, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* trexp = initprojectedVector(&P, &OpTimeReversal, n);
  calcexpectprojectedMBME(&P, &OpTimeReversal, trme, S, &E, trexp);

  // output

  fprintinfo(stdout);
  fprintProjectinfo(stdout, &P);

  showprojectedTimeReversal(stdout, &P, trexp, &E);

 cleanup:

  return 0;
}



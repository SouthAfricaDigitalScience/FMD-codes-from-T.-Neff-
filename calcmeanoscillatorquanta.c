/**

  \file calcmeanoscillatorquanta.c

  calculate mean number of oscillator quanta in many-body state


  (c) 2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/MeanOscillatorQuanta.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] mcstate"
	    "\n   -o OMEGA       Oscillator constant [MeV]"
	    "\n   -A             show all eigenstates\n", argv[0]);
    exit(-1);
  }

  int all=0;
  int hermit=0;
  double omega=0.0;

  char c;

  /* manage command-line options */

  while ((c = getopt(argc, argv, "Ao:")) != -1)
    switch (c) {
    case 'o':
      omega = atof(optarg)/hbc;
      break;
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

  void* quantame[n*n]; 

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      quantame[a+b*n] = initprojectedMBME(&P, &OpMeanOsciQuanta);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpMeanOsciQuanta, S[a], S[b], 
				    quantame[a+b*n])) {
	calcprojectedMBME(&P, &OpMeanOsciQuanta, &Q[a], &Q[b], 
			  S[a], S[b], quantame[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpMeanOsciQuanta, S[a], S[b], 
				 quantame[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpMeanOsciQuanta, quantame, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* quantaexp = initprojectedVector(&P, &OpMeanOsciQuanta, n);
  calcexpectprojectedMBME(&P, &OpMeanOsciQuanta, quantame, S, &E, quantaexp);


  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.meanosciquanta", stripstr(mcstatefile, tostrip));
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedMeanOscillatorQuanta(outfp, omega,
				     &P, quantaexp, &E);

  fclose(outfp);

 cleanup:

  return 0;
}



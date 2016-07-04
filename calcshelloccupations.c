/**

  \file calcshelloccupations.c

  calculate distribution into N-hw excitations


  (c) 2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/ProjectedShellOccupations.h"
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
	    "\n   -n NMAX        up to NMAX excitations"
	    "\n   -N NINT        use NINT points for theta integration"
	    "\n   -a             show also proton and neutron occupations"
	    "\n   -A             show all eigenstates\n", argv[0]);
    exit(-1);
  }

  int all=0;
  int hermit=0;
  int totalonly=1;

  double omega=0.0;
  int nmax = 10;
  int nint;

  char c;

  /* manage command-line options */

  while ((c = getopt(argc, argv, "Aao:n:N:")) != -1)
    switch (c) {
    case 'a':
      totalonly = 0;
      break;
    case 'o':
      omega = atof(optarg)/hbc;
      break;
    case 'n':
      nmax = atoi(optarg);
      break;
    case 'N':
      nint = atoi(optarg);
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

  if (!nint)
    nint = 2*nmax;

  ShellOccupationsPara par = {
    omega : omega,
    nmax  : nmax,
    nint  : nint
  };

  initOpShellOccupations(&par);

  void* occupationsme[n*n]; 

  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      occupationsme[a+b*n] = initprojectedMBME(&P, &OpShellOccupations);
    }

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, 
				    &OpShellOccupations, S[a], S[b], 
				    occupationsme[a+b*n])) {
	calcprojectedMBME(&P, &OpShellOccupations, &Q[a], &Q[b], 
			  S[a], S[b], occupationsme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpShellOccupations, S[a], S[b], 
				 occupationsme[a+b*n]);
      }
    }


  if (hermit) {
    hermitizeprojectedMBME(&P, &OpShellOccupations, occupationsme, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // calculate expectation values
  void* occupationsexp = initprojectedVector(&P, &OpShellOccupations, n);
  calcexpectprojectedMBME(&P, &OpShellOccupations, occupationsme, S, &E, occupationsexp);


  // output

  char outfile[255];
  char tostrip[255];
  FILE* outfp;

  snprintf(tostrip, 255, ".states");
  snprintf(outfile, 255, "%s.shelloccupations-%05.2f", 
	   stripstr(mcstatefile, tostrip), hbc*omega);
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    goto cleanup;
  }

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  writeprojectedShellOccupations(outfp, &par,
				 &P, occupationsexp, &E, totalonly);

  fclose(outfp);

 cleanup:

  return 0;
}



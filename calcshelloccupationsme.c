/**

  \file calcshelloccupationsme.c

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

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA [SYMMETRY:]MBSTATE [SYMMETRY:]MBSTATE\n"
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

  char* projpar = argv[optind];
  char** mbfile = &argv[optind+1];

  SlaterDet Q[2];
  Symmetry S[2];

  int i;
  for (i=0; i<2; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      goto cleanup;
  }

  // odd numer of nucleons ?
  int odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);


  ShellOccupationsPara par = {
    omega : omega,
    nmax  : nmax,
    nint  : nint
  };

  initOpShellOccupations(&par);

  void* occupationsme; 

  occupationsme = initprojectedMBME(&P, &OpShellOccupations);

  // read or calculate matrix element
  if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P, 
				&OpShellOccupations, S[0], S[1], 
				occupationsme)) {
    calcprojectedMBME(&P, &OpShellOccupations, &Q[0], &Q[1], 
		      S[0], S[1], occupationsme);
    writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
			     &OpShellOccupations, S[0], S[1], 
			     occupationsme);
  }

 cleanup:

  return 0;
}



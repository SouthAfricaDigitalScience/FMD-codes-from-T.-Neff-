/**

  \file calcformfactorsme.c

  calculate formfactor matrix elements for all multipoles
  not parallelized yet, see calcformfactormes.c
  

  (c) 2004-2010 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"
#include "fmd/Formfactors.h"

#include "misc/utils.h"
#include "misc/physics.h"


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  // enough arguments ?

  if (argc < 4) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA [SYMMETRY:]MBSTATE [SYMMETRY:]MBSTATE\n"
	    "\n     -q QMAX         calculate to max momentum"
	    "\n     -n NPOINTS      number of points"
	    "\n     -p NALPHA-NBETA number of angles"
	    "\n     -l MULTIPOLE    select only l-multipole transition"
	    "\n     -r              recoil"
	    "\n     -t              test only\n",
	    argv[0]);
    goto cleanup;
  }

  int npoints = 50+1;
  int nalpha = 5;
  int ncosb = 4;
  double qmax = 5.0;
  int recoil = 0;
  int test = 0;
  int l = -1;

  char* projparas;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "rq:n:p:l:t")) != -1)
    switch (c) {
    case 'r':
      recoil = 1;
      break;
    case 'q':
      qmax = atof(optarg);
      break;
    case 'n':
      npoints = atoi(optarg);
      break;
    case 'p':
      projparas = optarg;
      char* ang = strtok(projparas, "-");
      nalpha = atoi(ang);
      ang = strtok(NULL, "-");
      ncosb = atoi(ang);
      break;
    case 't':
      test = 1;
      break;
    case 'l':
      l = atoi(optarg);
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


  FormfactorPara FfP = {
    qmax : qmax,
    npoints : npoints,
    nalpha : nalpha,
    ncosb : ncosb, 
    recoil : recoil
  };

  // which multipoles to calculate
  int multi[NMULTIPOLES] = {0, 0, 0, 0};

  if (l != -1) {
    multi[l] = 1;
  } else {
    for (i=0; i<NMULTIPOLES; i++) multi[i] = 1;
  }

  // Project formfactor on angular momentum

  initOpFormfactors(&FfP);

  void* ffactorme[NMULTIPOLES];

  // read or calculate matrix elements

  for (i=0; i<NMULTIPOLES; i++)   
    ffactorme[i] = initprojectedMBME(&P, &OpMultipoleFormfactor[i]);

  int havetocalc = 0;
  for (i=0; i<NMULTIPOLES; i++)   
      if (readprojectedMBMEfromFile(mbfile[0], mbfile[1], &P,
				    &OpMultipoleFormfactor[i], S[0], S[1],
				    ffactorme[i]))
	havetocalc = 1;

  if (havetocalc) {

    // test matrix elements only ? return error code
    if (test) {
      fprintf(stderr, "... missing or corrupt matrix elements\n");
      return -1;
    }

    // calculate all multipoles
    calcprojectedMBMEs(&P, &OpMultipoleFormfactors, &Q[0], &Q[1],
			  S[0], S[1], ffactorme);

    // but only write the ones wanted
    for (i=0; i<NMULTIPOLES; i++)
      if (multi[i])
	writeprojectedMBMEtoFile(mbfile[0], mbfile[1], &P, 
				 &OpMultipoleFormfactor[i], S[0], S[1], 
				 ffactorme[i]);
  } else
    fprintf(stderr, "... matrix elements exist already\n");

 cleanup:
  return 0;
}

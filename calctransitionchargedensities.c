/**

  \file calctransitionchargedensities.c

  calculate transition charge densities


  (c) 2004, 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Formfactors.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"




int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] mcstate"
	    "\n     -f              also write Fourier-Bessel coefficients"
	    "\n     -j J            2j of initial state"
	    "\n     -p +1|-1        parity of initial state"
	    "\n     -l MULTIPOLE    multipole density"
	    "\n     -a INDEX        index of initial state"
	    "\n     -b INDEX        index of final state"
	    "\n     -r              with recoil"
	    "\n     -q QMAX         calculate to max momentum"
	    "\n     -n NPOINTS      number of points"
	    "\n     -P NALPHA-NBETA number of angles\n",
	    argv[0]);
    exit(-1);
  }

  int hermit = 0;

  // formfactor parameters

  int npoints = 51;
  double qmax = 5.0;
  int nalpha = 5; 
  int ncosb = 4;
  int recoil = 0;
  int l=0;
  int fb=0;

  // grid in coordinate space

  double rmax = 10.0;
  double deltar = 0.05;
  int npointsr = (int) (rmax/deltar) + 1;

  // project to this state
  int j=-1, pi=-1, alpha=0, beta=0;

  /* manage command-line options */

  char c; char* projparas;
  while ((c = getopt(argc, argv, "fj:p:a:b:rq:n:P:")) != -1)
    switch (c) {
    case 'f':
      fb = 1;
      break;
    case 'j':
      j = atoi(optarg);
      break;
    case 'p':
      pi = atoi(optarg) == -1 ? 1 : 0;
      break;
    case 'a':
      alpha = atoi(optarg);
      break;
    case 'b':
      beta = atoi(optarg);
      break;
    case 'r':
      recoil = 1;
      break;
    case 'q':
      qmax = atof(optarg);
      break;
    case 'n':
      npoints = atoi(optarg);
      break;
    case 'P':
      projparas = optarg;
      char* ang = strtok(projparas, "-");
      nalpha = atoi(ang);
      ang = strtok(NULL, "-");
      ncosb = atoi(ang);
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

  if (readMulticonfigfile(mcstatefile, &mbfile, &P, &Q, &S, &E, &n)) {
    fprintf(stderr, "couldn't open %s for reading\n", mcstatefile);
    exit(-1);
  }

  FormfactorPara FfP = {
    qmax : qmax,
    npoints : npoints,
    nalpha : nalpha,
    ncosb : ncosb,
    recoil : recoil
  };

  int i;
  int a,b;

  // J or Pi unset
  if (j == -1) j = (Q[0].A % 2 ? 1 : 0);
  if (pi == -1) pi = (Q[0].A % 2 ? 1 : 0); 

  fprintf(stderr, "J: %d, Pi: %c1, alpha: %d, beta: %d\n", 
	  j, pi ? '-' : '+', alpha, beta);

  // calculate certain multipoles

  int nmulti = 1;
  int multi[NMULTIPOLES] = {0, 0, 0, 0};
  multi[l] = 1;

  // Project formfactor on angular momentum

  initOpFormfactors(&FfP);

  void* ffactorme[NMULTIPOLES][n*n];

  // read or calculate matrix elements

  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      for (i=0; i<NMULTIPOLES; i++)   
	if (multi[i]) {
	  ffactorme[i][a+b*n] = initprojectedMBME(&P, &OpMultipoleFormfactor[i]);

	  if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P,
					&OpMultipoleFormfactor[i], S[a], S[b],
					ffactorme[i][a+b*n])) {
	    calcprojectedMBME(&P, &OpMultipoleFormfactor[i], &Q[a], &Q[b],
			      S[a], S[b], ffactorme[i][a+b*n]);
	    writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				     &OpMultipoleFormfactor[i], S[a], S[b], 
				     ffactorme[i][a+b*n]);
	  }
	}
    }

  if (hermit) {
    for (i=0; i<NMULTIPOLES; i++)
      if (multi[i])
	hermitizeprojectedMBME(&P, &OpMultipoleFormfactor[i], ffactorme[i], n);
  }

  fprintf(stderr, "calculate transition values\n");

  // calculate static formfactors
  void* ffactortrans[NMULTIPOLES];

  for (i=0; i<NMULTIPOLES; i++)   
    if (multi[i]) {
      ffactortrans[i] = initprojectedtransitionVectornull(&P, &OpMultipoleFormfactor[i], n, n);
      calctransitionprojectedMBMEipj(&P, &OpMultipoleFormfactor[i], ffactorme[i], 
				     S, S, &E, &E, 
				     j, pi, beta, j, pi, alpha, 
				     ffactortrans[i]);
    }

  // output

  char outfile[255];

  snprintf(outfile, 255, "%s-%s.%d-%s.%d.transchargedens", 
	   stripstr(mcstatefile, ".states"), 
	   AngmomtoStr(j, pi), alpha,
	   AngmomtoStr(j, pi), beta);

  // write densities to data files

  char datafile[255];
  FILE *datafp;

  snprintf(datafile, 255, "%s-%d--%05.2f-%d-%d-%d%s.data", 
	   outfile, l,
	   qmax, npoints, nalpha, ncosb, recoil ? "-recoil" : "");

  if (!(datafp = fopen(datafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", datafile);
    exit(-1);
  }

  for (i=0; i<NMULTIPOLES; i++)
    if (multi[i]) {
      writeTransitionChargeDensities(datafp, &P, &FfP,
				     rmax, npointsr,
				     Q[0].Z, i,
				     j, pi, beta, j, pi, alpha, 
				     ffactortrans[l], &E);
    }

  fclose(datafp);

  // write Fourier-Bessel coefficients

  if (fb) {

    char datafile[255];
    FILE *datafp;

    snprintf(datafile, 255, "%s-%d--%05.2f-%d-%d-%d%s.fb-data", 
	     outfile, l,
	     qmax, npoints, nalpha, ncosb, recoil ? "-recoil" : "");

    if (!(datafp = fopen(datafile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", datafile);
      exit(-1);
    }

    for (i=0; i<NMULTIPOLES; i++)
      if (multi[i]) {
      writeTransitionChargeFourierBessel(datafp, &P, &FfP,
					 rmax, npointsr,
					 Q[0].Z, i,
					 j, pi, beta, j, pi, alpha, 
					 ffactortrans[l], &E);
      }

    fclose(datafp);

  }

  return 0;
}

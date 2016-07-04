/**

  \file calctransitionpointdensities.c

  calculate transition point densities


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
	    "\n     -j jini:jfin    2j of initial and final state"
	    "\n     -p +1|-1        parity of initial state"
	    "\n     -l MULTIPOLE    multipole density"
	    "\n     -a INDEX        index of initial state"
	    "\n     -b INDEX        index of final state"
	    "\n     -q QMAX         calculate to max momentum"
	    "\n     -n NPOINTS      number of points"
	    "\n     -P NALPHA-NBETA number of angles"
	    "\n     -r              with recoil\n",
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

  // grid in coordinate space

  double rmax = 15.0;
  double deltar = 0.05;
  int npointsr = (int) (rmax/deltar) + 1;

  // project to this state
  int jfin=-1, jini=-1, pfin=-1, pini=-1, alpha=0, beta=0;

  /* manage command-line options */

  char c; char* jparas; char* projparas;
  while ((c = getopt(argc, argv, "j:p:a:b:l:rq:n:P:")) != -1)
    switch (c) {
    case 'j':
      jparas = optarg;
      char* j = strtok(jparas, ":");
      jini = atoi(j);
      j = strtok(NULL, ":");
      jfin = atoi(j);
      break;
    case 'p':
      pini = atoi(optarg) == -1 ? 1 : 0;
      break;
    case 'a':
      alpha = atoi(optarg);
      break;
    case 'b':
      beta = atoi(optarg);
      break;
    case 'l':
      l = atoi(optarg);
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

  // Pi unset
  if (pini == -1) pini = (Q[0].A % 2 ? 1 : 0); 
  pfin = (pini+l)%2;

  fprintf(stderr, "Jini: %d%c1, Jfin: %d%c1, alpha: %d, beta: %d\n", 
	  jini, pini ? '-' : '+', jfin, pfin ? '-' : '+', alpha, beta);

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
	    fprintf(stderr, "could not read matrix elements\n");
	    exit(-1);
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
				     jfin, pfin, beta, jini, pini, alpha, 
				     ffactortrans[i]);
    }

  // output

  char outfile[255];

  snprintf(outfile, 255, "%s-%s.%d-%s.%d.transpointdensred", 
	   stripstr(mcstatefile, ".states"), 
	   AngmomtoStr(jini, pini), alpha,
	   AngmomtoStr(jfin, pfin), beta);

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
      writeTransitionPointDensities(datafp, &P, &FfP,
				    rmax, npointsr,
				    Q[0].A, i,
				    jfin, pfin, beta, jini, pini, alpha, 
				    ffactortrans[l], &E);
    }

  fclose(datafp);

  return 0;
}

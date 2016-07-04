/**

  \file calcspindensities.c

  calc spin density cuts in coordinate and momentum space
  call IDL


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "fmd/SlaterDet.h"
#include "fmd/SpinDensities.h"

#include "misc/physics.h"
#include "misc/utils.h"

#define SQR(x)	(x)*(x)


int main(int argc, char* argv[])
{

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "Usage: %s [OPTIONS] state\n"
	    "\n   -g            debug"
	    "\n   -q            quiet"
	    "\n   -l LABEL	use label"
	    "\n   -v {xy|yz|xz} select view"
	    "\n   -c RANGE      coordinate range"
	    "\n   -m RANGE	momentum range\n",
	    argv[0]);
    exit(-1);
  }

  int debug=0, quiet=0;
  int coordinate=0; int momentum=0;
  int view[3] = {1, 1, 1};
  int nplots[2] = {2, 3};

  double xmax=5.5;
  double pmax=3.5;

  int npoints=39+1;

  char* label = NULL;

  char c;
  while ((c = getopt(argc, argv, "c:m:l:gqv:")) != -1)
    switch (c) {
    case 'g':
      debug=1;
      break;
    case 'q':
      quiet=1;
      break;
    case 'l':
      label=optarg;
      break;
    case 'c':
      coordinate = 1;
      xmax = atof(optarg);
      break;
    case 'm':
      momentum = 1;
      pmax = atof(optarg);
      break;
    case 'v':
      view[0]=0; view[1]=0; view[2]=0;
      nplots[1] = 1;
      if (!strcmp(optarg, "xy")) view[2]=1;
      if (!strcmp(optarg, "xz")) view[1]=1;
      if (!strcmp(optarg, "yz")) view[0]=1;
      break;
    }	

  if (!(coordinate || momentum)) {
    coordinate = 1;
    momentum = 1;
  }

  nplots[0] = (coordinate && momentum) ? 2 : 1;

  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  readSlaterDetfromFile(&Q, slaterdetfile);

  if (!label)
    label = nucleusIDLformat(nucleusname(Q.A, Q.Z));
  
  double* dens = (double*) malloc(3*SQR(npoints)*sizeof(double));

  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  // write densities to data files

  char datafile[255];
  FILE *datafp;
  int v;

  snprintf(datafile, 255, "%s.spindens", slaterdetfile);;
  if (!(datafp = fopen(datafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", datafile);
    exit(-1);
  }

  // coordinate space densities
  if (coordinate) {

    for (v=0; v<3; v++)
      if (view[v]) {

	calcSpinDensitiesCoordinate(&Q, &X, v, npoints, xmax, dens);

	int i,j,k;
	for (k=0; k<3; k++)
	  for (j=0; j<npoints; j++) {
	    for (i=0; i<npoints; i++)
	      fprintf(datafp, "%f  ", dens[k+(i+j*npoints)*3]/rho0);
	    fprintf(datafp, "\n");
	  }
      }
  }

  // momentum space densities
  if (momentum) {
    
    for (v=0; v<3; v++)
      if (view[v]) {

	calcSpinDensitiesMomentum(&Q, &X, v, npoints, pmax, dens);

	int i,j,k;
	for (k=0; k<3; k++)
	  for (j=0; j<npoints; j++) {
	    for (i=0; i<npoints; i++)
	      fprintf(datafp, "%f  ", dens[k+(i+j*npoints)*3]/Q.A);
	    fprintf(datafp, "\n");
	  }
      }
  }

  fclose(datafp);

  // writing IDL script

  char scriptfile[255];
  FILE* scriptfp;

  char* fmdhome;

  if (!(fmdhome = getenv("FMD"))) {
    fprintf(stderr, "environment variable FMD not defined\n");
    exit(-1);
  }

  snprintf(scriptfile, 255, "%s.spindens.script", slaterdetfile);
  if (!(scriptfp = fopen(scriptfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", scriptfile);
    exit(-1);
  }

  fprintf(scriptfp, 
	  "; plot spin densities of nuclei\n; written by calcspindensities\n\n");
  fprintf(scriptfp, "dens=fltarr(%d,%d,3,%d)\n", 
	  npoints, npoints, nplots[0]*nplots[1]);

  fprintf(scriptfp, "!path = '%s/lib:' + !path\n", fmdhome);
  fprintf(scriptfp, "openr, unit, '%s', /get_lun\n", datafile);
  fprintf(scriptfp, "readf, unit, dens\n");
  fprintf(scriptfp, "free_lun, unit\n");
  
  fprintf(scriptfp, "\n.run multipost, spindensityplot\n");
  fprintf(scriptfp, 
	  "pos = initpost('%s.spindens.eps',%d,%d, gapx=1.5, gapy=1.5, /color)\n",
	  slaterdetfile, nplots[0], nplots[1]);
  fprintf(scriptfp, "!x.minor = 4\n!y.minor = 4\n");
  fprintf(scriptfp, "loadct, 3\n");

  int i = -1;

  if (coordinate) {
    for (v=0; v<3; v++) {
      if (view[v]) {
	i++;
	fprintf(scriptfp, "!p.position = pos(%d,*)\n", i);
	fprintf(scriptfp, "spindensityplot, dens(*,*,*,%d), %d, %f, %f, %f, %f, $\n",
		i, v, -xmax, xmax, -xmax, xmax);
	fprintf(scriptfp, "\t/cut");
	fprintf(scriptfp, ", /coordinate, /vec, /dens, key = '%s'\n", 
		label); 
      }
    }
  }

  if (momentum) {
    for (v=0; v<3; v++) {
      if (view[v]) {
	i++;
	fprintf(scriptfp, "!p.position = pos(%d,*)\n", i);
	fprintf(scriptfp, "spindensityplot, dens(*,*,*,%d), %d, %f, %f, %f, %f, $\n",
		i, v, -pmax, pmax, -pmax, pmax);
	fprintf(scriptfp, "\t/cut");
	fprintf(scriptfp, ", /momentum, /vec, /dens, key = '%s'\n", 
		label);
      }
    }
  }

  fclose(scriptfp);

  // calling IDL

  char call[255];

  snprintf(call, 255, "idl < %s", scriptfile);
  system(call);

  char epsfile[255];
  snprintf(epsfile, 255, "%s.spindens.eps", slaterdetfile);

  // calling gv

  if (!quiet) {
    snprintf(call, 255, "gv %s &", epsfile);
    system(call);
  }

  // clean up

  if (!debug) {
    remove(datafile); 
    remove(scriptfile);
  }

  return 0;
}


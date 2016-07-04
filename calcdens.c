/**

  \file calcdensities.c

  calc density cuts in coordinate and momentum space
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
#include "fmd/Densities.h"

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
	    "\n   -s            write (radial) densities into file"
	    "\n   -l LABEL	use label"
	    "\n   -r RANGE      radial range"
	    "\n   -c RANGE      coordinate range\n",
	    argv[0]);
    exit(-1);
  }

  int debug=0, quiet=0;
  int save=0;
  int shape=1;
  int view=1;

  double xmax=5.5;
  double rmax=5.4;

  int npoints=40+1;
  int izcw=5;

  char* label = NULL;

  char c;
  while ((c = getopt(argc, argv, "c:r:l:gqsn")) != -1)
    switch (c) {
    case 'g':
      debug=1;
      break;
    case 'q':
      quiet=1;
      break;
    case 's':
      save=1;
      break;
    case 'n':
      shape=0;
      break;
    case 'l':
      label=optarg;
      break;
    case 'c':
      xmax = atof(optarg);
      break;
    case 'r':
      rmax = atof(optarg);
      break;
    }	


  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  readSlaterDetfromFile(&Q, slaterdetfile);

  if (!label)
    label = nucleusIDLformat(nucleusname(Q.A, Q.Z));
  
  double* dens = (double*) malloc(3*SQR(npoints)*sizeof(double));
  double* rdens = (double*) malloc(3*npoints*sizeof(double));

  // write densities to data files

  char datafile[255];
  FILE *datafp;
  int v;

  snprintf(datafile, 255, "%s.rdens", slaterdetfile);;
  if (!(datafp = fopen(datafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", datafile);
    exit(-1);
  }

  int i,j;
  // coordinate space densities
  if (shape) {
    calcDensitiesCoordinate(&Q, view, npoints, xmax, dens);

    int axis=2;
    for (j=0; j<npoints; j++) {
      for (i=0; i<npoints; i++)
	fprintf(datafp, "%f  ", dens[j*npoints+i]/rho0);
      fprintf(datafp, "\n");
    }
  }

  // radial densities
  calcRadialDensitiesCoordinate(&Q, izcw, npoints, rmax, rdens);
  
  for (i=0; i<npoints; i++)
    fprintf(datafp, "%15.8g    ", rdens[i]);
  fprintf(datafp, "\n");
  for (i=0; i<npoints; i++)
    fprintf(datafp, "%15.8g    ", rdens[i+  npoints]);
  fprintf(datafp, "\n");
  for (i=0; i<npoints; i++)
    fprintf(datafp, "%15.8g    ", rdens[i+2*npoints]);
  fprintf(datafp, "\n");

  fclose(datafp);

  // write radial densities into file
  if (save) {
    char rdensfile[255];
    FILE* rdensfp;

    snprintf(rdensfile, 255, "%s.dens.dat", slaterdetfile);
    if (!(rdensfp = fopen(rdensfile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", rdensfile);
      exit(-1);
    }  
    
    for (i=0; i<npoints; i++)
      fprintf(rdensfp, "%8.3f\t%15.8g\t%15.8g\t%15.8g\n",
	      rmax/(npoints-1)*i,
	      rdens[i],
	      rdens[i+  npoints],
	      rdens[i+2*npoints]);
    
    fclose(rdensfp);
  }
    

  // writing IDL script

  char scriptfile[255];
  FILE* scriptfp;

  char* fmdhome;

  if (!(fmdhome = getenv("FMD"))) {
    fprintf(stderr, "environment variable FMD not defined\n");
    exit(-1);
  }

  snprintf(scriptfile, 255, "%s.rdens.script", slaterdetfile);
  if (!(scriptfp = fopen(scriptfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", scriptfile);
    exit(-1);
  }

  fprintf(scriptfp, 
	  "; plot densities of nuclei\n; written by calcdens\n\n");
  fprintf(scriptfp, "densxy=fltarr(%d,%d,%d)\n", 
	  npoints, npoints, 1);\
  fprintf(scriptfp, "rdens=fltarr(%d,3,%d)\n",
	  npoints, 1);

  fprintf(scriptfp, "!path = '%s/lib:' + !path\n", fmdhome);
  fprintf(scriptfp, "openr, unit, '%s', /get_lun\n", datafile);
  if (shape)
    fprintf(scriptfp, "readf, unit, densxy\n");
  fprintf(scriptfp, "readf, unit, rdens\n");
  fprintf(scriptfp, "free_lun, unit\n");
  
  fprintf(scriptfp, "\n.run multipost, densityplot, radialdensityplot\n");
  fprintf(scriptfp, "!p.thick = 1.5\n");
  fprintf(scriptfp, "pos = initpost('%s.rdens.eps', %d, 1, gapx=2.0, /color)\n",
	  slaterdetfile, 1+shape);
  fprintf(scriptfp, "!x.minor = 0\n!y.minor = 0\n");
  fprintf(scriptfp, "loadct, 3\n");

  fprintf(scriptfp, "!p.position = pos(0,*)\n");
  if (shape) {
    fprintf(scriptfp, "densityplot, densxy(*,*,%d), %d, %f, %f, %f, %f, $\n",
	    0, view, -xmax, xmax, -xmax, xmax);
    fprintf(scriptfp, "\t/cut, /coordinate, /cont, /dens, key = '%s'\n", 
	    label);

    fprintf(scriptfp, "!p.position = pos(1,*)\n");
  }
  fprintf(scriptfp, "!x.minor = 4\n!y.minor = 0\n");
  fprintf(scriptfp, "radialdensityplot, rdens(*,*,0), %f, $\n", rmax);
  fprintf(scriptfp, "\t/coordinate, key = \"%s\", /legend, /col\n",
	  label);
  fprintf(scriptfp, "exit\n");

  fclose(scriptfp);

  // calling IDL

  char call[255];

  snprintf(call, 255, "idl < %s", scriptfile);
  system(call);

  char epsfile[255];
  snprintf(epsfile, 255, "%s.rdens.eps", slaterdetfile);
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


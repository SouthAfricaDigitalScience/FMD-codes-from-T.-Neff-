/**

  \file calcdensitiespn.c

  calc density cuts in coordinate space for protons and neutrons
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
	    "\n   -C LOGO       show logo"
	    "\n   -P            portrait orientation"
	    "\n   -o            reorient"
	    "\n   -f            no frame and axes"
	    "\n   -x            no x-axis"
            "\n   -n            no labels, no annotations"
	    "\n   -l LABEL	use label"
	    "\n   -v {xy|yz|xz} select view"
	    "\n   -c RANGE      coordinate range\n",
	    argv[0]);
    exit(-1);
  }

  int debug=0, quiet=0;
  int portrait = 0;
  int orient=0;
  int view[3] = {1, 1, 1};
  int nplots[2] = {2, 3};
  int frame=1;
  int noxaxis=0;
  int plain=0;

  double xmax=5.5;

  int npoints=39+1;

  char* logo = NULL;
  char* label = NULL;

  char c;
  while ((c = getopt(argc, argv, "c:nNl:ogC:Pqv:fx")) != -1)
    switch (c) {
    case 'g':
      debug=1;
      break;
    case 'q':
      quiet=1;
      break;
    case 'P':
      portrait=1;
      break;
    case 'C':
      logo = optarg;
      break;
    case 'o':
      orient=1;
      break;
    case 'f':
      frame=0;
      break;
    case 'x':
      noxaxis=1;
      break;
    case 'n':
      plain=1;
      break;
    case 'l':
      label=optarg;
      break;
    case 'c':
      xmax = atof(optarg);
      break;
    case 'v':
      view[0]=0; view[1]=0; view[2]=0;
      nplots[1] = 1;
      if (!strcmp(optarg, "xy")) view[2]=1;
      if (!strcmp(optarg, "xz")) view[1]=1;
      if (!strcmp(optarg, "yz")) view[0]=1;
      break;
    }	

  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  if (readSlaterDetfromFile(&Q, slaterdetfile))
    exit(-1);

  if (orient) {
    SlaterDetAux X;
    initSlaterDetAux(&Q, &X);
    calcSlaterDetAux(&Q, &X);
    orientSlaterDet(&Q, &X);
  }

  if (!label)
    label = nucleusIDLformat(nucleusname(Q.A, Q.Z));
  
  double* dens = (double*) malloc(3*SQR(npoints)*sizeof(double));

  // write densities to data files

  char datafile[255];
  FILE *datafp;
  int v;

  snprintf(datafile, 255, "%s.denspn", slaterdetfile);;
  if (!(datafp = fopen(datafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", datafile);
    exit(-1);
  }

  // coordinate space densities
  for (v=0; v<3; v++)
    if (view[v]) {

      calcDensitiesCoordinate(&Q, v, npoints, xmax, dens);

      int i,j;

      // protons
      for (j=0; j<npoints; j++) {
	for (i=0; i<npoints; i++)
	  fprintf(datafp, "%f  ", 
		  dens[i+j*npoints+npoints*npoints]/(0.5*rho0));
	fprintf(datafp, "\n");
      }

      // neutrons
      for (j=0; j<npoints; j++) {
	for (i=0; i<npoints; i++)
	  fprintf(datafp, "%f  ", 
		  dens[i+j*npoints+2*npoints*npoints]/(0.5*rho0));
	fprintf(datafp, "\n");
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

  snprintf(scriptfile, 255, "%s.denspn.script", slaterdetfile);
  if (!(scriptfp = fopen(scriptfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", scriptfile);
    exit(-1);
  }

  fprintf(scriptfp, 
	  "; plot densities of nuclei\n; written by calcdensities\n\n");
  fprintf(scriptfp, "denspn=fltarr(%d,%d,%d)\n", 
	  npoints, npoints, nplots[0]*nplots[1]);

  fprintf(scriptfp, "!path = '%s/lib:' + !path\n", fmdhome);
  fprintf(scriptfp, "openr, unit, '%s', /get_lun\n", datafile);
  fprintf(scriptfp, "readf, unit, denspn\n");
  fprintf(scriptfp, "free_lun, unit\n");
  
  fprintf(scriptfp, "\n.run multipost, densityplot\n");
  if (portrait)
    fprintf(scriptfp, 
	  "pos = initpost('%s.denspn.eps',%d,%d, gapx=1.5, gapy=-0.01, /color)\n",
	  slaterdetfile, nplots[1], nplots[0]);
  else
    fprintf(scriptfp, 
	  "pos = initpost('%s.denspn.eps',%d,%d, gapx=-0.01, gapy=1.5, /color)\n",
	  slaterdetfile, nplots[0], nplots[1]);
  fprintf(scriptfp, "!x.minor = 4\n!y.minor = 4\n");

  int i = 0;
  // proton
  for (v=0; v<3; v++) {
    if (view[v]) {
      fprintf(scriptfp, "!p.position = pos(%d,*)\n", portrait ? 2*i : i);
      fprintf(scriptfp, "loadct, 14\n");
      fprintf(scriptfp, "densityplot, denspn(*,*,%d), %d, %f, %f, %f, %f, $\n",
	      2*i, v, -xmax, xmax, -xmax, xmax);
      fprintf(scriptfp, "\t/cut, /coordinate, /cont, /dens");
      if (portrait) fprintf(scriptfp, ", /noxaxes");
      if (noxaxis) fprintf(scriptfp, ", /noxaxes");
      if (!frame) fprintf(scriptfp, ", /noframe");
      if (plain) fprintf(scriptfp, ", /noannot");
      else fprintf(scriptfp, ", key = '%s - p'", label);
      fprintf(scriptfp, "\n");

      i++;
    }
  }

  i = 0;
  // neutron
  for (v=0; v<3; v++) {
    if (view[v]) {
      fprintf(scriptfp, "!p.position = pos(%d,*)\n", portrait ? 2*i+1 : i+nplots[1] );
      fprintf(scriptfp, "loadct, 1\n");
      fprintf(scriptfp, "densityplot, denspn(*,*,%d), %d, %f, %f, %f, %f, $\n",
	      2*i+1, v, -xmax, xmax, -xmax, xmax);
      fprintf(scriptfp, "\t/cut, /coordinate, /cont, /dens");
      if (!portrait) fprintf(scriptfp, ", /noyaxes");
      if (noxaxis) fprintf(scriptfp, ", /noxaxes");
      if (!frame) fprintf(scriptfp, ", /noframe");
      if (plain) fprintf(scriptfp, ", /noannot");
      else {
	fprintf(scriptfp, ", key = '%s - n'", label);
	if (logo) fprintf(scriptfp, ", logo = '%s'", logo);
      }
      fprintf(scriptfp, "\n");

      i++;
    }
  }

  fclose(scriptfp);

  // calling IDL

  char call[255];

  snprintf(call, 255, "idl < %s", scriptfile);
  system(call);

  char epsfile[255];
  snprintf(epsfile, 255, "%s.denspn.eps", slaterdetfile);

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


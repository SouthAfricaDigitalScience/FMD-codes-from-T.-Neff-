/**

  \file calcdensities3d.c

  calc densitions and generate DataExplorer plot


  (c) 2003 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "fmd/SlaterDet.h"
#include "fmd/SpatialOrientation.h"
#include "fmd/Densities3d.h"

#include "misc/physics.h"


#define CUBE(x)	(x)*(x)*(x)


#define INDIVIDUAL 0
#define CORNER 1
#define PLANE 2

int main(int argc, char* argv[])
{

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "Usage: %s [OPTIONS] state\n"
	    "\n   -g            debug"
	    "\n   -q            quiet	    "
	    "\n   -C            cut out corner"
	    "\n   -P            cut out plane"
	    "\n   -o A,B,G      rotate using Euler angles"
	    "\n   -l LABEL	plot label"
	    "\n   -n POINTS     number of grid points"
	    "\n   -c RANGE      coordinate range\n",
	    argv[0]);
    exit(-1);
  }

  int debug=0, quiet=0;

  int cut=CORNER;
  int orientation=0;
  double xmax=5.5;
  int npoints=49+1;
  double alpha=0.0, beta=0.0, gamma=0.0;

  char* label = NULL;

  char c;
  while ((c = getopt(argc, argv, "c:n:o:CPl:gq")) != -1)
    switch (c) {
    case 'g':
      debug=1;
      break;
    case 'q':
      quiet=1;
      break;
    case 'C':
      cut=CORNER;
      orientation=CORNER;
      break;
    case 'P':
      cut=PLANE;
      orientation=PLANE;
      break;
    case 'o':
      orientation=INDIVIDUAL;
      sscanf(optarg, "%lf,%lf,%lf", &alpha, &beta, &gamma);
      alpha *= M_PI/180; beta *= M_PI/180; gamma *= M_PI/180;
      break;
    case 'l':
      label=optarg;
      break;
    case 'c':
      xmax = atof(optarg);
      break;
    case 'n':
      npoints = atoi(optarg);
      break;
    }	

  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  readSlaterDetfromFile(&Q, slaterdetfile);

  if (!label)
    label = nucleusname(Q.A, Q.Z);
  
  double* dens = (double*) malloc(CUBE(npoints)*sizeof(double));

  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  // spatial orientation 
  if (orientation==INDIVIDUAL) {
    rotateSlaterDet(&Q, alpha, beta, gamma);
  } else if (orientation==CORNER) {
    calcSpatialOrientation(&Q, &X, &alpha, &beta, &gamma);
    rotateSlaterDet(&Q, alpha, beta, gamma);
  } else if (cut==PLANE) {
    calcSpatialOrientation(&Q, &X, &alpha, &beta, &gamma);
    rotateSlaterDet(&Q, alpha, beta, gamma);
    // rotate z to x axis
    rotateSlaterDet(&Q, 0, M_PI/2, 0);
  }

  calcSlaterDetAux(&Q, &X);

  // write densities to data files

  char datafile[255];
  FILE *datafp;

  // coordinate space densities

  snprintf(datafile, 255, "%s.cdens3d", slaterdetfile);;
  if (!(datafp = fopen(datafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", datafile);
    exit(-1);
  }
    
  fprintf(datafp, "# coordinate space densities of %s"
	  "\n# data explorer format"
	  "\n# x|y|z : %5.2f - %5.2f"
	  "\n# %d grid points in each dimension\n\n", 
	  slaterdetfile, -xmax, xmax, npoints);

  fprintf(datafp, "object 1 class array items %d data follows\n\n", 
	  CUBE(npoints));

  calcDensitiesCoordinate3d(&Q, &X, npoints, xmax, dens);

  int i,j,k;
  for (k=0; k<npoints; k++) {
    for (j=0; j<npoints; j++) {
      for (i=0; i<npoints; i++)
	fprintf(datafp, "%f    ", dens[k*npoints*npoints+j*npoints+i]/rho0);
      fprintf(datafp, "\n");
    }
    fprintf(datafp, "\n");
  }

  fprintf(datafp, "attribute \"dep\" string \"positions\""
	  "\nobject 2 class gridpositions counts %d %d %d"
	  "\norigin %5.2f %5.2f %5.2f"
	  "\ndelta %5.2f 0.0 0.0"
	  "\ndelta 0.0 %5.2f 0.0"
	  "\ndelta 0.0 0.0 %5.2f\n", 
	  npoints, npoints, npoints, 
	  -xmax, -xmax, -xmax, 
	  2*xmax/(npoints-1), 2*xmax/(npoints-1), 2*xmax/(npoints-1));

  fprintf(datafp, "object 3 class gridconnections counts %d %d %d"
	  "\nattribute \"element type\" string \"cubes\""
	  "\nattribute \"ref\" string \"positions\""
	  "\n"
	  "\nobject \"density\" class field"
	  "\ncomponent \"data\" 1 component \"positions\" 2 component \"connections\" 3\n",
	  npoints, npoints, npoints);
  fclose(datafp);

  // writing the DataExplorer script

  char scriptfile[255];
  FILE* scriptfp;

  char* fmdhome;

  if (!(fmdhome = getenv("FMD"))) {
    fprintf(stderr, "environment variable FMD not defined\n");
    exit(-1);
  }

  snprintf(scriptfile, 255, "%s.dens3d.script", slaterdetfile);
  if (!(scriptfp = fopen(scriptfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", scriptfile);
    exit(-1);
  }

  // todo: specify plotrange and cutplanes

  fprintf(scriptfp, "dataimport=\"%s.cdens3d\";\n", slaterdetfile);
  fprintf(scriptfp, "nucleusdens=\"%s.cdens3d\";\n", slaterdetfile);
  fprintf(scriptfp, "nucleusname=\"%s\";\n", label);
  fprintf(scriptfp, "plotformat=\"tiff\";\n");
  if (cut==CORNER)
    fprintf(scriptfp, "include \"%s/lib/Cut3D-corner-PlotsDX.net\";\n", fmdhome);
  else if (cut==PLANE)
    fprintf(scriptfp, "include \"%s/lib/Cut3D-plane-PlotsDX.net\";\n", fmdhome);
  fclose(scriptfp);

  // calling DataExplorer

  char call[255];

  snprintf(call, 255, "dx -processors 1 -script -file %s", scriptfile);
  system(call);

  // calling convert

  char tifffile[255]; char jpgfile[255]; char epsfile[255];
  snprintf(tifffile, 255, "%s.cdens3d.tiff", slaterdetfile);
  snprintf(jpgfile, 255, "%s.cdens3d.jpg", slaterdetfile);
  snprintf(epsfile, 255, "%s.cdens3d.eps", slaterdetfile);

  snprintf(call, 255, "convert %s %s; jpeg2ps %s > %s", 
	   tifffile, jpgfile, jpgfile, epsfile);
  system(call);

  // calling gv

  if (!quiet) {
    snprintf(call, 255, "gv %s &", epsfile);
    system(call);
  }

  // clean up

  if (!debug) {
    remove(datafile); 
    remove(scriptfile);
    remove(tifffile);
    remove(jpgfile);
  }

  return 0;
}


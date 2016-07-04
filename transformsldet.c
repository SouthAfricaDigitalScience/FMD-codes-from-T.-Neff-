/**

  \file transformsldet.c

  move, boost, rotate, invert SlaterDet

  
  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "fmd/Parameterization.h"
#include "fmd/ParameterizationFMD.h"
#include "fmd/SlaterDet.h"
#include "fmd/CenterofMass.h"
#include "fmd/SpatialOrientation.h"

#include "misc/utils.h"
#include "misc/physics.h"


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  int origin=0;
  int invert=0;
  int timerevert=0;
  int orientinertia=0;
  int orientmaxinertia=0;
  int orientangular=0;

  double d[3] = {0.0, 0.0, 0.0};
  double v[3] = {0.0, 0.0, 0.0};
  double alpha=0.0, beta=0.0, gamma=0.0;
  double scale=0.0;
  
  int c, i;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s [OPTIONS] Nuc transformedNuc\n"
	    "\n  -o                      move to origin in phase-space"
	    "\n  -p                      invert Nuc"
	    "\n  -T			 timerevert Nuc"
	    "\n  -s SCALE		 scale lengths by SCALE"
	    "\n  -t X,Y,Z                move Nuc by (X,Y,Z)"
	    "\n  -b VX,VY,VZ             boost Nuc by (VX,VY,VZ)"
            "\n  -r ALPHA,BETA,GAMMA     rotate Nuc by ALPHA, BETA, GAMMA (degrees)"
	    "\n  -I                      orient along inertia tensor"
	    "\n  -M                      orient along maximal inertia axis"
	    "\n  -J                      orient along angular momenta\n",
	    argv[0]);
    exit(-1);
  }

  /* manage command line options */

  while ((c = getopt(argc, argv, "opTs:t:b:r:IMJ")) != -1)
    switch (c) {
    case 'o':
      origin = 1;
      break;
    case 'p':
      invert = 1;
      break;
    case 'T':
      timerevert = 1;
      break;
    case 's':
      scale = atof(optarg);
      break;
    case 't':
      sscanf(optarg, "%lf,%lf,%lf", &d[0], &d[1], &d[2]);
      break;
    case 'b':
      sscanf(optarg, "%lf,%lf,%lf", &v[0], &v[1], &v[2]);
      break;
    case 'r':
      sscanf(optarg, "%lf,%lf,%lf", &alpha, &beta, &gamma);
      alpha *= M_PI/180;
      beta *= M_PI/180;
      gamma *= M_PI/180;
      break;
    case 'I':
      orientinertia = 1;
      break;
    case 'M':
      orientmaxinertia = 1;
      break;
    case 'J':
      orientangular = 1;
      break;
    }
     
 
  char* Qfile = argv[optind];
  char* transQfile = argv[optind+1];

  SlaterDet Q;
  SlaterDetAux X;
  readSlaterDetfromFile(&Q, Qfile);

  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  if (origin) {
    calcCMPosition(&Q, &X, d);
    calcCMVelocity(&Q, &X, v);  

    for (i=0; i<3; i++)	
      v[i] *= -1;
    for (i=0; i<3; i++)
      d[i] *= -1;
  }

  if (orientinertia) {
    double alpha0, beta0, gamma0;
    calcSpatialOrientation(&Q, &X, &alpha0, &beta0, &gamma0);
    alpha=-gamma0; beta=-beta0; gamma=-alpha0;
  }

  if (orientmaxinertia) {
    double alpha0, beta0, gamma0;
    calcOrientedOrientation(&Q, &X, &alpha0, &beta0, &gamma0);
    alpha=-gamma0; beta=-beta0; gamma=-alpha0;
  }

  if (orientangular) {
    double alpha0, beta0, gamma0;
    calcSpinOrientation(&Q, &X, &alpha0, &beta0, &gamma0);
    alpha=-gamma0; beta=-beta0; gamma=-alpha0;
  }

  if (scale > 0.0)
    scaleSlaterDet(&Q, scale);

  if (invert)
    invertSlaterDet(&Q);

  if (timerevert)
    timerevertSlaterDet(&Q);

  rotateSlaterDet(&Q, alpha, beta, gamma);
  boostSlaterDet(&Q, v);
  moveSlaterDet(&Q, d);

  Parameterization P = ParameterizationFMD;
  Para q;
  SlaterDetinitFMD(&Q, &q);

  // write transformed SlaterDet parameters
  
  FILE* outfp;
 
  fprintf(stderr, "... writing transformed SlaterDet to file %s\n", transQfile);
  if (!(outfp = fopen(transQfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", transQfile);
    exit(-1);
  }

  fprintinfo(outfp);

  fprintf(outfp, "\n# Parameterization\n");
  fprintf(outfp, "<Parameterization %s>\n", P.name);
  P.Parawrite(outfp, &q);

  fprintf(outfp, "\n# SlaterDet\n");
  writeSlaterDet(outfp, &Q);

  fclose(outfp);

  return 0;
}

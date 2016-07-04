/**

  \file calcquadrupole.c

  calculate quadrupole deformation parameters


  (c) 2010 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Quadrupole.h"

#include "misc/utils.h"


int main(int argc, char *argv[])
{
  createinfo(argc, argv);
  
  if (argc < 2) {
    fprintf(stderr, "\nusage: %s slaterdetfile\n",
	    argv[0]);
    exit(-1);
  }
  
  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  readSlaterDetfromFile(&Q, slaterdetfile);

  fprintinfo(stdout);

  SlaterDetAux X;
    
  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  double Beta[3], gamma[3];

  calcQuadrupole(&Q, &X, Beta, gamma);

  fprintf(stdout, "\nQuadrupole deformation:\n");
  fprintf(stdout, "\t nucleon: Beta = %8.3f fm^2, gamma = %8.3f\n",
          Beta[0], gamma[0]);
  fprintf(stdout, "\t proton:  Beta = %8.3f fm^2, gamma = %8.3f\n",
          Beta[1], gamma[1]);
  fprintf(stdout, "\t neutron: Beta = %8.3f fm^2, gamma = %8.3f\n",
          Beta[2], gamma[2]);

  calcQuadrupolePrime(&Q, &X, Beta);

  fprintf(stdout, "\nQuadrupole deformation:\n");
  fprintf(stdout, "\t nucleon: Beta' = %8.3f\n",
          Beta[0]);
  fprintf(stdout, "\t proton:  Beta' = %8.3f\n",
          Beta[1]);
  fprintf(stdout, "\t neutron: Beta' = %8.3f\n",
          Beta[2]);

  return 0;
}

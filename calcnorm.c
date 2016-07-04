/**

  \file calcnorm.c

  calculate norm of Slater determinant


  (c) 2004 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/CenterofMass.h"

#include "misc/utils.h"
#include "misc/physics.h"


int main(int argc, char *argv[])
{
  createinfo(argc, argv);
  
  if (argc < 2) {
    fprintf(stderr, "\nusage: %s slaterdetfile\n", argv[0]);
    exit(-1);
  }
  
  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  readSlaterDetfromFile(&Q, slaterdetfile);

  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);


  calcSlaterDetAuxod(&Q, &Q, &X);

  double norm = sqrt(creal(X.ovlap));

  double tcm;
  calcTCM(&Q, &X, &tcm);
  double a = Q.A*0.75/(tcm*(mproton*Q.Z+mneutron*Q.N));

  fprintinfo(stdout);

  fprintf(stdout, "norm: %15.8e\n", norm);
  fprintf(stdout, "a   : %15.8e\n", a);

  return 0;
}

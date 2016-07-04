/**

  \file calcorientation.c

  calculate expectation value of angular momentum J
  and inertia tensor


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/AngularMomenta.h"
#include "fmd/SpatialOrientation.h"

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

  calcSlaterDetAux(&Q, &X);

  double L2, S2, J2;
  double L[3], S[3], J[3];
  double T[3][3];

  calcAngularMomenta(&Q, &X, &L2, &S2, &J2);
  calcL(&Q, &X, L);
  calcS(&Q, &X, S);
  calcJ(&Q, &X, J);
  calcInertiaTensor(&Q, &X, T);

  fprintinfo(stdout);

  fprintf(stdout, "J2: %8.5f\n\n", J2);

  fprintf(stdout, "L: (%8.5f, %8.5f, %8.5f)\n", L[0], L[1], L[2]);
  fprintf(stdout, "S: (%8.5f, %8.5f, %8.5f)\n", S[0], S[1], S[2]);
  fprintf(stdout, "J: (%8.5f, %8.5f, %8.5f)\n\n", J[0], J[1], J[2]);

  fprintf(stdout, "T:  ((%8.5f, %8.5f, %8.5f),\n", T[0][0], T[1][0], T[2][0]);
  fprintf(stdout, "     (%8.5f, %8.5f, %8.5f),\n", T[0][1], T[1][1], T[2][1]);
  fprintf(stdout, "     (%8.5f, %8.5f, %8.5f))\n\n", T[0][2], T[1][2], T[2][2]);

  return 0;
}

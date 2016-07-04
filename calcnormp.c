/**

  \file calcnormp.c

  calculate norm of Parity Projected Slater determinant


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

  SlaterDet Q, Qp;
  readSlaterDetfromFile(&Q, slaterdetfile);
  initSlaterDet(&Q, &Qp);

  copySlaterDet(&Q, &Qp);
  invertSlaterDet(&Qp);

  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);

  calcSlaterDetAuxod(&Q, &Q, &X);
  double ovld = creal(X.ovlap);

  calcSlaterDetAuxod(&Q, &Qp, &X);
  double ovle = creal(X.ovlap);

  double norm = sqrt(ovld);
  double normp = sqrt(0.5*(ovld+ovle));
  double normm = sqrt(0.5*(ovld-ovle));

  complex double tcmd;
  calcSlaterDetAuxod(&Q, &Q, &X);
  calcTCMod(&Q, &Q, &X, &tcmd);

  complex double tcme;
  calcSlaterDetAuxod(&Q, &Qp, &X);
  calcTCMod(&Q, &Qp, &X, &tcme);

  double tcm = tcmd/ovld;
  double tcmp = (tcmd+tcme)/(ovld+ovle);
  double tcmm = (tcmd-tcme)/(ovld-ovle);

  double a = Q.A*0.75/(tcm*(mproton*Q.Z+mneutron*Q.N));
  double ap = Q.A*0.75/(tcmp*(mproton*Q.Z+mneutron*Q.N));
  double am = Q.A*0.75/(tcmm*(mproton*Q.Z+mneutron*Q.N));

  fprintinfo(stdout);

  fprintf(stdout, "norm: %15.8e\n", norm);
  fprintf(stdout, "norm pi=+1: %15.8e\n", normp);
  fprintf(stdout, "norm pi=-1: %15.8e\n", normm);

  fprintf(stdout, "effective width a: %15.8e\n", a);
  fprintf(stdout, "effective width a pi=+1: %15.8e\n", ap);
  fprintf(stdout, "effective width a pi=-1: %15.8e\n", am);

  return 0;
}

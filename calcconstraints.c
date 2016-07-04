/**

  \file calcconstraints.c

  calculate values for constraint observables


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"

#include "fmd/Constraint.h"
#include "fmd/ConstraintR2.h"
#include "fmd/ConstraintDipole.h"
#include "fmd/ConstraintQuadrupole.h"
#include "fmd/ConstraintOctupole.h"
#include "fmd/ConstraintNOsci.h"
#include "fmd/ConstraintS2.h"
#include "fmd/ConstraintT2.h"


#include "misc/utils.h"


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

  double r2, pr2, nr2, ed2, q2, pq2, nq2, o2, po2, no2, nosci, pnosci, nnosci;
  double s2, ps2, ns2, t2, pt2, nt2;

  calcConstraintR2(&Q, &X, &r2);
  calcConstraintER2(&Q, &X, &pr2);
  calcConstraintNR2(&Q, &X, &nr2);

  calcConstraintED2(&Q, &X, &ed2);

  calcConstraintQ2(&Q, &X, &q2);
  calcConstraintEQ2(&Q, &X, &pq2);
  calcConstraintNQ2(&Q, &X, &nq2);

  calcConstraintO2(&Q, &X, &o2);
  calcConstraintEO2(&Q, &X, &po2);
  calcConstraintNO2(&Q, &X, &no2);

  calcConstraintNOsci(&Q, &X, &nosci);
  calcConstraintPNOsci(&Q, &X, &pnosci);
  calcConstraintNNOsci(&Q, &X, &nnosci);

  calcConstraintS2(&Q, &X, &s2);
  calcConstraintPS2(&Q, &X, &ps2);
  calcConstraintNS2(&Q, &X, &ns2);

  calcConstraintT2(&Q, &X, &t2);
  calcConstraintPT2(&Q, &X, &pt2);
  calcConstraintNT2(&Q, &X, &nt2);

  fprintf(stdout, "\nR:     %8.3f    PR:    %8.3f    NR:    %8.3f\n", sqrt(r2), sqrt(pr2), sqrt(nr2));
  fprintf(stdout, "                   ED2:   %8.3f\n", ed2);
  fprintf(stdout, "Q2:    %8.3f    PQ2:   %8.3f    NQ2:   %8.3f\n", q2, pq2, nq2);
  fprintf(stdout, "O2:    %8.3f    PO2:   %8.3f    NO2:   %8.3f\n", o2, po2, no2);
  fprintf(stdout, "NOsci: %8.3f    PNOsci:%8.3f    NNOsci:%8.3f\n", sqrt(nosci), sqrt(pnosci), sqrt(nnosci));
  fprintf(stdout, "S2:    %8.3f    PS2:   %8.3f    NS2:   %8.3f\n", s2, ps2, ns2);
  fprintf(stdout, "T2:    %8.3f    PT2:   %8.3f    NT2:   %8.3f\n", t2, pt2, nt2);


  return 0;
}

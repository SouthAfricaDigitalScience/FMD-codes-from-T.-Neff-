/**

  \file calctwobodydensitiesmom

  calculate two-body densities in momentum space


  (c) 2011 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/TwoBodyDensity.h"

#include "misc/utils.h"
#include "misc/physics.h"




int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s [OPTIONS] mcstate"
	    "\n     -j J            2j of state"
	    "\n     -p +1|-1        parity of state"
	    "\n     -a INDEX        index of state"
	    "\n     -q QMAX         calculate to max separation"
	    "\n     -n NPOINTS      number of points\n",
	    argv[0]);
    exit(-1);
  }

  int hermit = 0;

  // grid in coordinate space

  int npoints = 51;
  double qmax = 10.0;

  // project to this state
  int j=-1, pi=-1, alpha=0;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "j:p:a:q:n:")) != -1)
    switch (c) {
    case 'j':
      j = atoi(optarg);
      break;
    case 'p':
      pi = (atoi(optarg) == -1 ? 1 : 0);
      break;
    case 'a':
      alpha = atoi(optarg);
      break;
    case 'q':
      qmax = atof(optarg);
      break;
    case 'n':
      npoints = atoi(optarg);
      break;
    }

  char* mcstatefile = argv[optind];
  char** mbfile;

  // open multiconfigfile
  Projection P;
  SlaterDet* Q;
  Symmetry* S;
  Eigenstates E;
  int n;

  if (readMulticonfigfile(mcstatefile, &mbfile, &P, &Q, &S, &E, &n)) {
    fprintf(stderr, "couldn't read %s\n", mcstatefile);
    exit(-1);
  }

  TBDensQPara TBDP = {
    qmax : qmax,
    npoints : npoints
  };

  int a,b;

  // J or Pi unset
  if (j == -1) j = (Q[0].A % 2 ? 1 : 0);
  if (pi == -1) pi = (Q[0].A % 2 ? 1 : 0); 

  fprintf(stderr, "J: %d, Pi: %c1, alpha: %d\n", j, pi ? '-' : '+', alpha);

  // initialize operator

  initOpTwoBodyDensityQ(&TBDP);

  void* tbdensme[n*n];

  // read or calculate matrix elements

  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {
      tbdensme[a+b*n] = initprojectedMBME(&P, &OpTwoBodyDensityQ);

      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P,
				    &OpTwoBodyDensityQ, S[a], S[b],
				    tbdensme[a+b*n])) {
	calcprojectedMBME(&P, &OpTwoBodyDensityQ, &Q[a], &Q[b],
			  S[a], S[b], tbdensme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpTwoBodyDensityQ, S[a], S[b], 
				 tbdensme[a+b*n]);
      }
    }

  if (hermit) {
    hermitizeprojectedMBME(&P, &OpTwoBodyDensityQ, tbdensme, n);
  }

  fprintf(stderr, "calculate expectation values\n");

  // expectation values

  void* tbdensexp;

  tbdensexp = initprojectedVectornull(&P, &OpTwoBodyDensityQ, n);
  calcexpectprojectedMBMEipj(&P, &OpTwoBodyDensityQ, tbdensme, 
				 S, &E, j, pi, alpha, 
				 tbdensexp);
  // output

  char outfile[255];

  snprintf(outfile, 255, "%s-%s.%d.tbdensmom", 
	   stripstr(mcstatefile, ".states"), AngmomtoStr(j, pi), alpha);

  // write densities to file

  char datafile[255];
  FILE *datafp;

  snprintf(datafile, 255, "%s--%05.2f-%d.data", 
	   outfile,
	   qmax, npoints);

  if (!(datafp = fopen(datafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", datafile);
    exit(-1);
  }

  writeTBDensQ(datafp, &P, &TBDP, 
               j, pi, alpha, tbdensexp, &E);

  fclose(datafp);


  return 0;
}

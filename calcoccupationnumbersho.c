/**

  \file calcoccupationsnumbersho.c

  calculate occupation numbers in harmonic oscillator basis


  (c) 2006 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "fmd/SlaterDet.h"
#include "fmd/CenterofMass.h"
#include "fmd/HOBasis.h"
#include "fmd/DensityMatrixHO.h"

#include "misc/physics.h"
#include "misc/utils.h"

#define SQR(x)	(x)*(x)
#define ABS(x)  ((x)<0 ? -(x) : (x))

int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  // enough arguments ?

  if (argc < 2) {
    fprintf(stderr, "Usage: %s [OPTIONS] state\n"
	    "\n   -n NMAX	oscillator cut off"
	    "\n   -o OMEGA	Oscillator constant [MeV]\n",
	    argv[0]);
    exit(-1);
  }

  int nmax = HONMAX;
  double omega = 0.0;

  char c;
  while ((c = getopt(argc, argv, "n:o:")) != -1)
    switch (c) {
    case 'n':
      nmax = atoi(optarg);
      break;
    case 'o':
      omega = atof(optarg)/hbc;
      break;
    }


  char* slaterdetfile = argv[optind];

  SlaterDet Q;
  readSlaterDetfromFile(&Q, slaterdetfile);


  // single-particle overlap matrix
  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  // Tcm
  double tcm;
  calcTCM(&Q, &X, &tcm);

  // derive oscillator constant from center of mass motion
  double omegacm = 4.0/3.0*tcm;

  if (omega == 0.0)
    omega = omegacm;

  
  // density matrix

  DensityMatrixHOPar par = {
    nmax : nmax,
    omega : omega,
    dim : dimHOBasis(nmax)
  };

  int dim = par.dim;

  complex double* rho = malloc(SQR(dim)*sizeof(complex double));

  initHOBasis(nmax);

  
  calcDensityMatrixHO(&par, &Q, &X, rho);

  fprintinfo(stdout);
  
  // occupation numbers

  fprintf(stdout, "\nusing oscillator constant: hbar Omega = %8.3f MeV\n",
          omega*hbc);

  fprintf(stdout, "\nOccupation numbers: \n");

  int i;
  int ixi, xi, N, n, l, twoj, twom, orbit, idx;
  double nshell, norbit[HONMAX+1];

  for (N=0; N<=nmax; N++) {

    fprintf(stdout, "\n   N=%d\t\t", N);
    for (n=N/2; n>=0; n--) {
      l=N-2*n;
      for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2)
	fprintf(stdout, "%s\t", labelHOorbit(n, l, twoj));
    }
    fprintf(stdout, "\n");

    for (ixi=0; ixi<=1; ixi++) {
      xi = ixi ? -1 : +1;
      fprintf(stdout, "%c  ", ixi ? 'n' : 'p');

      nshell = 0.0;
      orbit = 0;
      for (n=N/2; n>=0; n--) {
	l=N-2*n;
	for (twoj=ABS(2*l-1); twoj<=2*l+1; twoj=twoj+2) {
	  norbit[orbit] = 0.0;
	  for (twom=-twoj; twom<=twoj; twom=twom+2) {
	    idx = HOxinljmidx(xi, n, l, twoj, twom);
	    norbit[orbit] += creal(rho[idx+idx*dim]);
	  }
	  nshell += norbit[orbit];
	  orbit++;
	}
      }
      fprintf(stdout, "%6.3f\t", nshell);
      for (orbit=0; orbit<=N; orbit++)
	fprintf(stdout, "%6.3f\t", norbit[orbit]);
      fprintf(stdout, "\n");
    }
  }
  
  double ntot = 0.0;
  for (i=0; i<dim; i++)
    ntot += creal(rho[i+i*dim]);

  fprintf(stdout, "\nTotal Occupation: %6.3f\n", ntot);
  
}


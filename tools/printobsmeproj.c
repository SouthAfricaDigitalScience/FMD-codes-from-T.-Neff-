/**

  \file printobsmeproj.c

  print matrix elements with projected Slater determinants
  matrix elements must have been calculated with no explicit cm projection


  (c) 2005, 2006 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/CenterofMass.h"
#include "fmd/Projection.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"
#include "numerics/zcw.h"
#include "numerics/cmat.h"


void extractmatrices(const char* nucsfile,
		     const Projection* P,
		     const Symmetry* S,
		     const Interaction* Int,
		     const void *obs,
		     const double* norm,
		     const double cmfactor,
		     int n, int diagonal,
		     int j, int p)
{
  int odd=P->odd;
  int jmax=P->jmax;
  complex double H[n*n];
  complex double VC[n*n];
  complex double T[n*n];
  complex double R2[n*n];
  complex double N[n*n];

  complex double phase[n];
  
  const Observablesod ***obsme = obs;

  int a, b;
  int m, k;
  int nj, ipj;

  nj = n*(j+1);
  ipj = idxpij(jmax,p,j);

  // calculate phase relation with respect to first Slater det
  a=0;
  for (b=0; b<n; b++)
    for (m=-j; m<=j; m += 2)
      for (k=-j; k<=j; k+= 2)
	if (SymmetryAllowed(S[a], p, j, m) && SymmetryAllowed(S[b], p, j, k))
	  phase[b] = cexp(I*carg(obsme[a+b*n][ipj][idxjmk(j,m,k)].n));


  char outfilename[255];
  FILE* outfp;

  // Overlap matrix

  for (b=0; b<n; b++) {
    for (a=diagonal ? b : 0; a<n; a += diagonal ? n : 1)
      for (m=-j; m<=j; m += 2)
	for (k=-j; k<=j; k+= 2)
	  if (SymmetryAllowed(S[a], p, j, m) && SymmetryAllowed(S[b], p, j, k)) {
	    N[a+b*n] = obsme[a+b*n][ipj][idxjmk(j,m,k)].n/cmfactor*
	      norm[a]*norm[b]/(conj(phase[a])*phase[b]);
	  }
  }			

  snprintf(outfilename, 255, "%s.N%s-%d%c.bin", 
	   filepart(nucsfile), diagonal ? "matrixd" : "matrix", 
	   j, p ? '-' : '+');
  outfp = fopen(outfilename, "w");
  fwritecmatbin(outfp, n, N);
  fclose(outfp);

  // Hamiltonian matrix

  for (b=0; b<n; b++) {
    for (a=diagonal ? b : 0; a<n; a += diagonal ? n : 1)
      for (m=-j; m<=j; m += 2)
	for (k=-j; k<=j; k+= 2)
	  if (SymmetryAllowed(S[a], p, j, m) && SymmetryAllowed(S[b], p, j, k)) {
	    H[a+b*n] = obsme[a+b*n][ipj][idxjmk(j,m,k)].h/cmfactor*
	      norm[a]*norm[b]/(conj(phase[a])*phase[b]);
	  }
  }			

  snprintf(outfilename, 255, "%s.H%s-%d%c.bin", 
	   filepart(nucsfile), diagonal ? "matrixd" : "matrix", 
	   j, p ? '-' : '+');
  outfp = fopen(outfilename, "w");
  fwritecmatbin(outfp, n, H);
  fclose(outfp);

  // T matrix

  for (b=0; b<n; b++) {
    for (a=diagonal ? b : 0; a<n; a += diagonal ? n : 1)
      for (m=-j; m<=j; m += 2)
	for (k=-j; k<=j; k+= 2)
	  if (SymmetryAllowed(S[a], p, j, m) && SymmetryAllowed(S[b], p, j, k)) {
	    T[a+b*n] = obsme[a+b*n][ipj][idxjmk(j,m,k)].t/cmfactor*
	      norm[a]*norm[b]/(conj(phase[a])*phase[b]);
	  }
  }			

  snprintf(outfilename, 255, "%s.T%s-%d%c.bin", 
	   filepart(nucsfile), diagonal ? "matrixd" : "matrix", 
	   j, p ? '-' : '+');
  outfp = fopen(outfilename, "w");
  fwritecmatbin(outfp, n, T);
  fclose(outfp);

  // Coulomb matrix
  // assume Coulomb is latest potential component

  for (b=0; b<n; b++) {
    for (a=diagonal ? b : 0; a<n; a += diagonal ? n : 1)
      for (m=-j; m<=j; m += 2)
	for (k=-j; k<=j; k+= 2)
	  if (SymmetryAllowed(S[a], p, j, m) && SymmetryAllowed(S[b], p, j, k)) {
	    VC[a+b*n] = obsme[a+b*n][ipj][idxjmk(j,m,k)].v[Int->n-1]/cmfactor*
	      norm[a]*norm[b]/(conj(phase[a])*phase[b]);
	  }
  }			

  snprintf(outfilename, 255, "%s.VC%s-%d%c.bin", 
	   filepart(nucsfile), diagonal ? "matrixd" : "matrix", 
	   j, p ? '-' : '+');
  outfp = fopen(outfilename, "w");
  fwritecmatbin(outfp, n, VC);
  fclose(outfp);


  // R2 matrix

  for (b=0; b<n; b++) {
    for (a=diagonal ? b : 0; a<n; a += diagonal ? n : 1)
      for (m=-j; m<=j; m += 2)
	for (k=-j; k<=j; k+= 2)
	  if (SymmetryAllowed(S[a], p, j, m) && SymmetryAllowed(S[b], p, j, k)) {
	    R2[a+b*n] = obsme[a+b*n][ipj][idxjmk(j,m,k)].r2m/cmfactor*
	      norm[a]*norm[b]/(conj(phase[a])*phase[b]);
	  }
  }			

  snprintf(outfilename, 255, "%s.R2%s-%d%c.bin", 
	   filepart(nucsfile), diagonal ? "matrixd" : "matrix", 
	   j, p ? '-' : '+');
  outfp = fopen(outfilename, "w");
  fwritecmatbin(outfp, n, R2);
  fclose(outfp);

}
  
  
  
#define MAXSTATES 100

int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA INTERACTION NUCSFILE"
	    "\n   -h             hermitize matrix elements"
	    "\n   -d             diagonal matrix elements only\n",
	    argv[0]);
    exit(-1);
  }

  int hermit=0;
  int diagonal=0;
  int odd;


  char c;

  /* manage command-line options */

  while ((c = getopt(argc, argv, "dh")) != -1)
    switch (c) {
    case 'd':
      diagonal=1;
      break;
    case 'h':
      hermit=1;
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* nucsfile = argv[optind+2];

  char* mbfile[MAXSTATES];
  int n;

  if (readstringsfromfile(nucsfile, &n, mbfile))
    return -1;

  SlaterDet Q[n];
  Symmetry S[n];

  int i;
  for (i=0; i<n; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      exit(-1);;
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    exit(-1);
  Int.cm = 1;

  // odd numer of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  // check that no cm-projection was used
  if (P.cm != CMNone) {
    fprintf(stderr, "You have to use cm-none! for Projection\n");
    exit(-1);
  }
    
  initOpObservables(&Int);

  int a,b; 

  // calculate norms of Slater determinants
  SlaterDetAux X;
  double norm[n];

  initSlaterDetAux(&Q[0], &X);

  for (i=0; i<n; i++) {
    calcSlaterDetAuxod(&Q[i], &Q[i], &X);
    norm[i] = sqrt(creal(X.ovlap));
  }

  /* 
  // calculate cm factor
  double tcm[n];
  for (i=0; i<n; i++) {
    calcSlaterDetAux(&Q[i], &X);
    calcTCM(&Q[i], &X, &tcm[i]);
  }

  double meantcm = 0.0;
  for (i=0; i<n; i++)
    meantcm += tcm[i]/n;

  double alpha = 0.25/0.75*(meantcm*(mproton*Q[0].Z+mneutron*Q[0].N));
  double cmfactor = 0.125*pow(M_PI*alpha,-1.5);
  */

  // initialize space for matrix elements
  Observablesod** obsme[n*n];
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++)
      obsme[a+b*n] = initprojectedMBME(&P, &OpObservables);

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=diagonal ? b : 0; a<n; a += diagonal ? n : 1)
	if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, &OpObservables,
				    S[a], S[b], obsme[a+b*n])) {
	fprintf(stderr, "Matrix elements between %s and %s missing\n",
		mbfile[a], mbfile[b]);
	exit(-1);
      }

  if (hermit)
    hermitizeprojectedMBME(&P, &OpObservables, obsme, n);

  int pi, j;
  for (pi=0; pi<=1; pi++) 
    for (j=odd; j<P.jmax; j=j+2)
      extractmatrices(nucsfile, &P, S, &Int, obsme, norm, 1.0, n, diagonal, j, pi);

}

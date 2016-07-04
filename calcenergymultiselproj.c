/**

  \file calcenergymultiselproj.c

  multiconfiguration calculations with projected Slater determinants
 
  use subsets of Slater determinants for multiconfiguration calculations

  (c) 2004 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"


#define MAXSTATES 100


// generatate (n over m) ordered permutations in i[m]  

int nextround(int i[], int n, int m, int k)
{
  int j;

  if (k == -1) {
    return 0;
  }
  else if (i[k] == n-(m-k))
    return nextround(i, n, m, k-1);
  else {
    i[k]++;
    for (j=1; j<=m-k-1; j++)
      i[k+j] = i[k]+j;
    return 1;
  }
}

int nextpermutation(int i[], int n, int m)
{
  if (i[m-1] == n-1)
    return nextround(i, n, m, m-2);
  else {
    i[m-1]++;
    return 1;
  }
}

// is k in i[n] ?
int inidx(int i[], int n, int k)
{
  int j;
  for (j=0; j<n; j++)
    if (i[j] == k) return 1;
  
  return 0;
}

void intoidx(int idx[], int idxfix[], int iidxsel[], int ridxsel[], 
	    int n, int m, int m0)
{
  int k;

  for (k=0; k<m0; k++)
    idx[k] = idxfix[k];

  for (k=0; k<m; k++)
    idx[m0+k] = ridxsel[iidxsel[k]];
}



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPAR INTERACTION M NUCSFILE" 
	    "\n   -h             hermitize matrix elements"
	    "\n   -A             show all eigenstates"
	    "\n   -f FIXNUCSFILE select these configs in addition"
	    "\n                  has to be a subset of NUCSFILE"
	    "\n   -t THRESH      threshold for SVD\n", argv[0]);
    exit(-1);
  }

  int all=0;
  int hermit=0;
  int odd;
  int threshkmix=0.01;
  double minnormkmix=0.001;
  int threshmulti=0.0000001;
  double minnorm=0.001;

  char* fixmbfiles=NULL;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "hAt:f:")) != -1)
    switch (c) {
    case 'A':
      all=1;
      break;
    case 'h':
      hermit=1;
      break;
    case 't':
      threshmulti = atof(optarg);
      break;
    case 'f':
      fixmbfiles = optarg;
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  int m1 = atoi(argv[optind+2]);
  char* nucsfile = argv[optind+3];

  char* mbfile[MAXSTATES];
  int n;

  if (readstringsfromfile(nucsfile, &n, mbfile)) {
    fprintf(stderr, "couldn't open %s\n", nucsfile);
    exit(-1);
  }

  SlaterDet Q[n]; 
  Symmetry S[n];

  int i;
  for (i=0; i<n; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i])) {
      fprintf(stderr, "couldn't read from %s\n", mbfile[i]);
      exit(-1);
    }
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile)) {
    fprintf(stderr, "couldn't read Interaction from %s\n", interactionfile);
    exit(-1);
  }
  Int.cm = 1; 


  char* fixmbfile[MAXSTATES];
  int m0 = 0;

  // if specified read set of fixed configs
  if (fixmbfiles)
    if (readstringsfromfile(fixmbfiles, &m0, fixmbfile)) {
      fprintf(stderr, "couldn't read %s\n", fixmbfiles);
      exit(-1);
    }

  // find indices of fixed files
  
  int idxfix[m0];

  int k=0;
  for (i=0; i<m0; i++)
    for (k=0; k<n; k++)
      if (!strcmp(fixmbfile[i], mbfile[k])) {
	idxfix[i] = k;
	break;
      }

  // reverse indices for not fixed files

  int ridxsel[n-m0];
  int r;

  for (i=0; i<n-m0; i++) {
    r = (i == 0 ? 0 : ridxsel[i-1]+1);
    while (inidx(idxfix, m0, r))
      r++;
    ridxsel[i] = r;
  }

  int iidxsel[m1];

  for (i=0; i<m1; i++)
    iidxsel[i] = i;

  // odd number of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  initOpObservables(&Int);

  int a,b; 

  // initialize space for matrix elements
  Observablesod** obsme[n*n];
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++)
      obsme[a+b*n] = initprojectedMBME(&P, &OpObservables);

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++)
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P, &OpObservables, 
				    S[a], S[b], obsme[a+b*n])) {
	calcprojectedMBME(&P, &OpObservables, 
			  &Q[a], &Q[b], S[a], S[b], obsme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], 
				 &P, &OpObservables, S[a], S[b], obsme[a+b*n]);
      }

  // scale matrix elements
  if (Int.mescaling) {
    fprintf(stderr, "... scaling matrix elements\n");
    for (b=0; b<n; b++)
      for (a=0; a<n; a++)
        scaleprojectedObservablesMBME(&P, &Int, obsme[a+b*n]);
  }



  if (hermit)
    hermitizeprojectedMBME(&P, &OpObservables, obsme, n);

  // read or calculate the Eigenstates
      
  double minnormp=0.001; int allp=0;
  Eigenstates Ep[n];
  Observablesod** obsp = initprojectedVector(&P, &OpObservables, 1);
  
  for (i=0; i<n; i++) {
    if (readEigenstatesfromFile(mbfile[i], &P, &Ep[i], 1)) {
      fprintf(stderr, "... calculating Eigenstates for %s\n", mbfile[i]);
      initEigenstates(&P, &Ep[i], 1);
      calcEigenstates(&P, &Int, &obsme[i+i*n], &Ep[i], threshkmix);
      calcexpectprojectedMBME(&P, &OpObservables, &obsme[i+i*n], &S[i], &Ep[i], obsp);
      sortEigenstates(&P, &Int, obsp, &Ep[i], minnormkmix, allp);
    }
  }

  // output
  FILE* outfp = stdout;

  fprintinfo(outfp);
  fprintProjectinfo(outfp, &P);

  fprintf(outfp, "\n# calculating subsets of dimension %d+%d\n", m0, m1);

  // Now we have the eigenstates and matrixelements for n Slater determinants

  // next step is to select subsets of m Slater determinants

  int m = m0+m1;
  int idx[m];

  Symmetry Ssel[m];
  Eigenstates Eselp[m] ;
  Observablesod** obsmesel[m*m];
  
  Eigenstates multiE;
  Amplitudes multiA;
  initEigenstates(&P, &multiE, m);
  initAmplitudes(&P, &multiA, m);

  Observablesod** obs = initprojectedVector(&P, &OpObservables, m);

  for (i=0; i<m; i++)
    idx[i] = i;

  do {
    intoidx(idx, idxfix, iidxsel, ridxsel, n, m1, m0);

    // Eselp and obsmesel are pointers to the selected
    // Eigenstates and matrixelements

    for (a=0; a<m; a++) {
      Ssel[a] = S[idx[a]];
      Eselp[a] = Ep[idx[a]];
    }
		      
    for (b=0; b<m; b++)
      for (a=0; a<m; a++)
	obsmesel[a+b*m] = obsme[idx[a]+idx[b]*n];

    // solve the eigenvalue problem
    calcMultiEigenstates(&P, &Int, obsmesel, &Eselp, &multiE, &multiA, threshmulti);

    // calculate expectation values
    calcexpectprojectedMBME(&P, &OpObservables, obsmesel, Ssel, &multiE, obs);

    // sort Eigenstates
    sortEigenstates(&P, &Int, obs, &multiE, minnorm, all);

    // 
    fprintf(stdout, "\n###  ");
    for (a=0; a<m; a++)
      fprintf(stdout, "%s ", mbfile[idx[a]]);
    fprintf(stdout, "\n");

    showprojectedObservables(outfp, &P, &Int, &Q[0], obs, &multiE, &multiA, "");

  } while (nextpermutation(iidxsel, n-m0, m1));

  return 0;
}



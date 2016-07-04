/**

  \file calcbasisovlap.c

  calculate overlap of Many-Body eigenstates with 
  Basis spanned by a set of projected Many-Body states

  uses full set of basis states (all J,M,K projections)
  probably better to use only proper eigenstates


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"

#include "numerics/cmat.h"


#define MAXBASISSTATES 500


inline int SQR(int i) { return i*i; }

void writeprojectedOvlaps(FILE* fp,
			  const Projection* P,
			  const complex double** ovlap,
			  const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i;

  char prefix[8];

  int* idx;
  int ngood;
  complex double *norm;
  complex double *H, *Ovlap;
  

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

	if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
	else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

	idx = E->index[idxpij(jmax,p,j)];
	norm = E->norm[idxpij(jmax,p,j)];
	H = E->v[idxpij(jmax,p,j)];

	Ovlap = ovlap[idxpij(jmax,p,j)];

	fprintf(fp, "\n%s N     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
	fprintf(fp, "\n%s H     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));
	fprintf(fp, "\n%s Ovlap = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", creal(Ovlap[idx[i]]/norm[idx[i]]));
	fprintf(fp, "\n");	  
      }
    }
}


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 2) {
    fprintf(stderr, "\nusage: %s mcstate basistates\n", argv[0]);
    exit(-1);
  }

  char* multiconfigfile = argv[optind];
  char** mbfile;

  // open multiconfigfile
  Projection P;
  SlaterDet* Q;
  Symmetry* S;
  Eigenstates E;
  int n;

  if (readMulticonfigfile(multiconfigfile, &mbfile, &P, &Q, &S, &E, &n))
    exit(-1);

  char* basisfile = argv[optind+1];
  char* mbbfile[MAXBASISSTATES];
  int nb;

  if (readstringsfromfile(basisfile, &nb, mbbfile))
    exit(-1);

  SlaterDet Qb[nb]; 
  Symmetry Sb[nb];

  int i;
  for (i=0; i<nb; i++) {
    extractSymmetryfromString(&mbbfile[i], &Sb[i]);
    if (readSlaterDetfromFile(&Qb[i], mbbfile[i]))
      exit(-1);
  }


  // calculate overlaps for basis states
 
  complex double** ovlbbme[nb*nb];

  int a,b;
  for (b=0; b<nb; b++)
    for (a=0; a<nb; a++)
      ovlbbme[a+b*nb] = initprojectedMBME(&P, &OpOvlap);

  for (b=0; b<nb; b++)
    for (a=0; a<nb; a++) {
      if (readprojectedMBMEfromFile(mbbfile[a], mbbfile[b], &P, 
				    &OpOvlap, Sb[a], Sb[b], ovlbbme[a+b*nb])) {
	calcprojectedMBME(&P, &OpOvlap, &Qb[a], &Qb[b], 
			  Sb[a], Sb[b], ovlbbme[a+b*nb]);
	writeprojectedMBMEtoFile(mbbfile[a], mbbfile[b], &P, 
				 &OpOvlap, Sb[a], Sb[b], ovlbbme[a+b*nb]);
      }
    }

  // calculate overlaps between eigenstate and basis states

  complex double** ovlbme[n*nb];
  for (b=0; b<nb; b++)
    for (a=0; a<n; a++)
      ovlbme[a+b*n] = initprojectedMBME(&P, &OpOvlap);

  for (b=0; b<nb; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMBMEfromFile(mbfile[a], mbbfile[b], &P, 
				    &OpOvlap, S[a], Sb[b], ovlbme[a+b*n])) {
	calcprojectedMBME(&P, &OpOvlap, &Q[a], &Qb[b], 
			  S[a], Sb[b], ovlbme[a+b*n]);
	writeprojectedMBMEtoFile(mbfile[a], mbbfile[b], &P, 
				 &OpOvlap, S[a], Sb[b], ovlbme[a+b*n]);
      }
    }

  
  // evaluate matrix elements of one-operator
  // implent Symmetries !

  ManyBodyOperator OpOne = {
    name : "One",
    rank : 0,
    pi : 0,
    dim : 1,
    size : 1,
    par : NULL,
    me : NULL
  };

  // use full set of basis states (all J,M,K projections)
  // using K mixing eigenstates would require to solve the
  // eigenvalue problem for each basis state

  int jmax = P.jmax;
  
  complex double** oneme[n*n];
  for (b=0; b<n; b++)
    for (a=0; a<n; a++)
      oneme[a+b*n] = initprojectedMBME(&P, &OpOne);

  complex double *nbbwork = malloc(SQR((jmax+1)*nb)*sizeof(complex double));
  complex double *obbwork = malloc(SQR((jmax+1)*nb)*sizeof(complex double));
  complex double *obbcovlwork = malloc(SQR(jmax+1)*nb*n*sizeof(complex double));

  // loop over parity and angular momenta

  fprintf(stderr, "... preparing basis-one operator\n");
      
  int p,j,idx;
  for (p=0; p<=1; p++)
    for (j=P.odd; j<P.jmax; j=j+2) {
      idx = idxpij(P.jmax, p, j);

      // calculate overlap matrix

      fprintf(stderr, "\tj: %d, p: %d basis overlap matrix\n", j, p);

      complex double (*nbb)[(j+1)*nb] = nbbwork;
      int ab,bb,kab,kbb;

      for (bb=0; bb<nb; bb++)
	for (kbb=-j; kbb<=j; kbb=kbb+2)
	  for (ab=0; ab<nb; ab++)
	    for (kab=-j; kab<=j; kab=kab+2)
		 
	      nbb[idxnjm(ab,j,kab)][idxnjm(bb,j,kbb)] =
		ovlbbme[ab+nb*bb][idx][idxjmk(j,kab,kbb)];

      // invert overlap matrix
	  
      double thresh = 1e-5;
      complex double (*obb)[(j+1)*nb] = obbwork;
      pseudoinverse(nbb, obb, (j+1)*nb, thresh);
 
      // calculate one operator matrix element

      fprintf(stderr, "\tj: %d, p: %d basis one matrix\n", j, p);

      complex double (*obbcovl)[n*(j+1)] = obbcovlwork;

      int ka,kb;

      for (b=0; b<n; b++)
	for (ab=0; ab<nb; ab++)
	  for (kb=-j; kb<=j; kb=kb+2)
	    for (kab=-j; kab<=j; kab=kab+2) {
	      obbcovl[idxnjm(ab,j,kab)][idxnjm(b,j,kb)] = 0.0;

	      for (bb=0; bb<nb; bb++)
		for (kbb=-j; kbb<=j; kbb=kbb+2)
		  obbcovl[idxnjm(ab,j,kab)][idxnjm(b,j,kb)] +=
		    obb[idxnjm(ab,j,kab)][idxnjm(bb,j,kbb)]*
		    conj(ovlbme[b+n*bb][idx][idxjmk(j,kb,kbb)]);
	    }

      for (b=0; b<n; b++)
	for (a=0; a<n; a++)
	  for (kb=-j; kb<=j; kb=kb+2)
	    for (ka=-j; ka<=j; ka=ka+2) {
	      oneme[a+b*n][idx][idxjmk(j,ka,kb)] = 0.0;

	      for (ab=0; ab<nb; ab++)
		for (kab=-j; kab<=j; kab=kab+2)
		  oneme[a+b*n][idx][idxjmk(j,ka,kb)] +=
		    ovlbme[a+n*ab][idx][idxjmk(j,ka,kab)]*
		    obbcovl[idxnjm(ab,j,kab)][idxnjm(b,j,kb)];
	    }		  

    }

  // expectation values of one operator

  fprintf(stderr, "... calculate expectation values\n");
  
  complex double **oneexp = initprojectedVector(&P, &OpOne, n);
  calcexpectprojectedMBME(&P, &OpOne, oneme, S, &E, oneexp);

  // output

  char outfile[1024];
  FILE* outfp;

  snprintf(outfile, 1024, "%s.basisovlap-%s", 
	   stripstr(multiconfigfile, ".states"), basisfile);
  if (!(outfp = fopen(outfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", outfile);
    exit(-1);
  }

  fprintinfo(outfp);

  writeprojectedOvlaps(outfp, &P, oneexp, &E);

  fclose(outfp);

}

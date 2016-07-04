/**

  \file calcbasisovlapmulti.c

  calculate overlap of Many-Body eigenstates with 
  Basis spanned by a set of projected Many-Body states

  uses full set of basis states (all J,M,K projections)
  probably better to use only proper eigenstates


  (c) 2005-2009 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/MultiSlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Projection.h"

#include "misc/utils.h"
#include "misc/physics.h"

#include "numerics/cmat.h"

#define SQR(x) ((x)*(x))


#define MAXBASISSTATES 100


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
	fprintf(fp, "\n%s Specf = ", prefix); 
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
    fprintf(stderr, "\nusage: %s mcstate basistates\n"
            "-p ANGPARA         projection parameters for overlap matrix elements"
	    "-t THRESH		threshold for inversion\n", 
	    argv[0]);
    exit(-1);
  }

  char* angparovl = NULL;
  double thresh=1e-5;

  int c;
  while((c = getopt(argc, argv, "p:t:")) != -1)
    switch (c) {
    case 'p':
      angparovl = optarg;
      break;
    case 't':   
      thresh = atof(optarg);
      break;
    }

  char* multiconfigfile = argv[optind];
  char** mbfile;

  // open multiconfigfile
  Projection Pstate;
  Projection Povl;
  MultiSlaterDet* Q;
  Indices* In;
  Eigenstates E;
  int n;

  if (readMultiMulticonfigfile(multiconfigfile, &mbfile, &Pstate, &Q, &In, &n, &E))
    exit(-1);

  if (angparovl) {
    int odd = Q[0].A % 2;
    initProjection(&Povl, odd, angparovl);
  } else
    Povl = Pstate;

  char* basisfile = argv[optind+1];
  char* mbbfile[MAXBASISSTATES];
  int nb;

  if (readstringsfromfile(basisfile, &nb, mbbfile))
    exit(-1);

  MultiSlaterDet Qb[nb]; 
  Indices Inb[nb];

  int i;
  for (i=0; i<nb; i++) {
    extractIndicesfromString(&mbbfile[i], &Inb[i]);
    if (readMultiSlaterDetfromFile(&Qb[i], &Inb[i], mbbfile[i]))
      exit(-1);
  }


  // calculate overlaps for basis states
 
  complex double*** ovlbbme[nb*nb];

  int a,b;
  for (b=0; b<nb; b++)
    for (a=0; a<nb; a++)
      ovlbbme[a+b*nb] = initprojectedMultiMBME(&Povl, &OpOvlap, &Qb[a], &Qb[b]);

  for (b=0; b<nb; b++)
    for (a=0; a<nb; a++) {
      if (readprojectedMultiMBMEfromFile(mbbfile[a], mbbfile[b], 
					 &Qb[a], &Qb[b],
					 &Povl, &OpOvlap, ovlbbme[a+b*nb])) {
	calcprojectedMultiMBME(&Povl, &OpOvlap, &Qb[a], &Qb[b], ovlbbme[a+b*nb]);
	writeprojectedMultiMBMEtoFile(mbbfile[a], mbbfile[b], &Qb[a], &Qb[b],
				      &Povl, &OpOvlap, ovlbbme[a+b*nb]);
      }
    }

  // calculate overlaps between eigenstate and basis states

  complex double*** ovlbme[n*nb];
  for (b=0; b<nb; b++)
    for (a=0; a<n; a++)
      ovlbme[a+b*n] = initprojectedMultiMBME(&Povl, &OpOvlap, &Q[a], &Qb[b]);

  for (b=0; b<nb; b++)
    for (a=0; a<n; a++) {
      if (readprojectedMultiMBMEfromFile(mbfile[a], mbbfile[b], 
					 &Q[a], &Qb[b],
					 &Povl, &OpOvlap, ovlbme[a+b*n])) {
	calcprojectedMultiMBME(&Povl, &OpOvlap, &Q[a], &Qb[b], ovlbme[a+b*n]);
	writeprojectedMultiMBMEtoFile(mbfile[a], mbbfile[b], &Q[a], &Qb[b],
				      &Povl, &OpOvlap, ovlbme[a+b*n]);
      }
    }


  // evaluate matrix elements of one-operator
  // implement Symmetries !

  ManyBodyOperator OpOne = {
    name : "One",
    rank : 0,
    pi : 0,
    dim : 1,
    size: 1,
    par : NULL,
    me : NULL
  };

  // use full set of basis states (all J,M,K projections)
  // using K mixing eigenstates would require to solve the
  // eigenvalue problem for each basis state

  int jmax = Povl.jmax;

  complex double*** oneme[n*n];
  for (b=0; b<n; b++)
    for (a=0; a<n; a++)
      oneme[a+b*n] = initprojectedMultiMBME(&Povl, &OpOne, &Q[a], &Q[b]);

  int dim=0;
  for (i=0; i<n; i++)
    dim += In[i].N;

  // use only selected basis states
  int dimb=0;
  for (i=0; i<nb; i++)
    dimb += Inb[i].n;

  complex double *nbbwork = malloc(SQR((jmax+1)*dimb)*sizeof(complex double));
  complex double *obbwork = malloc(SQR((jmax+1)*dimb)*sizeof(complex double));
  complex double *obbcovlwork = malloc(SQR(jmax+1)*dimb*dim*sizeof(complex double));

  // loop over parity and angular momenta

  fprintf(stderr, "... preparing basis-one operator\n");

  int p,j,idx;
  for (p=0; p<=1; p++)
    for (j=Povl.odd; j<Povl.jmax; j=j+2) {
      idx = idxpij(Povl.jmax, p, j);

      // calculate overlap matrix

      fprintf(stderr, "\tj: %d, p: %d basis overlap matrix\n", j, p);

      complex double (*nbb)[(j+1)*dimb] = nbbwork;
      int ab,abp,idxab,bb,bbp,idxbb,kab,kbb;

      idxbb=-1;
      for (bb=0; bb<nb; bb++)
        for (bbp=0; bbp<Inb[bb].n; bbp++){
          idxbb++;
          for (kbb=-j; kbb<=j; kbb=kbb+2) {
            idxab=-1;
            for (ab=0; ab<nb; ab++)
              for (abp=0; abp<Inb[ab].n; abp++) {
                idxab++;
                for (kab=-j; kab<=j; kab=kab+2) {

                  nbb[idxjm(j,kab)+idxab*(j+1)][idxjm(j,kbb)+idxbb*(j+1)] = 
                    ovlbbme[ab+bb*nb][Inb[ab].idx[abp]+Inb[bb].idx[bbp]*Inb[ab].N][idx][idxjmk(j,kab,kbb)];

                }
              }
          }
        }
      
      // invert overlap matrix

      complex double (*obb)[(j+1)*dimb] = obbwork;
      pseudoinverse(nbb, obb, (j+1)*dimb, thresh);

      // calculate one operator matrix element

      fprintf(stderr, "\tj: %d, p: %d basis one matrix\n", j, p);

      complex double (*obbcovl)[dim*(j+1)] = obbcovlwork;

      int ap,bp,idxa,idxb,ka,kb;

      idxb=-1;
      for (b=0; b<n; b++)
        for (bp=0; bp<In[b].N; bp++) {
          idxb++;

          idxab=-1;
          for (ab=0; ab<nb; ab++)
            for (abp=0; abp<Inb[ab].n; abp++) {
              idxab++;
              
              for (kb=-j; kb<=j; kb=kb+2)
                for (kab=-j; kab<=j; kab=kab+2) {
                  obbcovl[idxjm(j,kab)+idxab*(j+1)][idxjm(j,kb)+idxb*(j+1)] = 0.0;

                  idxbb=-1;
                  for (bb=0; bb<nb; bb++)
                    for (bbp=0; bbp<Inb[bb].n; bbp++) {
                      idxbb++;
                      for (kbb=-j; kbb<=j; kbb=kbb+2)
                        obbcovl[idxjm(j,kab)+idxab*(j+1)][idxjm(j,kb)+idxb*(j+1)] +=
                          obb[idxjm(j,kab)+idxab*(j+1)][idxjm(j,kbb)+idxbb*(j+1)]*
                          conj(ovlbme[b+n*bb][bp+Inb[bb].idx[bbp]*In[b].N][idx][idxjmk(j,kb,kbb)]);
                    }
                }
            }
        }

      idxb=-1;
      for (b=0; b<n; b++)
        for (bp=0; bp<In[b].N; bp++) {
          idxb++;

          for (a=0; a<n; a++)
            for (ap=0; ap<In[a].N; ap++) {
              idxa++;

              for (kb=-j; kb<=j; kb=kb+2)
                for (ka=-j; ka<=j; ka=ka+2) {

                  oneme[a+b*n][ap+bp*In[a].N][idx][idxjmk(j,ka,kb)] = 0.0;

                  idxab=-1;
                  for (ab=0; ab<nb; ab++)
                    for (abp=0; abp<Inb[ab].n; abp++) {
                      idxab++;
                      for (kab=-j; kab<=j; kab=kab+2)

                        oneme[a+b*n][ap+bp*In[a].N][idx][idxjmk(j,ka,kb)] +=
                          ovlbme[a+n*ab][ap+Inb[ab].idx[abp]*In[a].N][idx][idxjmk(j,ka,kab)]*
                          obbcovl[idxjm(j,kab)+idxab*(j+1)][idxjm(j,kb)+idxb*(j+1)];
                    }
                }
            }
        }

    }


  // expectation values of one operator

  fprintf(stderr, "... calculate expectation values\n");

  complex double **oneexp = initprojectedVector(&Povl, &OpOne, dim);
  calcexpectprojectedMultiMBME(&Povl, &OpOne, oneme, Q, In, n, &E, oneexp);

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

  writeprojectedOvlaps(outfp, &Povl, oneexp, &E);

  fclose(outfp);

  return 0;
}

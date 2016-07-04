/**

  \file printovlmeproj-fmd-frozen.c

  multiconfiguration calculations with projected states


  (c) 2008 Thomas Neff

*/

/*
  07/14/09 do not multiply matrix elements by norm, Slater determinants are assumed to be
           properly normalized
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Observables.h"
#include "fmd/Symmetry.h"
#include "fmd/Projection.h"
#include "fmd/Ovlap.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"
#include "numerics/cmat.h"


#define MAXSTATES 250

char* angmomtostr(int j)
{
  char* str = malloc(4*sizeof(char));

  if (j%2) 
    sprintf(str, "%d/2", j);
  else
    sprintf(str, "%d", j/2);

  return str;
}


char* pitostr(int pi)
{
  char* str = malloc(2*sizeof(char));

  if (pi==0)
    sprintf(str, "+1");
  else
    sprintf(str, "-1");

  return str;
}


void extractmatrices(const char* nucsfile,
		     const Projection* P,
		     const Interaction* Int,
		     const complex double ***ovlme,
		     const Eigenstates* Ep,
		     int n,
		     int j, int p,
                     int diagonal)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int a, aa, b, bb;
  int ipj, i, m, k;
  int smalldim;
  int ai, bi, iai, ibi, idxa, idxb, idxai, idxbi;

  ipj=idxpij(jmax,p,j);

  // how many good basis states do we really have ?
  smalldim=0;
  idxa=-1;
  for (a=0; a<n; a++) {
    idxa++;
    smalldim += Ep[idxa].ngood[ipj];
  }

  // do we have at least a one-dimensional space
  if (smalldim == 0) {
    fprintf(stderr, "no basis states for j: %d, pi: %d\n", j, p);
  } else {
    fprintf(stderr, "%d basis states for j: %d, pi: %d\n", smalldim, j, p);

    complex double N[smalldim*smalldim];

    idxb=-1; idxbi=-1;
    for (b=0; b<n; b++) {
      idxb++;
      for (bi=0; bi<Ep[idxb].ngood[ipj]; bi++) {
        idxbi++;
        ibi=Ep[idxb].index[ipj][bi];
        idxa=-1; idxai=-1;
        for (a=0; a<n; a++) {

          if (diagonal && a != b)
            continue;

          idxa++;
          for (ai=0; ai<Ep[idxa].ngood[ipj]; ai++) {
            idxai++;
            iai=Ep[idxa].index[ipj][ai];
										
            N[idxai+idxbi*smalldim] = 0.0;

            for (k=-j; k<=j; k=k+2)
              for (m=-j; m<=j; m=m+2) {
                      
                N[idxai+idxbi*smalldim] +=
                  conj(Ep[idxa].V[ipj][idxjm(j,m)+iai*(j+1)])*
                  ovlme[a+b*n][ipj][idxjmk(j,m,k)]*
                  Ep[idxb].V[ipj][idxjm(j,k)+ibi*(j+1)];
			 
              }
          }
        }
      }	
    }

    char outfilename[255];
    FILE* outfp;

    snprintf(outfilename, 255, "%s.N%s-%d%c.bin", 
	     filepart(nucsfile), diagonal ? "matrixd" : "matrix", j, p ? '-' : '+');
    outfp = fopen(outfilename, "w");
    fwritecmatbin(outfp, smalldim, N);
    fclose(outfp);
  }		      
}



int main(int argc, char* argv[])
{
  createinfo(argc, argv);

  /* enough arguments ? */

  if (argc < 6) {
    fprintf(stderr, "\nusage: %s PROJPAR INTERACTION FMDNUCSFILE CLUSTERNUCSFILE OUTFILENAME\n"
            "   -d              diagonal matrix elements only\n"
	    "   -n MINNORM\n"
	    "   -t THRESHOLD\n",
	    argv[0]);
    return -1;
  }

  int odd;
  int diagonal=0;
  double threshkmix=0.1;
  double minnormp=0.001;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "dn:t:")) != -1)
    switch (c) {
    case 'd':
      diagonal = 1;
      break;
    case 'n':
      minnormp = atof(optarg);
      break;
    case 't':
      threshkmix = atof(optarg);
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* fmdnucsfile = argv[optind+2];
  char* clusternucsfile = argv[optind+3];
  char* outfilename = argv[optind+4];

  char* mbfile[MAXSTATES];
  int nfmd, nfrozen;
  int n;

  if (readstringsfromfile(fmdnucsfile, &nfmd, &mbfile[0]))
    return -1;

  if (readstringsfromfile(clusternucsfile, &nfrozen, &mbfile[nfmd]))
    return -1;

  n = nfmd+nfrozen;

  if (n>MAXSTATES) {
    fprintf(stderr, "abort: too many states, increase MAXSTATES\n");
    exit(-1);
  }

  SlaterDet Q[n]; 
  Symmetry S[n];
  
  int i;
  for (i=0; i<n; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      return -1;
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    return -1;
  Int.cm = 1;

  // odd number of nucleons ?
  odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);
    
  initOpObservables(&Int);

  int a,b; 


  // initialize space for matrix elements
  complex double** ovlme[n*n];
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++)
      ovlme[a+b*n] = initprojectedMBME(&P, &OpOvlap);

  // read or calculate matrix elements
  for (b=0; b<n; b++)
    for (a=0; a<n; a++) {

      if (diagonal && a != b)
        continue;

      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], 
				    &P, &OpOvlap, 
				    S[a], S[b], ovlme[a+b*n])) {
	calcprojectedMBME(&P, &OpOvlap, 
			  &Q[a], &Q[b], S[a], S[b], ovlme[a+b*n]);

	writeprojectedMBMEtoFile(mbfile[a], mbfile[b],
				 &P, &OpOvlap, 
				 S[a], S[b], ovlme[a+b*n]);
      }
    }

  // observables are need to determine eigenstates
  Observablesod** obsme[n];
  for (a=0; a<n; a++)
    obsme[a] = initprojectedMBME(&P, &OpObservables);

  // read or calculate matrix elements
  for (a=0; a<n; a++) {

    if (readprojectedMBMEfromFile(mbfile[a], mbfile[a], 
                                  &P, &OpObservables, 
                                  S[a], S[a], obsme[a])) {
      calcprojectedMBME(&P, &OpObservables, 
                        &Q[a], &Q[a], S[a], S[a], obsme[a]);

      writeprojectedMBMEtoFile(mbfile[a], mbfile[a],
                               &P, &OpObservables, 
                               S[a], S[a], obsme[a]);
    }
  }

  // scale matrix elements
  if (Int.mescaling) {
    for (a=0; a<n; a++)
      scaleprojectedObservablesMBME(&P, &Int, obsme[a]);
  }

  // dimension of many-body basis
  // beware: matrix elements are calculated in larger basis

  int dimfmd=0, dimfrozen=0, dim=0;
  dimfmd = nfmd;
  dimfrozen = nfrozen;
  dim = dimfmd+dimfrozen;


  // read or calculate the Eigenstates
      
  int allp=0;
  Eigenstates Ep[dim];
  Observablesod** obsp = initprojectedVector(&P, &OpObservables, 1);
  
  for (i=0; i<n; i++) {
    if (readEigenstatesfromFile(mbfile[i], &P, &Ep[i], 1)) {
      fprintf(stderr, "... calculating Eigenstates for %s\n", mbfile[i]);
      initEigenstates(&P, &Ep[i], 1);
      calcEigenstates(&P, &Int, &obsme[i], &Ep[i], threshkmix);
      calcexpectprojectedMBME(&P, &OpObservables, &obsme[i], &S[i], &Ep[i], obsp);
      sortEigenstates(&P, &Int, obsp, &Ep[i], minnormp, allp);
    }
  }


  // Mathematica input
  char mmafilename[255];
  snprintf(mmafilename, 255, "%s.mma", outfilename); 
  
  FILE* mmafp;
  mmafp = fopen(mmafilename, "w");

  fprintf(mmafp, "(* parameters written by printobsmeproj-fmd-frozen *)\n\n");

  fprintf(mmafp, "nfmd=%d;", nfmd); 
  fprintf(mmafp, "nfrozen=%d;\n", nfrozen);

  fprintf(mmafp, "\n");

  // filenames and md5s and symmetries
  
  for (i=0; i<n; i++)
    fprintf(mmafp, "mbfile[%d]=\"%s\"; md5sum[%d] = \"%s\"; sym[%d] = %d;\n", 
	    i+1, mbfile[i], i+1, md5hash(mbfile[i]), i+1, S[i]);

  fprintf(mmafp, "\n");

  // dimensions and K-mixing coefficients

  int pi, j, ipj, k;
  int adim, totaldim;
  int aa, idxa, ai, iai, idxai;

  fprintf(mmafp, "dfmd=%d;", dimfmd);
  fprintf(mmafp, "dfrozen=%d;\n\n", dimfrozen);

  for (pi=0; pi<=1; pi++)
    for (j=P.odd; j<P.jmax; j=j+2) {
      ipj = idxpij(P.jmax, pi, j);

      totaldim=0;
      idxa=-1;
      for (a=0; a<nfmd; a++) {

	  idxa++;
	  adim = Ep[idxa].ngood[ipj];
	  fprintf(mmafp, "dim[%d][%s][%s]=%d;\n",
		  idxa+1, angmomtostr(j), pitostr(pi), adim); 

	totaldim += adim;
      }
      fprintf(mmafp, "dimfmd[%s][%s]=%d;\n", angmomtostr(j), pitostr(pi), totaldim); 

      totaldim=0;
      for (a=0; a<nfrozen; a++) {

	  idxa++;
	  adim = Ep[idxa].ngood[ipj];
	  if (adim != 1)
	    fprintf(stderr, "warning: %d eigenstates for frozen config %d\n",
		    adim, a+1);

	totaldim += adim;
      }
      fprintf(mmafp, "dimfrozen[%s][%s]=%d;\n\n", angmomtostr(j), pitostr(pi), totaldim); 

    }  

  fprintf(mmafp, "\n");

  // print K-mixing coefficients for Mathematica 

  for (pi=0; pi<=1; pi++)
    for (j=P.odd; j<P.jmax; j=j+2) {
      ipj = idxpij(P.jmax, pi, j);

      idxa=-1;
      for (a=0; a<nfmd; a++) {

	  idxa++;
	  for (ai=0; ai<Ep[idxa].ngood[ipj]; ai++) {
	    iai=Ep[idxa].index[ipj][ai];	  

	    fprintf(mmafp, "Kmixfmd[%d][%d][%s][%s]={", 
		    idxa+1, ai+1, angmomtostr(j), pitostr(pi));	

	    for (k=-j; k<=j; k=k+2)
	      fprintf(mmafp, " %+.8f %+.8f I%s",
		      creal(Ep[idxa].V[ipj][idxjm(j,k)+iai*(j+1)]),
		      cimag(Ep[idxa].V[ipj][idxjm(j,k)+iai*(j+1)]),
		      k==j ? "};\n" : ",");		    
	  }

      }

      // frozen configs K-mixing, assume to be identical for all frozen states

      a = nfmd;

	idxa++;
	for (ai=0; ai<Ep[idxa].ngood[ipj]; ai++) {
	  idxai++;
	  iai=Ep[idxa].index[ipj][ai];	  

	  fprintf(mmafp, "Kmixfrozen[%s][%s]={", angmomtostr(j), pitostr(pi));
	  for (k=-j; k<=j; k=k+2)
	    fprintf(mmafp, " %+.8f %+.8f I%s",
		    creal(Ep[idxa].V[ipj][idxjm(j,k)+iai*(j+1)]),
		    cimag(Ep[idxa].V[ipj][idxjm(j,k)+iai*(j+1)]),
		    k==j ? "};\n\n" : ",");		    
	}

    }
  
  fclose(mmafp);


  // extract the matrix elements and write to file

  for (pi=0; pi<=1; pi++)
    for (j=odd; j<P.jmax; j=j+2)
      extractmatrices(outfilename, &P, &Int, ovlme, Ep, n, j, pi, diagonal);

  return 0;
}

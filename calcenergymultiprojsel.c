/**

  \file calcenergymultiprojsel.c

  iteratively select subset of many-body states
  that give the lowest energy

  we might work in a large set of states selecting only a small number of them
  memory consumption becomes a problem, 

  (c) 2006 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Ovlap.h"
#include "fmd/Observables.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"
#include "fmd/ProjectedObservables.h"

#include "misc/utils.h"
#include "misc/physics.h"

#include "numerics/cmat.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/Projectionmpi.h"
#include "fmdmpi/ProjectionSlave.h"
#endif

#define MAXSTATES 500

#define EUNDEFINED 1000.0

inline int SQR(int i) { return i*i; }


// extract energy of eigenstate
double energy(const Projection* P, const Eigenstates* E,
	      int j, int p, int i)
{
  int ipj = idxpij(P->jmax,p,j);
  int ngood = E->ngood[ipj];

  fprintf(stderr, "\n");

  if (!ngood)
    return EUNDEFINED;

  if (i>=ngood) 
    i=ngood-1;

  return E->v[ipj][E->index[ipj][i]];
}


// calculate the ovlap of the n-th config with the n-1 configs for j,p
double calcbasisovlap(const Projection* P, const Eigenstates *E, 
		      const complex double ***ovlmes, 
		      int n, int j, int p)
{
  int ipj = idxpij(P->jmax,p,j);
 
  complex double* nbb = malloc(SQR((n-1)*(j+1))*sizeof(complex double)); 
  complex double* obb = malloc(SQR((n-1)*(j+1))*sizeof(complex double));

  int a,b,ia,idxa,iai,ib,idxb,ibi,ka,kb,dim;

  dim=0;
  for (a=0; a<n-1; a++)
    dim += E[a].ngood[ipj];

  for (idxb=0; idxb<dim; idxb++)
    for (idxa=0; idxa<dim; idxa++)
      nbb[idxa+idxb*dim] = 0.0;
      
  // n-1 overlap matrix

  idxb=-1;
  for (b=0; b<n-1; b++)
    for (ib=0; ib<E[b].ngood[ipj]; ib++) {
      ibi = E[b].index[ipj][ib];
      idxb++;
      idxa=-1;
      for (a=0; a<n-1; a++)
	for (ia=0; ia<E[a].ngood[ipj]; ia++) {
	  iai = E[a].index[ipj][ia];
	  idxa++;
	  for (kb=0; kb<j+1; kb++)
	    for (ka=0; ka<j+1; ka++)
	      nbb[idxa+idxb*dim] += 
		conj(E[a].V[ipj][ka+iai*(j+1)])*
		ovlmes[a+b*n][ipj][ka+kb*(j+1)]*
		E[b].V[ipj][kb+ibi*(j+1)];
	}
    }

  // one operator

  double thresh=1e-9;
  pseudoinverse(nbb, obb, dim, thresh);

  // calculate one operator matrix element
	  
  int i,ii,k,kp;
  double ovlmax = 0.0, ovl;

  // return biggest overlap
  for (i=0; i<E[n-1].ngood[ipj]; i++) {

    ovl = 0.0;

    ii = E[n-1].index[ipj][i];
    idxb=-1;
    for (b=0; b<n-1; b++)
      for (ib=0; ib<E[b].ngood[ipj]; ib++) {
	ibi = E[b].index[ipj][ib];
	idxb++;
	idxa=-1;
	for (a=0; a<n-1; a++)
	  for (ia=0; ia<E[a].ngood[ipj]; ia++) {
	    iai = E[a].index[ipj][ia];
	    idxa++;

	    for (kb=0; kb<j+1; kb++)
	      for (ka=0; ka<j+1; ka++)
		for (kp=0; kp<j+1; kp++)
		  for (k=0; k<j+1; k++)
		    
		    ovl += conj(E[n-1].V[ipj][k+ii*(j+1)])*
		      ovlmes[(n-1)+a*n][ipj][k+ka*(j+1)]*
		      E[a].V[ipj][ka+iai*(j+1)]*
		      obb[idxa+idxb*dim]*
		      conj(E[b].V[ipj][kb+ibi*(j+1)])*
		      ovlmes[b+(n-1)*n][ipj][kb+kp*(j+1)]*
		      E[n-1].V[ipj][kp+ii*(j+1)]/
		      E[n-1].norm[ipj][ii];

	  }
      }

    if (ovl > ovlmax)
      ovlmax = ovl;
  }
			 
  free(nbb);
  free(obb);

  return ovlmax;
}


static int cmpmerit(void* ap, void* bp) 
{
  double a=*(double*) ap;
  double b=*(double*) bp;
  
  return (a<=b ? (a<b ? 1 : 0) : -1);
}

double calcenergyeigenvalue(const Projection* P, const Eigenstates *E, 
			    const Observablesod*** obsme, 
			    int n, int j, int p, int ei,
                            double thresh, double minnorm)
{
  int ipj=idxpij(P->jmax,p,j);
  int dim, i;
  double norma2, normb2;

  dim=0;
  for (i=0; i<n; i++)
    dim += E[i].ngood[ipj];

  // do we have at least a one-dimensional space
  if (dim == 0) {
    return EUNDEFINED;
  } else {

    complex double *H = malloc(SQR(n*(j+1))*sizeof(complex double));
    complex double *N = malloc(SQR(n*(j+1))*sizeof(complex double));
    complex double *v = malloc(n*(j+1)*sizeof(complex double));
    complex double *V = malloc(SQR(n*(j+1))*sizeof(complex double));
  
    int a, b;
    int m, k, l;
    int d;
    int ai, bi, iai, ibi, idxa, idxb;
 
    idxb=-1;
    for (b=0; b<n; b++)
      for (bi=0; bi<E[b].ngood[ipj]; bi++) {
	ibi=E[b].index[ipj][bi];
        normb2=E[b].norm[ipj][ibi];
	idxb++;
	idxa=-1;
	for (a=0; a<n; a++)	
	  for (ai=0; ai<E[a].ngood[ipj]; ai++) {
	    iai=E[a].index[ipj][ai];
            norma2=E[a].norm[ipj][iai];
	    idxa++;
	    N[idxa+idxb*dim] = 0.0;
	    H[idxa+idxb*dim] = 0.0;

	    for (k=-j; k<=j; k=k+2)
	      for (m=-j; m<=j; m=m+2) {
		N[idxa+idxb*dim] +=
		  conj(E[a].V[ipj][idxjm(j,m)+iai*(j+1)])*
		  obsme[a+b*n][ipj][idxjmk(j,m,k)].n*
		  E[b].V[ipj][idxjm(j,k)+ibi*(j+1)]/
                  sqrt(norma2*normb2);
			 
		H[idxa+idxb*dim] +=
		  conj(E[a].V[ipj][idxjm(j,m)+iai*(j+1)])*
		  obsme[a+b*n][ipj][idxjmk(j,m,k)].h*
		  E[b].V[ipj][idxjm(j,k)+ibi*(j+1)]/
                  sqrt(norma2*normb2);

	      }	
	  }	
      }      

    generalizedeigensystem(H, N, dim, thresh, v, V, &d);

    // we got d eigenstates, filter out the noise
    // see sortEigenstates

    complex double norm[d];
    for (i=0; i<d; i++) {
      norm[i] = 0.0;
      for (l=0; l<dim; l++)
        for (k=0; k<dim; k++)
          norm[i] += conj(V[k+i*dim])*N[k+l*dim]*V[l+i*dim];
    }

    double maxnorm = 0.0;
    for (i=0; i<d; i++)
      maxnorm = fmax(maxnorm, creal(norm[i]));

    // for debugging
    /*
    for (i=0; i<d; i++)
      fprintf(stderr, "i: %d, e[i] = %8.3f MeV, n[i] = %7.5f\n", i, hbc*creal(v[i]), norm[i]);
    */

    // merits

    struct merit {
      double val;
      int idx;
    } merits[d];

    for (i=0; i<d; i++) {
      merits[i].idx = i;
      if (creal(norm[i]) < minnorm*maxnorm)
        merits[i].val = -100000.0;
      else
        merits[i].val = -creal(v[i]);
    }

    // sort according to merits value
    qsort(merits, d, sizeof(struct merit), cmpmerit);

    i=0; while (i<d && merits[i].val > -100000.0) i++;

    if (i < ei)
      return EUNDEFINED;

    double e = -merits[ei].val;

    free(H);
    free(N);
    free(v);
    free(V);

    return e;
  }
}


void readorcalcandwriteprojectedMBMEfromtoFile(const char* mbfile, const char* mbfilep,
					       const Projection* P, 
					       const ManyBodyOperator* Op,
					       const SlaterDet* Q, const SlaterDet* Qp, 
					       Symmetry S, Symmetry Sp,
					       const void** mbme)
{
  // allocate space
  *mbme = initprojectedMBME(P, Op);

  if (readprojectedMBMEfromFile(mbfile, mbfilep, P, Op, S, Sp, *mbme)) {
#ifdef MPI
    calcprojectedMBMEmpi(P, Op, Q, Qp, S, Sp, *mbme); 
#else
    calcprojectedMBME(P, Op, Q, Qp, S, Sp, *mbme); 
#endif
    writeprojectedMBMEtoFile(mbfile, mbfilep, P, Op, S, Sp, *mbme); 
  }

  // scale interaction matrix elements

  if (!strncmp(Op->name, "Observables", 11)) {
    Interaction* Int = Op->par;

    if (Int->mescaling) {
      scaleprojectedObservablesMBME(P, Int, *mbme);
    }
  }
        
}


void cleanup(int ret)
{
#ifdef MPI
  int task=TASKFIN;
  BroadcastTask(&task);

  MPI_Finalize();
#endif

  exit(ret);
}


int main(int argc, char* argv[])
{
  createinfo(argc, argv);

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  // fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    ProjectionSlave();

    MPI_Finalize();
  } else {
#endif


  /* enough arguments ? */

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPAR INTERACTION NUCSFILE"
	    "\n   -f FIXNUCSFILE"
	    "\n   -n N           search only for up to N states"
	    "\n   -j JSEL        select for angular momentum JSEL"
	    "\n   -p PSEL        select for parity PSEL [0|1]"
	    "\n   -i ISEL        select for #ISEL"
	    "\n   -o OVLTHRESH   overlap threshold"
	    "\n   -e ENTHRESH    energy threshold [MeV]"
	    "\n   -t THRESH      threshold for SVD\n", argv[0]);
    exit(-1);
  }

  int odd;
  double threshkmix=0.01;
  double minnormkmix=0.001;
  int all=0;
  double threshmulti=0.0000001;
  double minnormmulti=0.001;
  double ovlthresh=0.95;
  double enthresh=0.050/hbc;
  int nmax=0;
  // optimize for which quantum numbers
  int jsel=0, psel=0, isel=0;

  char* fixnucsfile = NULL;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "f:t:n:j:p:i:o:e:")) != -1)
    switch (c) {
    case 'f':
      fixnucsfile = optarg;
      break;
    case 't':
      threshmulti = atof(optarg);
      break;
    case 'n':
      nmax = atoi(optarg);
      break;
    case 'j':
      jsel = atoi(optarg);
      break;
    case 'p':
      psel = atoi(optarg);
      break;
    case 'i':
      isel = atoi(optarg);
      break;
    case 'o':
      ovlthresh = atof(optarg);
      break;
    case 'e':
      enthresh = atof(optarg)/hbc;
      break;
    }

  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];
  char* nucsfile = argv[optind+2];

  char* mbfile[MAXSTATES];
  int fixn=0, searchn, n, nsel;
  int i;

  // are there fixed configs ?
  if (fixnucsfile) {
    if (readstringsfromfile(fixnucsfile, &fixn, mbfile)) {
      fprintf(stderr, "couldn't open %s\n", fixnucsfile);
      cleanup(-1);
    }
  }
  nsel = fixn;

  // possibly there are identical states in nucsfile and fixnucsfile
  // they should be eliminated by the overlap threshold
  if (readstringsfromfile(nucsfile, &searchn, &mbfile[fixn])) {
    fprintf(stderr, "couldn't open %s\n", nucsfile);
    cleanup(-1);
  }
 
  // total number of states
  n = fixn + searchn;

  // which many-body states are selected
  int idx[n];
  
  // selected (best) energies for set of i states
  double esel[n];

  // maximum number of states to select
  if (nmax==0) 
    nmax=n;
  else 
    nmax=fixn+nmax;


  SlaterDet Q[n]; 
  Symmetry S[n];

  for (i=0; i<n; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i])) {
      fprintf(stderr, "couldn't read from %s\n", mbfile[i]);
      cleanup(-1);
    }
  }

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile)) {
    fprintf(stderr, "couldn't read Interaction from %s\n", interactionfile);
    cleanup(-1);
  }
  Int.cm = 1; 

  // odd number of nucleons ?
  odd = Q[0].A % 2;

  // integer or half-integer jsel - compatible with odd ?
  if (odd != jsel %2) {
    fprintf(stderr, "not possible to project on j=%d\n", jsel);
    cleanup(-1);
  }

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastInteraction(&Int);
  BroadcastA(&Q[0].A);
#endif

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);

  initOpObservables(&Int);

  int a,b; 

  // matrix elements
  complex double** ovlme[n*n];
  Observablesod** obsme[n*n];

  // keep track of already calculated matrix elements
  int ovlmedone[n*n];
  int obsmedone[n*n];

  // initialize space for matrix elements only on demand
  for (b=0; b<n; b++)	
    for (a=0; a<n; a++) {
      ovlmedone[a+b*n] = 0; obsmedone[a+b*n] = 0;
      ovlme[a+b*n] = NULL; obsme[a+b*n] = NULL;	
      // ovlme[a+b*n] = initprojectedMBME(&P, &OpOvlap);
      // obsme[a+b*n] = initprojectedMBME(&P, &OpObservables);
    }

  // we read/calculate the diagonal matrix elements right at the beginning

  for (a=0; a<n; a++) {

    readorcalcandwriteprojectedMBMEfromtoFile(mbfile[a], mbfile[a], &P, &OpOvlap, 
					      &Q[a], &Q[a], S[a], S[a], &ovlme[a+a*n]);
    ovlmedone[a+a*n] = 1;
    
    readorcalcandwriteprojectedMBMEfromtoFile(mbfile[a], mbfile[a], &P, &OpObservables, 
					      &Q[a], &Q[a], S[a], S[a], &obsme[a+a*n]);
    obsmedone[a+a*n] = 1;
  }

  // read or calculate the Eigenstates
      
  int allp=0;
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


  // common root for all output files
  char outfile[1024];
  if (fixn > 0)
    sprintf(outfile, "%s--%s--sel-%d-%d-%d", fixnucsfile, nucsfile, jsel, psel, isel);
  else
    sprintf(outfile, "%s--sel-%d-%d-%d", nucsfile, jsel, psel, isel);

  // the log file
  char nucssellogfile[1024];
  sprintf(nucssellogfile, "%s.log", outfile);

  FILE* logfp;
  if (!(logfp = fopen(nucssellogfile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", nucssellogfile);
    cleanup(-1);
  }  
  setlinebuf(logfp);
  fprintinfo(logfp);

  fprintf(logfp, "\n");
  fprintf(logfp, "ovlap threshold :  %5.3f\n", ovlthresh);
  fprintf(logfp, "energy threshold:  %5.3f MeV\n", hbc*enthresh);


  // Eselp and obsmesel are pointers to the selected
  // Eigenstates and matrixelements
  Symmetry Ssel[n];
  Eigenstates Eselp[n] ;
  complex double** ovlmesel[n*n];
  Observablesod** obsmesel[n*n];
  Observablesod** obs;
  
  Eigenstates multiE;
  Amplitudes multiA;

  initEigenstates(&P, &multiE, n);
  initAmplitudes(&P, &multiA, n);
  obs = initprojectedVector(&P, &OpObservables, n);

  int imin;
  double e, ovl, emin = EUNDEFINED;

  if (fixn) {

    // are there fixed states ?

    for (i=0; i<fixn; i++)
      idx[i] = i;

    multiE.n = fixn;
    multiA.n = fixn;

    // read/calc overlap matrix elements 
    for (b=0; b<fixn; b++)
      for (a=0; a<fixn; a++) {
	if (!ovlmedone[idx[a]+idx[b]*n]) {
	  readorcalcandwriteprojectedMBMEfromtoFile(mbfile[idx[a]], mbfile[idx[b]], 
						    &P, &OpOvlap, 
						    &Q[idx[a]], &Q[idx[b]], 
						    S[idx[a]], S[idx[b]], 
						    &ovlme[idx[a]+idx[b]*n]);
	  ovlmedone[idx[a]+idx[b]*n] = 1;	  
	}
	if (!obsmedone[idx[a]+idx[b]*n]) {
	  readorcalcandwriteprojectedMBMEfromtoFile(mbfile[idx[a]], mbfile[idx[b]], 
						    &P, &OpObservables, 
						    &Q[idx[a]], &Q[idx[b]], 
						    S[idx[a]], S[idx[b]], 
						    &obsme[idx[a]+idx[b]*n]);
	  obsmedone[idx[a]+idx[b]*n] = 1;	  
	}
      }    

    for (a=0; a<fixn; a++) {
      Ssel[a] = S[idx[a]];
      Eselp[a] = Ep[idx[a]];
    }
		      
    for (b=0; b<fixn; b++)
      for (a=0; a<fixn; a++) {
	ovlmesel[a+b*fixn] = ovlme[idx[a]+idx[b]*n];
	obsmesel[a+b*fixn] = obsme[idx[a]+idx[b]*n];
      }

    e = calcenergyeigenvalue(&P, Eselp, obsmesel, fixn, jsel, psel, isel, threshmulti, minnormmulti);

    esel[fixn-1] = e;

    fprintf(logfp, "\n ... taking %2d fix configs\n\n", fixn);
    fprintf(logfp, "\n                      energy: %8.3f MeV\n", hbc*e);

  } else {

    // find state with lowest energy

    imin=-1;
    for (a=0; a<n; a++) {
      e=energy(&P, &Ep[a], jsel, psel, isel);

      fprintf(stderr, "[%2d] - %s - energy:  %8.3f MeV\n",
	      a, mbfile[a], hbc*e);

      if (e < emin) {
	emin = e;
	imin = a;
      }
    }
    
    idx[0] = imin;
    esel[0] = emin;
    nsel = 1;

    fprintf(logfp, "\n\n ... selecting  1 config\n\n");
    fprintf(logfp, "selected [%2d] - %s - energy:  %8.3f MeV\n", 
	    imin, mbfile[imin], hbc*emin);
  }
   
  // now select additional configs

  int done=0;
  while (!done && nsel < nmax) {

    fprintf(logfp, "\n\n ... selecting %2d configs\n\n", nsel+1);

    // set correct dimension for eigenstates and amplitudes
    multiE.n = nsel+1;
    multiA.n = nsel+1;

    for (a=0; a<nsel; a++) {
      Ssel[a] = S[idx[a]];
      Eselp[a] = Ep[idx[a]];
    }
		      
    for (b=0; b<nsel; b++)
      for (a=0; a<nsel; a++) {
	ovlmesel[a+b*(nsel+1)] = ovlme[idx[a]+idx[b]*n];
	obsmesel[a+b*(nsel+1)] = obsme[idx[a]+idx[b]*n];
      }

    // find the config which lowers the energy the most and has not too much overlap
    // with already selected configurations

    imin = -1;
    emin = EUNDEFINED;
    for (i=fixn; i<n; i++) {
      
      // is i already selected ?
      for (int j=0; j<nsel; j++)
	if (idx[j] == i)
	  goto next;

      // read/calc overlap matrix elements 
      for (a=0; a<nsel; a++) {
	if (!ovlmedone[idx[a]+i*n]) {
	  readorcalcandwriteprojectedMBMEfromtoFile(mbfile[idx[a]], mbfile[i], 
						    &P, &OpOvlap, 
						    &Q[idx[a]], &Q[i], S[idx[a]], S[i], 
						    &ovlme[idx[a]+i*n]);
	  ovlmedone[idx[a]+i*n] = 1;	  
	}
	if (!ovlmedone[i+idx[a]*n]) {
	  readorcalcandwriteprojectedMBMEfromtoFile(mbfile[i], mbfile[idx[a]], 
						    &P, &OpOvlap, 
						    &Q[i], &Q[idx[a]], S[i], S[idx[a]], 
						    &ovlme[i+idx[a]*n]);
	  ovlmedone[i+idx[a]*n] = 1;	  
	}
      }

      Ssel[nsel] = S[i];
      Eselp[nsel] = Ep[i];
      for (a=0; a<nsel; a++) {
	ovlmesel[a+nsel*(nsel+1)] = ovlme[idx[a]+i*n];
	ovlmesel[nsel+a*(nsel+1)] = ovlme[i+idx[a]*n];
      }
      ovlmesel[nsel+nsel*(nsel+1)] = ovlme[i+i*n];

      // if overlap with already selected basis states to big skip
      ovl = calcbasisovlap(&P, Eselp, ovlmesel, nsel+1, jsel, psel);
      fprintf(logfp, "[%2d]     ovlap:    %6.3f\t", i, ovl);
      if (ovl > ovlthresh) {
	fprintf(logfp, "\n");
	goto next;
      }

      // read/calc observables matrix elements 
      for (a=0; a<nsel; a++) {
	if (!obsmedone[idx[a]+i*n]) {
	  readorcalcandwriteprojectedMBMEfromtoFile(mbfile[idx[a]], mbfile[i], 
						    &P, &OpObservables, 
						    &Q[idx[a]], &Q[i], S[idx[a]], S[i], 
						    &obsme[idx[a]+i*n]);
	  obsmedone[idx[a]+i*n] = 1;	  
	}
	if (!obsmedone[i+idx[a]*n]) {
	  readorcalcandwriteprojectedMBMEfromtoFile(mbfile[i], mbfile[idx[a]], 
						    &P, &OpObservables, 
						    &Q[i], &Q[idx[a]], S[i], S[idx[a]], 
						    &obsme[i+idx[a]*n]);
	  obsmedone[i+idx[a]*n] = 1;	  
	}
      }
      
      for (a=0; a<nsel; a++) {
	obsmesel[a+nsel*(nsel+1)] = obsme[idx[a]+i*n];
	obsmesel[nsel+a*(nsel+1)] = obsme[i+idx[a]*n];
      }
      obsmesel[nsel+nsel*(nsel+1)] = obsme[i+i*n];

      e = calcenergyeigenvalue(&P, Eselp, obsmesel, nsel+1, jsel, psel, isel, threshmulti, minnormmulti);

      fprintf(logfp, "        energy:  %8.3f MeV\n", hbc*e);

      if (e < emin) {
	imin = i;
	emin = e;
      }
      
    next:
      ;
    }
    
    if (emin < esel[nsel-1]-enthresh) {
      idx[nsel] = imin;
      esel[nsel] = emin;

      fprintf(logfp, "\nselected [%2d] - %s - energy: %8.3f MeV\n",
	      imin, mbfile[imin], hbc*emin);

      Ssel[nsel+1] = S[imin];
      Eselp[nsel+1] = Ep[imin];
      for (a=0; a<nsel; a++) {
	ovlmesel[a+nsel*(nsel+1)] = ovlme[idx[a]+imin*n];
	ovlmesel[nsel+a*(nsel+1)] = ovlme[imin+idx[a]*n];
      }
      ovlmesel[nsel+nsel*(nsel+1)] = ovlme[imin+imin*n];
      for (a=0; a<nsel; a++) {
	obsmesel[a+nsel*(nsel+1)] = obsme[idx[a]+imin*n];
	obsmesel[nsel+a*(nsel+1)] = obsme[imin+idx[a]*n];
      }
      obsmesel[nsel+nsel*(nsel+1)] = obsme[imin+imin*n];

      nsel++;
      
      // solve the eigenvalue problem
      calcMultiEigenstates(&P, &Int, obsmesel, &Eselp, &multiE, &multiA, threshmulti);

      // calculate expectation values
      calcexpectprojectedMBME(&P, &OpObservables, obsmesel, Ssel, &multiE, obs);

      // sort Eigenstates
      sortEigenstates(&P, &Int, obs, &multiE, minnormmulti, all);


      // create nucsfile with selected configs
      
      char nucsselfile[1024];
      sprintf(nucsselfile, "%s--%d", outfile, nsel);
      FILE *outfp;
      if (!(outfp = fopen(nucsselfile, "w"))) {
	fprintf(stderr, "couldn't open %s for writing\n", nucsselfile);
	cleanup(-1);
      }
      for (a=0; a<nsel; a++)
	if (S[idx[a]])
	  fprintf(outfp, "%s:%s\n", SymmetrytoStr(S[idx[a]]), mbfile[idx[a]]);
	else
	  fprintf(outfp, "%s\n", mbfile[idx[a]]);
      fclose(outfp);

      char nucsselprojfile[1024];
      sprintf(nucsselprojfile, "%s.multi-%s", nucsselfile, ProjectiontoStr(&P));
      if (!(outfp = fopen(nucsselprojfile, "w"))) {
	fprintf(stderr, "couldn't open %s for writing\n", nucsselprojfile);
	cleanup(-1);
      }
      fprintinfo(outfp);
      fprintProjectinfo(outfp, &P);
      showprojectedObservables(outfp, &P, &Int, &Q[0], obs, &multiE, &multiA, "");
      fclose(outfp);
      
    } else {

      fprintf(logfp, "\n\n... no many-body state lowers energy enough, stopping here\n");
      // indicate we have finished
      char nucsselfinfile[1024];
      sprintf(nucsselfinfile, "%s.finished", outfile);

      FILE* finfp;
      if (!(finfp = fopen(nucsselfinfile, "w"))) {
	break;
      }
      fclose(finfp);

      done = 1;
    }

  }

  fclose(logfp);

  cleanup(0);

#ifdef MPI
  }
#endif

}



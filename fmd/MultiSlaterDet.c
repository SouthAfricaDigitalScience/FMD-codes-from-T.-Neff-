/**

   \file MultiSlaterDet.c

   set of Many-body states described by multiple Slaterdets

   
   (c) 2005

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#include "SlaterDet.h"
#include "Symmetry.h"
#include "Observables.h"

#include "MultiSlaterDet.h"
#include "SimpleSlaterDet.h"
#include "SymmetricSlaterDet.h"
#include "SymmetricMultiSlaterDet.h"
#include "DiClusterProjSlaterDet.h"
#include "DiClusterMultiProjSlaterDet.h"
#include "DiClusterMulticonfig.h"

#include "misc/utils.h"
#include "numerics/wignerd.h"
#include "numerics/clebsch.h"
#include "numerics/zcw.h"
#include "numerics/gaussquad.h"
#include "numerics/cmat.h"


#define SQR(x) ((x)*(x))

void extractIndicesfromString(char** str, Indices* Ind)
{
  int i;

  Ind->N = 0;
  Ind->n = 0;
  for (i=0; i<NMAX; i++)
    Ind->idx[i] = 0;

  char *c, *s = *str;

  i = 0;
  while(1) {
    c = strtok(s, ": >");

    if (c == NULL) {
      *str = NULL;
      return;
    }

    if (isdigit(*c)) 
      Ind->idx[i++] = atoi(c);
    else {
      Ind->n = i;
      *str = c;
      return;
    }
    s = NULL;
  }
}
 

char* IndicestoStr(const Indices* Ind)
{
  char* str;

  str = malloc(20*sizeof(char));
  str[0] = '\0';

  // all indices means empty string
  if (Ind->N == Ind->n)
    return str;
    
  int i;
  for (i=0; i<Ind->n; i++)
    sprintf(str, "%s%d:", str, i);
  
  str[strlen(str)-1] = '\0';
  return str;
}



#define BUFSIZE 4096

int readMultiSlaterDetfromFile(MultiSlaterDet* MB, Indices* In, 
			       const char* fname)
{
  FILE* fp;

  if (!(fp = fopen(fname, "r"))) {  
    fprintf(stderr, "couldn't open %s for reading\n", fname);
    return -1;
  }

  fprintf(stderr, "... reading MultiSlaterDet from file %s\n", fname);

  char buf[BUFSIZE];

  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<SlaterDet", 10) && 
	 strncmp(buf, "<MultiSlaterDet", 14) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <MultiSlaterDet ...>\n");
    return -1;
  }

  // simple SlaterDet
  if (!strncmp(buf, "<SlaterDet>", 11)) {
    rewind(fp);
    if (SimpleSlaterDetRead(fp, MB))
      return -1;
  }

  else {
    char Mname[80];
    sscanf(buf, "<MultiSlaterDet %s>", Mname);
    stripstr(Mname, ">");

    if (!strcmp(Mname, "SymmetricSlaterDet")) {
      if (SymmetricSlaterDetRead(fp, MB))
	return -1;
    }    
    else if (!strcmp(Mname, "SymmetricMultiSlaterDet")) {
      if (SymmetricMultiSlaterDetRead(fp, MB))
	return -1;
    }    
    else if (!strcmp(Mname, "DiClusterProjSlaterDet")) {
      if (DiClusterProjSlaterDetRead(fp, MB))
	return -1;
    }
    else if (!strcmp(Mname, "DiClusterMultiProjSlaterDet")) {
      if (DiClusterMultiProjSlaterDetRead(fp, MB))
	return -1;
    }
    else if (!strcmp(Mname, "DiClusterMulticonfig")) {
      if (DiClusterMulticonfigRead(fp, MB))
	return -1;
    }
    else {
      fprintf(stderr, "MultiSlaterDet %s not known !\n", Mname);
      return -1;
    }
  }

  In->N = MB->N;

  // no indices selected, then take all possible ones
  if (In->n == 0) {
    In->n = In->N;
    int i;
    for (i=0; i<In->N; i++)
      In->idx[i] = i;
  }
  
  if (In->n > In->N) {
    fprintf(stderr, "too many Indices\n");
    return -1;
  }

  return 0;
}
    

void* initprojectedMultiMBME(const Projection* P, const ManyBodyOperator* Op,
			     const MultiSlaterDet* MBA, const MultiSlaterDet* MBB)
{
  int NA = MBA->N;
  int NB = MBB->N;

  void** me = malloc(NA*NB*sizeof(void*));
  
  int iA, iB;
  for (iB=0; iB<NB; iB++)
    for (iA=0; iA<NA; iA++)
      me[iA+iB*NA] = initprojectedMBME(P, Op);

  return me;
}


int readprojectedMultiMBMEfromFile(const char* mbfilea, const char* mbfileb,
				   const MultiSlaterDet* MBA,
				   const MultiSlaterDet* MBB,
				   const Projection* P,
				   const ManyBodyOperator* Op,
				   void** mbme)
{
  gzFile mefp;
  char mefilename[255];
  char buf[BUFSIZE];
  char fnam[255], md5fnam[33];
  char fileS[80];

  snprintf(mefilename, 255, "ME/%s--%s--%s--%s.gz", 
	   Op->name, 
	   filepart(mbfilea), 
	   filepart(mbfileb), 
	   ProjectiontoStr(P));

  if (fileexists(mefilename)) {
      fprintf(stderr, "... %s does not exist\n", mefilename);
      return -1;
    }

  if (!(mefp = gzopen(mefilename, "r"))) {
    fprintf(stderr, "... couldn't open %s for reading\n", mefilename);
    return -1;
  }
  fprintf(stderr, "... reading matrix elements from file %s\n", mefilename);

  // do many-body files match ?
  // agreement of Symmetries not checked
  do
    gzgets(mefp, buf, BUFSIZE);
  while (strncmp(buf, "<MBFile ", 7) && !gzeof(mefp));
  if (gzeof(mefp)) {
    fprintf(stderr, "...   did't find <MBFile ...>\n");
    gzclose(mefp);
    return -2;
  }
  sscanf(buf, "<MBFile %s %s %s>", fnam, md5fnam, fileS);
  if (strncmp(md5fnam, md5hash(mbfilea), 32)) {
    fprintf(stderr, "...    MBFile does not match with existing ME\n");
    return -2;
  }

  do
    gzgets(mefp, buf, BUFSIZE);
  while (strncmp(buf, "<MBFile ", 7) && !gzeof(mefp));
  if (gzeof(mefp)) {
    fprintf(stderr, "...    did't find <MBFile ...>\n");
    gzclose(mefp);
    return -2;
  }
  sscanf(buf, "<MBFile %s %s %s>", fnam, md5fnam, fileS);
  if (strncmp(md5fnam, md5hash(mbfileb), 32)) {
    fprintf(stderr, "...    MBFile does not match with existing ME\n");
    gzclose(mefp);
    return -2;
  }

  // now read the matrix elements
  int NA = MBA->N; int NB = MBB->N;
  int iA, iB;
  Symmetry SA, SB;
  int res = 0;

  for (iB=0; iB<NB; iB++) {
    SB = MBB->symmetry(MBB, iB);
    for (iA=0; iA<NA; iA++) {
      SA = MBA->symmetry(MBA, iA); 
      res &= readprojectedMBME(mefp, P, Op, SA, SB, mbme[iA+iB*NA]); 
    }
  }

  gzclose(mefp);
  return res;
}

static inline double dmin(double a, double b)
{
  return (a < b ? a : b);
}

static inline double dmax(double a, double b)
{
  return (a > b ? a : b);
}

static double dminarray(double a[], int n)
{
  double m = a[0];
  int i;
  for (i=1; i<n; i++)
    m = dmin(m, a[i]);

  return m;
}

static double dmaxarray(double a[], int n)
{
  double m = a[0];
  int i;
  for (i=1; i<n; i++)
    m = dmax(m, a[i]);

  return m;
}


#define MAX(a,b) ((a)>(b) ? (a) : (b))

// 07/13/09 imortant change: do not normalize individual SlaterDets anymore

void calcprojectedMultiMBME(const Projection* P, const ManyBodyOperator* Op,
			    const MultiSlaterDet* MBA, 
			    const MultiSlaterDet* MBB,
			    void** mbme)
{
  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  complex double (***val)[(rank+1)*size] = mbme;

  int jmax = P->jmax;
  int odd = P->odd;

  // loop over the individual SlaterDets in the MultiSlaterDets
  int nA = MBA->n; int nB = MBB->n;
  int iA, iB;

  int NA = MBA->N; int NB = MBB->N;
  int IA, IB;

  SlaterDet Q, Qp, Qpp;
  allocateSlaterDet(&Q, MBA->A);
  allocateSlaterDet(&Qp, MBB->A);
  allocateSlaterDet(&Qpp, MBB->A);

  SlaterDetAux X;
  allocateSlaterDetAux(&X, MAX(MBA->A, MBB->A));

  Symmetry S, Sp;
  int axialsym;

  int l, r;
  int p, j, m, k;


  // set matrix elements to zero
  for (IB=0; IB<NB; IB++)
    for (IA=0; IA<NA; IA++)
      for (p=0; p<=1; p++)
	for (j=odd; j<jmax; j=j+2) {
	  for (k=-j; k<=j; k=k+2)
	    for (m=-j; m<=j; m=m+2)
	      for (l=0; l<dim; l++)
		for (r=0; r<=rank; r++)
		  val[IA+IB*NA][idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = 0.0;
	}			


  // we need a common denominaotor for the NA*NB matrix elements
  // all states have axial symmetry ? use axial as indicator
  // also spherical will be tretated like axial

  S=0;
  axialsym=1;
  for (IA=0; IA<NA; IA++)
    axialsym &= hasAxialSymmetry(MBA->symmetry(MBA, IA));
  if (axialsym)
    setSymmetry(&S, axial);

  Sp=0;
  axialsym=1;
  for (IB=0; IB<NB; IB++)
    axialsym &= hasAxialSymmetry(MBB->symmetry(MBB, IB));
  if (axialsym)
    setSymmetry(&Sp, axial);

  fprintf(stderr, "Symmetry - MBA: %s, MBB: %s\n", SymmetrytoStr(S), SymmetrytoStr(Sp)); 

  double kappaA[nA], kappaB[nB];
  double acmA[nA], acmB[nB];
  
  for (iA=0; iA<nA; iA++) {
    MBA->get(MBA, iA, &Q);
    kappaA[iA] = _estimateangkappa(&Q);
    acmA[iA] = _estimateacm(&Q);
  }

  for (iB=0; iB<nB; iB++) {
    MBB->get(MBB, iB, &Qp);
    kappaB[iB] = _estimateangkappa(&Qp);
    acmB[iB] = _estimateacm(&Qp);
  }

  fprintf(stderr, "kappacrit: %6.2f\n", _getangkappacrit());
  fprintf(stderr, "kappa - MBA: [%6.2f - %6.2f], MBB: [%6.2f - %6.2f]\n",
          dminarray(kappaA, nA), dmaxarray(kappaA, nA),
          dminarray(kappaB, nB), dmaxarray(kappaB, nB));


  int c=0; 
  int cmax=nA*nB; 

  for (iB=0; iB<nB; iB++) {
    MBB->get(MBB, iB, &Qp);

    for (iA=0; iA<nA; iA++) {
      MBA->get(MBA, iA, &Q);

      // progress indicator
      c++;
      if (c%100==0) fprintf(stderr, "%d%%", (100*c)/cmax);
      if (c%10==0) fprintf(stderr, ".");

      // norms of SlaterDets
      // calcSlaterDetAuxod(&Q, &Q, &X);
      // double norm = sqrt(creal(X.ovlap));

      // calcSlaterDetAuxod(&Qp, &Qp, &X);
      // double normp = sqrt(creal(X.ovlap));  

      // set up cm integration
      cmintegrationpara cmpara;
      double cmalpha = 0.5/(acmA[iA]+acmB[iB]);
      _initcmintegration(P, cmalpha, &cmpara);
      int ncm = cmpara.n;

      // set up ang integration
      angintegrationpara angpara;
      double angkappa = dmin(kappaA[iA], kappaB[iB]);
      _initangintegration(P, angkappa, S, Sp, &angpara);
      int nang = angpara.n;

      int icm; 
      double xcm[3]; double weightcm;

      int iang;
      double alpha, beta, gamma; double weightang;

      complex double sval[(rank+1)*size];
      double weight;
      int ip;

      for (icm=0; icm<ncm; icm++) {
	getcmintegrationpoint(icm, &cmpara, xcm, &weightcm);

	for (iang=0; iang<nang; iang++) {
	  getangintegrationpoint(iang, &angpara, &alpha, &beta, &gamma, &weightang);
	  // weight = 1.0/(2*norm*normp)*weightcm*weightang;
          weight = 0.5*weightcm*weightang;
      
	  copySlaterDet(&Qp, &Qpp);
	  moveSlaterDet(&Qpp, xcm);
	  rotateSlaterDet(&Qpp, alpha, beta, gamma);

	  for (ip=0; ip<=1; ip++) {
	    if (ip) invertSlaterDet(&Qpp);

	    // can only calculate Auxiliaries if Sldets are compatible
	    if (Q.A == Qp.A && Q.Z == Qp.Z && Q.N == Qp.N)
	      calcSlaterDetAuxod(&Q, &Qpp, &X);
	    Op->me(Op->par, &Q, &Qpp, &X, sval);

	    complex double w, wA, wB;
	    Symmetry SA, SB;
	    for (IB=0; IB<NB; IB++) {
	      SB = MBB->symmetry(MBB, IB);
	      wB = MBB->weight(MBB, IB, iB);
	      for (IA=0; IA<NA; IA++) {
		SA = MBA->symmetry(MBA, IA);
		wA = MBA->weight(MBA, IA, iA);
		for (p=0; p<=1; p++)
		  for (j=odd; j<jmax; j=j+2)
		    for (k=-j; k<=j; k=k+2)
		      for (m=-j; m<=j; m=m+2) {
			if ((Op->rank != 0 || SymmetryAllowed(SA, p, j, m)) &&
			    SymmetryAllowed(SB, p, j, k)) {
			  w = conj(wA)*wB*	 
			    weight * (p && ip%2 ? -1 : 1)*
			    (j+1)/(8*SQR(M_PI))*Djmkstar(j,m,k,alpha,beta,gamma);
			  for (l=0; l<dim; l++)
			    for (r=0; r<=rank; r++)
			      val[IA+IB*NA][idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] += w*sval[r+l*(rank+1)];
		    }	
		  }
	      }
	    }	
	  }

	}

      }	

      freeAngintegration(&angpara);
      freecmintegration(&cmpara);
      
    }
  }
}


int writeprojectedMultiMBMEtoFile(const char* mbfilea, const char* mbfileb,
				  const MultiSlaterDet* MBA,
				  const MultiSlaterDet* MBB,
				  const Projection* P,
				  const ManyBodyOperator* Op,
				  void** mbme)
{
  gzFile mefp;
  char mefilename[255];

  ensuredir("ME");

  snprintf(mefilename, 255, "ME/%s--%s--%s--%s.gz", 
	   Op->name, 
	   filepart(mbfilea), 
	   filepart(mbfileb), 
	   ProjectiontoStr(P));

  if (!(mefp = gzopen(mefilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", mefilename);
    return -1;
  }
    
  gzprintinfo(mefp);

  int NA = MBA->N; int NB = MBB->N;
  int IA, IB;
  Symmetry SA, SB;
  char syma[80] = "", symb[80] = "";
  
  for (IA=0; IA<NA; IA++) {
    SA = MBA->symmetry(MBA, IA);
    sprintf(syma, "%s%d:", syma, SA);
  }
  syma[strlen(syma)-1] = '\0';

  for (IB=0; IB<NB; IB++) {
    SB = MBB->symmetry(MBB, IB);
    sprintf(symb, "%s%d:", symb, SB);
  }
  symb[strlen(symb)-1] = '\0';

  gzprintf(mefp, "<MBFile %s %s %s>\n", 
	   mbfilea, md5hash(mbfilea), syma);
  gzprintf(mefp, "<MBFile %s %s %s>\n", 
	   mbfileb, md5hash(mbfileb), symb);

  // write matrix elements

  int res = 0;

  for (IB=0; IB<NB; IB++) {
    SB = MBB->symmetry(MBB, IB);
    for (IA=0; IA<NA; IA++) {
      SA = MBA->symmetry(MBA, IA); 

      gzprintf(mefp, "\n# Matrixelement :%d  :%d\n", IA, IB);
      res &= writeprojectedMBME(mefp, P, Op, SA, SB, mbme[IA+IB*NA]); 
    }
  }

  gzclose(mefp);
  return res;
} 


void calcMultiEigenstatesMulti(const Projection* P,
			       const Interaction* Int,
			       const Observablesod ****obsme,
			       const Eigenstates* Ep,
			       const Indices* In,
			       int n,
			       Eigenstates* multiE, Amplitudes* multiA, 
			       double thresh)
{
  int dim=multiE->n;
  int odd=P->odd;
  int jmax=P->jmax;

  complex double* H = malloc(SQR(dim*(jmax+1))*sizeof(complex double));
  complex double* N = malloc(SQR(dim*(jmax+1))*sizeof(complex double));
  complex double* v = malloc(dim*(jmax+1)*sizeof(complex double));
  complex double* V = malloc(SQR(dim*(jmax+1))*sizeof(complex double));
  
  int a, aa, b, bb;
  int p, j, ipj, i, m, k;
  int smalldim, fulldim, d;
  int ai, bi, iai, ibi, idxa, idxb, idxai, idxbi;
  double norma2, normb2, normi2;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      fulldim = dim*(j+1);

      ipj=idxpij(jmax,p,j);

      // how many good basis states do we really have ?
      smalldim=0;
      idxa=-1;
      for (a=0; a<n; a++)
	for (aa=0; aa<In[a].n; aa++) {
	  idxa++;
	  smalldim += Ep[idxa].ngood[ipj];
	}

      // do we have at least a one-dimensional space
      if (smalldim == 0) {
	multiE->dim[ipj] = 0;
      } else {

	idxb=-1; idxbi=-1;
	for (b=0; b<n; b++)
	  for (bb=0; bb<In[b].n; bb++) {
	    idxb++;
	    for (bi=0; bi<Ep[idxb].ngood[ipj]; bi++) {
	      ibi=Ep[idxb].index[ipj][bi];
              normb2=Ep[idxb].norm[ipj][ibi];
	      idxbi++;
	      idxa=-1; idxai=-1;
	      for (a=0; a<n; a++)
		for (aa=0; aa<In[a].n; aa++) {
		  idxa++;
		  for (ai=0; ai<Ep[idxa].ngood[ipj]; ai++) {
		    iai=Ep[idxa].index[ipj][ai];
                    norma2=Ep[idxa].norm[ipj][iai];
		    idxai++;
		    N[idxai+idxbi*smalldim] = 0.0;
		    H[idxai+idxbi*smalldim] = 0.0;

		    for (k=-j; k<=j; k=k+2)
		      for (m=-j; m<=j; m=m+2) {
			N[idxai+idxbi*smalldim] +=
			  conj(Ep[idxa].V[ipj][idxjm(j,m)+iai*(j+1)])*
			  obsme[a+b*n][In[a].idx[aa]+In[b].idx[bb]*In[a].N][ipj][idxjmk(j,m,k)].n*
			  Ep[idxb].V[ipj][idxjm(j,k)+ibi*(j+1)]/
                          sqrt(norma2*normb2);
			 
			H[idxai+idxbi*smalldim] +=
			  conj(Ep[idxa].V[ipj][idxjm(j,m)+iai*(j+1)])*
			  obsme[a+b*n][In[a].idx[aa]+In[b].idx[bb]*In[a].N][ipj][idxjmk(j,m,k)].h*
			  Ep[idxb].V[ipj][idxjm(j,k)+ibi*(j+1)]/
                          sqrt(norma2*normb2);	

		      }
		  }
		}
	    }	
	  }		      

	generalizedeigensystem(H, N, smalldim, thresh, 
			       v, V, 
			       &d);

	// embed solution into full space

	// real dimension of eigenspace, may be smaller than smalldim
	multiE->dim[ipj] = d;

	for (i=0; i<d; i++) {
	  multiE->v[ipj][i] = v[i];

	  idxa=-1;
	  for (a=0; a<n; a++)
	    for (aa=0; aa<In[a].n; aa++) {
	      idxa++;
	      for (k=-j; k<=j; k=k+2)
		multiE->V[ipj][idxjm(j,k)+idxa*(j+1)+i*fulldim] = 0.0;
	    }

	  idxa=-1; idxai=-1;
	  for (a=0; a<n; a++)
	    for (aa=0; aa<In[a].n; aa++) {
	      idxa++;
	      for (ai=0; ai<Ep[idxa].ngood[ipj]; ai++) {
		idxai++;
		iai = Ep[idxa].index[ipj][ai];
                norma2 = Ep[idxa].norm[ipj][iai];

		for (k=-j; k<=j; k=k+2) {

		  multiE->V[ipj][idxjm(j,k)+idxa*(j+1)+i*fulldim] += 
                    Ep[idxa].V[ipj][idxjm(j,k)+iai*(j+1)]/sqrt(norma2)*
                    V[idxai+i*smalldim];
		}
	      }				      
	    }
	}

	// calculate Amplitudes

        for (i=0; i<d; i++) {

	  normi2 = 0.0;
	  for (idxbi=0; idxbi<smalldim; idxbi++)
	    for (idxai=0; idxai<smalldim; idxai++)
	      normi2 += conj(V[idxai+i*smalldim])*
		N[idxai+idxbi*smalldim]*V[idxbi+i*smalldim];

	  multiE->norm[ipj][i] = normi2;

	  idxa=-1; idxai=-1;
	  for (a=0; a<n; a++)
	    for (aa=0; aa<In[a].n; aa++) {
	      idxa++;

	      multiA->ngood[ipj][idxa] = Ep[idxa].ngood[ipj];
	  
	      for (ai=0; ai<Ep[idxa].ngood[ipj]; ai++) {
		idxai++;
		norma2 = N[idxai+idxai*smalldim];
		multiA->amp[ipj][ai+idxa*(j+1)+i*fulldim] = 0.0;
		for (idxbi=0; idxbi<smalldim; idxbi++)
		  multiA->amp[ipj][ai+idxa*(j+1)+i*fulldim] +=
		    N[idxai+idxbi*smalldim]*V[idxbi+i*smalldim]/
		    sqrt(norma2*normi2);
	      }
	    }
	}
        
      }

   }

  free(H); free(N);
  free(v); free(V);
}


// calculates reduced matrix element divided by sqrt(2j+1)

void calcexpectprojectedMultiMBME(const Projection* P,
				  const ManyBodyOperator* Op,
				  const void* mbme,
				  const MultiSlaterDet* Q,
				  const Indices* In,
				  int n,
				  const Eigenstates* E,
				  void* expectmbme)
{
  int p,j,i;
  int k, k1, k2;
  int l, r, nu;
  int idx;
  int rank=Op->rank;
  int dim=Op->dim;
  int size=Op->size;
  int odd=P->odd;
  int jmax=P->jmax;
  complex double (****me)[(rank+1)*size] = mbme;
  complex double (**expme)[size] = expectmbme;
 
  int dimj;
  int a, b, aa, bb;
  int idxa, idxb;
  int Idxaa, Idxbb;
  Symmetry Sa, Sb;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {

      idx = idxpij(jmax,p,j);
      dimj=E->n*(j+1);
      for (i=0; i<E->dim[idx]; i++) {

        for (l=0; l<dim; l++)
          expme[idx][i][l] = 0.0;
        
	idxb=-1;
	for (b=0; b<n; b++)
	  for (bb=0; bb<In[b].n; bb++) {
	    idxb++;
	    Idxbb = In[b].idx[bb];
	    Sb = Q[b].symmetry(&Q[b], Idxbb);
	    for (k2=-j; k2<=j; k2=k2+2) {
	      idxa=-1;
	      for (a=0; a<n; a++)
		for (aa=0; aa<In[a].n; aa++) {
		  idxa++;
		  Idxaa = In[a].idx[aa];
		  Sa = Q[a].symmetry(&Q[a], Idxaa);
		  for (k1=-j; k1<=j; k1=k1+2)	    
		    for (k=max(-j,k1-rank); k<=min(j,k1+rank); k=k+2) {
		      nu=k1-k;
		      if (SymmetryAllowed(Sa, p, j, k1) &&
			  SymmetryAllowed(Sb, p, j, k2)) {
			for (l=0; l<dim; l++)
			  expme[idx][i][l] +=
			    clebsch(j, rank, j, k, nu, k1)*
			    me[a+b*n][Idxaa+Idxbb*In[a].N][idx][idxjmk(j,k,k2)][(nu+rank)/2+l*(rank+1)]*
			    conj(E->V[idx][idxjm(j,k1) + idxa*(j+1) + i*dimj])*
			    E->V[idx][idxjm(j,k2)+ idxb*(j+1) + i*dimj];
		      }
		    }
		}		
	    }
	  }
      }
    }			
}


// calculate transition matrix elements for Eigenstates Efin <- Eini

// calculates reduced matrix element divided by sqrt(2jfin+1)

// slow but transparent
/*
void calctransitionprojectedMultiMBME(const Projection* P,
				      const ManyBodyOperator* Op,
				      const void* mbme,
				      const MultiSlaterDet* Qfin,
				      const MultiSlaterDet* Qini,
				      const Indices* Infin,
				      const Indices* Inini,
				      int nfin, int nini,
				      const Eigenstates* Efin,
				      const Eigenstates* Eini,
				      void* transmbme)
{
  int pini,jini,iini, pfin,jfin,ifin;
  int k, kfin, kini;
  int l, nu;
  int p=Op->pi;
  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  int oddini=P->odd;
  int jmax=P->jmax;
  complex double (****me)[(rank+1)*size] = mbme;
  complex double (****transme)[size] = transmbme;

  int afin, aini;
  int aafin, aaini;
  int dimini, dimfin;
  int idxini, idxfin;
  int idxaini, idxafin;
  int Idxaaini, Idxaafin;
  Symmetry Sini, Sfin;

  for (pini=0; pini<=1; pini++)
    for (jini=oddini; jini<jmax; jini=jini+2) {
      
      idxini = idxpij(jmax,pini,jini);
      dimini=Eini->n*(jini+1);

      pfin=(pini+p)%2;
      for (jfin=abs(jini-rank); jfin<=min(jmax-1,jini+rank); jfin=jfin+2) {

	idxfin = idxpij(jmax,pfin,jfin);
	dimfin=Efin->n*(jfin+1);
	for (iini=0; iini<Eini->dim[idxini]; iini++)
	  for (ifin=0; ifin<Efin->dim[idxfin]; ifin++) {

	    for (l=0; l<dim; l++)
	      transme[idxfin][idxini][ifin][iini][l] = 0.0;

	    idxaini=-1;
	    for (aini=0; aini<nini; aini++)
	      for (aaini=0; aaini<Inini[aini].n; aaini++) {
		idxaini++;
		Idxaaini = Inini[aini].idx[aaini];
		Sini = Qini[aini].symmetry(&Qini[aini], Idxaaini);

		for (kini=-jini; kini<=jini; kini=kini+2) {
		  
		  idxafin=-1;
		  for (afin=0; afin<nfin; afin++)
		    for (aafin=0; aafin<Infin[afin].n; aafin++) {
		      idxafin++;
		      Idxaafin = Infin[afin].idx[aafin];
		      Sfin = Qfin[afin].symmetry(&Qfin[afin], Idxaafin);
		      
		      for (kfin=-jfin; kfin<=jfin; kfin=kfin+2)
			for (k=max(-jini,kfin-rank); k<=min(jini,kfin+rank); k=k+2) {
			  nu = kfin-k;
			  if (SymmetryAllowed(Sfin, pfin, jfin, kfin) &&
                              SymmetryAllowed(Sini, pini, jini, kini)) {
			    for (l=0; l<dim; l++)
			      transme[idxfin][idxini][ifin][iini][l] +=
				clebsch(jini, rank, jfin, k, nu, kfin)*
				me[afin+aini*nfin][Idxaafin+Idxaaini*Infin[afin].N][idxini][idxjmk(jini,k,kini)][(nu+rank)/2+l*(rank+1)]*
				conj(Efin->V[idxfin][idxjm(jfin,kfin)+idxafin*(jfin+1)+ifin*dimfin])*
				Eini->V[idxini][idxjm(jini,kini)+idxaini*(jini+1)+iini*dimini];
			  }
			}
		    }
		}			
	      }
	  }
      }
    }
}
*/

// much faster, calculate first product of matrix elements with initial vectors
// then product with final vectors in a second step

void calctransitionprojectedMultiMBME(const Projection* P,
				      const ManyBodyOperator* Op,
				      const void* mbme,
				      const MultiSlaterDet* Qfin,
				      const MultiSlaterDet* Qini,
				      const Indices* Infin,
				      const Indices* Inini,
				      int nfin, int nini,
				      const Eigenstates* Efin,
				      const Eigenstates* Eini,
				      void* transmbme)
{
  int pini,jini,iini, pfin,jfin,ifin;
  int k, kfin, kini;
  int l, nu;
  int p=Op->pi;
  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  int oddini=P->odd;
  int jmax=P->jmax;
  complex double (****me)[(rank+1)*size] = mbme;
  complex double (****transme)[size] = transmbme;

  int afin, aini;
  int aafin, aaini;
  int dimini, dimfin;
  int idxini, idxfin;
  int idxaini, idxafin;
  int Idxaaini, Idxaafin;
  Symmetry Sini, Sfin;

  complex double* meVwork = malloc(Efin->n*Eini->n*(jmax+1)*(jmax+1)*size*(rank+1)*sizeof(complex double));

  for (pini=0; pini<=1; pini++)
    for (jini=oddini; jini<jmax; jini=jini+2) {
      
      idxini = idxpij(jmax,pini,jini);
      dimini=Eini->n*(jini+1);

      pfin=(pini+p)%2;
      for (jfin=abs(jini-rank); jfin<=min(jmax-1,jini+rank); jfin=jfin+2) {

	idxfin = idxpij(jmax,pfin,jfin);
	dimfin=Efin->n*(jfin+1);

        complex double (*meV)[dimini][(jini+1)][size*(rank+1)] = meVwork;

	for (iini=0; iini<Eini->dim[idxini]; iini++) {
          idxafin=-1;
	  for (afin=0; afin<nfin; afin++)
            for (aafin=0; aafin<Infin[afin].n; aafin++) {
              idxafin++;
              for (k=-jini; k<=jini; k=k+2)
                for (l=0; l<dim; l++)
                  for (nu=-rank; nu<=rank; nu=nu+2)
                    meV[idxafin][iini][idxjm(jini,k)][(nu+rank)/2+l*(rank+1)] = 0.0;
            }
        }

	for (iini=0; iini<Eini->dim[idxini]; iini++) {
          idxafin=-1;
	  for (afin=0; afin<nfin; afin++)
            for (aafin=0; aafin<Infin[afin].n; aafin++) {
              idxafin++;
              Idxaafin = Infin[afin].idx[aafin];
              Sfin = Qfin[afin].symmetry(&Qfin[afin], Idxaafin);

              for (k=-jini; k<=jini; k=k+2) {
                idxaini=-1;
                for (aini=0; aini<nini; aini++)
                  for (aaini=0; aaini<Inini[aini].n; aaini++) {
                    idxaini++;
                    Idxaaini = Inini[aini].idx[aaini];
                    Sini = Qini[aini].symmetry(&Qini[aini], Idxaaini);

                    for (kini=-jini; kini<=jini; kini=kini+2)
                      if (SymmetryAllowed(Sini, pini, jini, kini)) {
	      
                        for (l=0; l<dim; l++)
                          for (nu=-rank; nu<=rank; nu=nu+2)
                            meV[idxafin][iini][idxjm(jini,k)][(nu+rank)/2+l*(rank+1)] +=
                              me[afin+aini*nfin][Idxaafin+Idxaaini*Infin[afin].N][idxini][idxjmk(jini,k,kini)][(nu+rank)/2+l*(rank+1)]*
                              Eini->V[idxini][idxjm(jini,kini)+idxaini*(jini+1)+iini*dimini];
                      } 
		  }     
              }
            }
        }

	for (iini=0; iini<Eini->dim[idxini]; iini++)
	  for (ifin=0; ifin<Efin->dim[idxfin]; ifin++) {

	    for (l=0; l<dim; l++)
	      transme[idxfin][idxini][ifin][iini][l] = 0.0;

            idxafin=-1;
            for (afin=0; afin<nfin; afin++)
              for (aafin=0; aafin<Infin[afin].n; aafin++) {
                idxafin++;
                Idxaafin = Infin[afin].idx[aafin];
                Sfin = Qfin[afin].symmetry(&Qfin[afin], Idxaafin);
		      
                for (kfin=-jfin; kfin<=jfin; kfin=kfin+2)
                  for (k=max(-jini,kfin-rank); k<=min(jini,kfin+rank); k=k+2) {
                    nu = kfin-k;
                    if (SymmetryAllowed(Sfin, pfin, jfin, kfin)) {
                      for (l=0; l<dim; l++)
                        transme[idxfin][idxini][ifin][iini][l] +=
                          clebsch(jini, rank, jfin, k, nu, kfin)*
                          conj(Efin->V[idxfin][idxjm(jfin,kfin)+idxafin*(jfin+1)+ifin*dimfin])*
                          meV[idxafin][iini][idxjm(jini,k)][(nu+rank)/2+l*(rank+1)];
                    }     
                  }
	      }  
	  }     
      }
    }

  free(meVwork);
}


int writeMultiEigenstatestoFile(const char* fname,
				const Projection* P,
				const char* mbfile,
				const Indices* In,
				const Eigenstates* E)
{
  FILE* fp;
  char filename[255];

  snprintf(filename, 255, "%s.estates", fname);

  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", filename);
    return -1;
  }	

  fprintinfo(fp);
  fprintProjectinfo(fp, P);
  fprintf(fp, "\n");
  fprintf(fp, "<Projected %s>\n", ProjectiontoStr(P));

  fprintf(fp, "<MultiSlaterDet %s %s %s>\n",
	  mbfile, md5hash(mbfile), IndicestoStr(In));
  fprintf(fp, "\n");
    
  writeEigenstates(fp, P, E);

  fprintf(fp, "</Projected>\n");
  fclose(fp);

  return 0;
}


int writeMultiMulticonfigfile(const char* fname,
			      const Projection* P,
			      const char** mbfile,
			      const Indices* In,
			      int n,
			      const Eigenstates* E)
{
  FILE* fp;
  char filename[255];

  snprintf(filename, 255, "%s.states", fname);

  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", filename);
    return -1;
  }	

  fprintinfo(fp);
  fprintProjectinfo(fp, P);
  fprintf(fp, "\n");
  fprintf(fp, "<Multiconfprojected %d %d %s>\n", n, E->n, ProjectiontoStr(P));

  int i;
  for (i=0; i<n; i++)
    fprintf(fp, "<MBFile %s %s %s>\n",
	    mbfile[i], md5hash(mbfile[i]), IndicestoStr(&In[i]));
  fprintf(fp, "\n");
    
  writeEigenstates(fp, P, E);

  fprintf(fp, "</Multiconfprojected>\n");
  fclose(fp);

  return 0;
}


int readMultiMulticonfigfile(const char* fname,
			     char*** mbfilep,
			     Projection* P,
			     MultiSlaterDet** Qp,
			     Indices** Inp,
			     int* nstates,
			     Eigenstates* E)
{
  FILE* fp;
  char buf[BUFSIZE];
  int i, n;
  char msldetfname[255], md5msldetfname[33];

  if (!(fp = fopen(fname, "r"))) {
    fprintf(stderr, "couldn't open %s for reading\n", fname);
    return -1;
  }

  // get the Projection parameters
  do	
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Multiconfprojected ", 19) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "didn't find <Multiconfprojected ...>\n");
    return -1;
  }

  int dim;
  char projpar[40];
  sscanf(buf, "<Multiconfprojected %d %d %s>", &n, &dim, projpar);
  *nstates = n;

  // read the Slater determinants
  *mbfilep = (char**) malloc(n*sizeof(char *));
  for (i=0; i<n; i++)
    (*mbfilep)[i] = (char*) malloc(255*sizeof(char));
  *Qp = malloc(n*sizeof(MultiSlaterDet));
  *Inp = malloc(n*sizeof(Indices));
  char* fileIn = malloc(255*sizeof(char));
  for (i=0; i<n; i++) {
    do	
      fgets(buf, BUFSIZE, fp);
    while (strncmp(buf, "<MBFile ", 7) && !feof(fp));
    if (feof(fp)) {
      fprintf(stderr, "...   did't find <SlaterDetFile ...>\n");
      fclose(fp);
      return -1;
    }
    sscanf(buf, "<MBFile %s %s %s>", msldetfname, md5msldetfname, fileIn);
    extractIndicesfromString(&fileIn, &(*Inp)[i]);
    if (strncmp(md5msldetfname, md5hash(msldetfname), 32)) {
      fprintf(stderr, "...    MultiSlaterDetFile %s changed\n", msldetfname);
      return -1;
    }
    strcpy((*mbfilep)[i], msldetfname);
    if (readMultiSlaterDetfromFile(&(*Qp)[i], &(*Inp)[i], msldetfname))
      return -1;
  }
  
  // odd or even ?
  int odd = Qp[0]->A % 2;

  // initialize Projection
  initProjection(P, odd, projpar);

  // read the Eigenstates
  if (readEigenstates(fp, P, E, dim))
    return -1;

  return 0;
}

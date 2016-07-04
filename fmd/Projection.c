/**

  \file Projection.c

  angular momentum projection of matrix elements


  (c) 2003-2007 Thomas Neff

*/


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <zlib.h>

#include "SlaterDet.h"
#include "Observables.h"
#include "CenterofMass.h"

#include "numerics/zcw.h"
#include "numerics/wignerd.h"
#include "numerics/clebsch.h"
#include "numerics/cmat.h"

#include "misc/physics.h"
#include "misc/utils.h"

#include "Projection.h"
#include "Symmetry.h"


#define SQR(x) ((x)*(x))



// pi=0 : positive parity, pi=1 : negative parity
char* AngmomtoStr(int j, int pi)
{
  char* str = malloc(5*sizeof(char));

  if (j%2) 
    sprintf(str, "%d2%c", j, pi ? '-' : '+');
  else
    sprintf(str, "%d%c", j/2, pi ? '-' : '+');

  return str;
}


char* ProjectiontoStr(const Projection* P)
{
  char* suffix = malloc(36*sizeof(char));

  char jmaxsuf[10] = "";
  if (P->jmax != JMAX)
    sprintf(jmaxsuf, "jmax-%d-", P->jmax-1);

  char angsuf[12];
  if (P->ang == AngNone)
    sprintf(angsuf, "ang-none");
  else if (P->ang == AngProd)
    sprintf(angsuf, "ang-%d-%d", P->angprod.nbeta, 
	    P->angprod.nazimuth);
  else if (P->ang == AngProdA)
    sprintf(angsuf, "ang-a%d-%d", P->angprod.nbeta,
	    P->angprod.nazimuth);
  else if (P->ang == AngZCW)
    sprintf(angsuf, "ang-zcw-%d", P->angzcw.idx);

  char cmsuf[14] = "";
  if (P->cm == CMNone)
    sprintf(cmsuf, "-cm-none");
  else if (P->cm == CMSimple)
    ;
  else if (P->cm == CMProd)
    sprintf(cmsuf, "-cm-%d-%d-%d", P->cmprod.nr, P->cmprod.ntheta, 
	    P->cmprod.nphi);
  else if (P->cm == CMPoly) {
    if (P->cmpoly.npoly == 4) 
      sprintf(cmsuf, "-cm-%d-tet", P->cmpoly.nr);
    else if (P->cmpoly.npoly == 6) 
      sprintf(cmsuf, "-cm-%d-oct", P->cmpoly.nr);
    else if (P->cmpoly.npoly == 8) 
      sprintf(cmsuf, "-cm-%d-cbe", P->cmpoly.nr);
  }

  sprintf(suffix, "%s%s%s", jmaxsuf, angsuf, cmsuf);

  return suffix;
}


// ang-[none|nbeta-nazimuth|zcw-izcw][-cm-nr-[ntheta-nphi|tet|oct|cbe]]

// should return error code if could analyze parameters

int initProjection(Projection* P, int odd, const char* projpar)
{
  P->odd = odd;

  P->jmax = JMAX;
  P->ang = AngNone; P->cm = CMSimple;

  char projparcpy[strlen(projpar)];
  strcpy(projparcpy, projpar);
  char* c = projparcpy;

  // different jmax
  c = strtok(c, "-");
  if (!strncmp(c, "jmax", 4)) {
    c=strtok(NULL, "-");
    P->jmax = atoi(c)+1;
    c=strtok(NULL, "-");
  }

  // angular momentum projection
  if (!strncmp(c, "ang", 3)) {
    c=strtok(NULL, "-");
    if (!strncmp(c, "none", 4))
      P->ang = AngNone;
    else if (!strncmp(c, "zcw", 3)) {
      P->ang = AngZCW;
      c=strtok(NULL, "-");
      P->angzcw.idx = atoi(c);
    } else if (!strncmp(c, "a", 1)) {
      P->ang = AngProdA;
      P->angprod.nbeta = atoi(++c);
      c=strtok(NULL, "-");
      P->angprod.nazimuth = atoi(c);
    } else {
      P->ang = AngProd;
      P->angprod.nbeta = atoi(c);
      c=strtok(NULL, "-");
      P->angprod.nazimuth = atoi(c);
    }
  }
    
  // center of mass projection
  c = strtok(NULL, "-");
  if (c && !strncmp(c, "cm", 2)) {
    c=strtok(NULL, "-");
    if (!strncmp(c, "none", 4)) {
      P->cm = CMNone; }
    else {
      int nr = atoi(c);
      c=strtok(NULL, "-");
      if (!strncmp(c, "tet", 3)) {
	P->cm = CMPoly;
	P->cmpoly.nr = nr;
	P->cmpoly.npoly = 4;
      } else if (!strncmp(c, "oct", 3)) {
	P->cm = CMPoly;
	P->cmpoly.nr = nr;
	P->cmpoly.npoly = 6;
      } else if (!strncmp(c, "cbe", 3)) {
	P->cm = CMPoly;
	P->cmpoly.nr = nr;
	P->cmpoly.npoly = 8;
      } else {
	P->cm = CMProd;
	P->cmprod.nr = nr;
	P->cmprod.ntheta = atoi(c);
	c = strtok(NULL, "-");
	P->cmprod.nphi = atoi(c);
      }
    }
  }

  return 0;
}


void fprintProjectinfo(FILE* fp, const Projection* P)
{
  if (P->ang == AngNone)
    fprintf(fp, "# angular momentum projeciton - none\n");
  else if (P->ang == AngProd)
    fprintf(fp, "# angular momentum projection - product integration with %d,%d,%d points (%d angles)\n", 
	    P->angprod.nazimuth, 
	    P->angprod.nbeta, P->angprod.nazimuth,
	    SQR(P->angprod.nazimuth)*P->angprod.nbeta);
  else if (P->ang == AngProdA)
    fprintf(fp, "# angular momentum projection - adaptive product integration with %d,%d,%d points (%d angles)\n", 
	    P->angprod.nazimuth, 
	    P->angprod.nbeta, P->angprod.nazimuth,
	    SQR(P->angprod.nazimuth)*P->angprod.nbeta);
  else if (P->ang == AngZCW)
    fprintf(fp, "# angular momentum projection - integration using ZCW set %d (%d angles)\n",
	    P->angzcw.idx, nangles3(P->angzcw.idx));

  if (P->cm == CMNone)
    fprintf(fp, "# center of mass projection - none\n");
  else if (P->cm == CMSimple)
    fprintf(fp, "# center of mass projection - simple\n");
  else if (P->cm == CMProd)
    fprintf(fp, "# center of mass projection - product integration with %d,%d,%d points\n",
	    P->cmprod.nr, P->cmprod.ntheta, P->cmprod.nphi);
  else if (P->cm == CMPoly)
    fprintf(fp, "# center of mass projection - integration using %d %s points\n",
	    P->cmpoly.nr, (P->cmpoly.npoly == 4 ? "tetrahedral" : "octahedral"));
}


void* initprojectedvector(const Projection* P, int size, int n)
{
  int p,j;
  char (**ppv)[size];
  int jmax=P->jmax;

  ppv = malloc((jmax+1)*sizeof(void*));
  for (p=0; p<=1; p++)
    for (j=P->odd; j<jmax; j=j+2) {
      ppv[idxpij(jmax,p,j)] = malloc(n*(j+1)*size);
    }

  return ppv;
}

// don't really allocate space for matrix elements
void* initprojectedvectornull(const Projection* P, int size, int n)
{
  int p,j;
  char (**ppv)[size];
  int jmax=P->jmax;

  ppv = malloc((jmax+1)*sizeof(void*));
  for (p=0; p<=1; p++)
    for (j=P->odd; j<jmax; j=j+2) {
      ppv[idxpij(jmax,p,j)] = NULL;
    }

  return ppv;
}


void* initprojectedmatrix(const Projection* P, int size, int n)
{
  int p,j;
  char (**ppm)[size];
  int jmax=P->jmax;

  ppm = malloc((jmax+1)*sizeof(void*));
  for (p=0; p<=1; p++)
    for (j=P->odd; j<jmax; j=j+2)
      ppm[idxpij(jmax,p,j)] = malloc(SQR(n*(j+1))*size);

  return ppm;
}


void* initprojectedtransitionvector(const Projection* P, 
				    int rank, int pi, int size, 
				    int nfin, int nini)
{
  int jmax=P->jmax;
  int oddfin = (P->odd+rank)%2;
  int pini,jini, pfin,jfin;
  int ifin;
  char (****v)[size];

  v = malloc((jmax+1)*sizeof(void*));
  for (pfin=0; pfin<=1; pfin++)
    for (jfin=oddfin; jfin<jmax; jfin=jfin+2) {
      v[idxpij(jmax,pfin,jfin)] = malloc((jmax+1)*sizeof(void*));
      pini = (pfin+pi)%2;
      for (jini=abs(jfin-rank); jini<=min(jmax-1, jfin+rank); jini=jini+2) {
	v[idxpij(jmax,pfin,jfin)][idxpij(jmax,pini,jini)] = 
	  malloc(nfin*(jfin+1)*sizeof(void*));
	for (ifin=0; ifin<nfin*(jfin+1); ifin++)
	  v[idxpij(jmax,pfin,jfin)][idxpij(jmax,pini,jini)][ifin] = 
	    malloc(nini*(jini+1)*size);
      }
    }
  return v;
}

// don't really allocate space for matrix elements
void* initprojectedtransitionvectornull(const Projection* P, 
					int rank, int pi, int size, 
					int nfin, int nini)
{
  int jmax=P->jmax;
  int oddfin = (P->odd+rank)%2;
  int pini,jini, pfin,jfin;
  int ifin;
  char (****v)[size];

  v = malloc((jmax+1)*sizeof(void*));
  for (pfin=0; pfin<=1; pfin++)
    for (jfin=oddfin; jfin<jmax; jfin=jfin+2) {
      v[idxpij(jmax,pfin,jfin)] = malloc((jmax+1)*sizeof(void*));
      pini = (pfin+pi)%2;
      for (jini=abs(jfin-rank); jini<=min(jmax-1, jfin+rank); jini=jini+2) {
	v[idxpij(jmax,pfin,jfin)][idxpij(jmax,pini,jini)] = 
	  malloc(nfin*(jfin+1)*sizeof(void*));
	for (ifin=0; ifin<nfin*(jfin+1); ifin++)
	  v[idxpij(jmax,pfin,jfin)][idxpij(jmax,pini,jini)][ifin] = NULL; 
      }
    }
  return v;
}

/// store projected MEs between two ManyBody states
void* initprojectedMBME(const Projection* P, const ManyBodyOperator* Op)
{
  return initprojectedmatrix(P, (Op->rank+1)*Op->size*sizeof(complex double), 1);
}


void* initprojectedVector(const Projection* P, const ManyBodyOperator* Op, 
			  int n)
{
  return initprojectedvector(P, (Op->rank+1)*Op->size*sizeof(complex double), n);
}

void* initprojectedVectornull(const Projection* P, const ManyBodyOperator* Op, 
			      int n)
{
  return initprojectedvectornull(P, (Op->rank+1)*Op->size*sizeof(complex double), n);
}


void* initprojectedtransitionVector(const Projection* P, 
				    const ManyBodyOperator* Op, 
				    int nfin, int nini)
{
  return initprojectedtransitionvector(P, Op->rank, Op->pi, 
				       Op->size*sizeof(complex double), 
				       nfin, nini);
}

void* initprojectedtransitionVectornull(const Projection* P, 
				    const ManyBodyOperator* Op, 
				    int nfin, int nini)
{
  return initprojectedtransitionvectornull(P, Op->rank, Op->pi, 
					   Op->size*sizeof(complex double), 
					   nfin, nini);
}



int writeprojectedMBME(gzFile fp, const Projection* P, 
		       const ManyBodyOperator* Op, 
		       Symmetry S, Symmetry Sp, 
		       const void* mbme)
{
  int jmax = P->jmax;
  int odd = P->odd;
  int size = Op->size;
  int dim = Op->dim;
  int rank = Op->rank;
  complex double (**me)[(rank+1)*size] = mbme;

  int p, j, m, k, l, r;

  gzprintf(fp, "<ProjectedMBME %s %d %d>\n", Op->name, S, Sp);
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {

      gzprintf(fp, "<ProjectedME %d %d>\n", j, p); 
      for (m=-j; m<=j; m=m+2)
	for (k=-j; k<=j; k=k+2) {
	  // Sp fixes j, k
	  // if Op is scalar also j,m fixed
	  if ((Op->rank != 0 || SymmetryAllowed(S, p, j, m)) && 
	      SymmetryAllowed(Sp, p, j, k)) {
	    for (l=0; l<dim; l++)
	      for (r=0; r<=rank; r++)
		gzprintf(fp, "(%15.8e,%15.8e) ", 
			creal(me[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)]),
			cimag(me[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)]));
	    gzprintf(fp, "\n");
	  }
	}
	gzprintf(fp, "</ProjectedME>\n");

    }	
  gzprintf(fp, "</ProjectedMBME>\n");

  return 0;
}


int writeprojectedMBMEtoFile(const char* mbfilea, const char* mbfileb,
			     const Projection* P,
			     const ManyBodyOperator* Op,
			     Symmetry Sa, Symmetry Sb,
			     const void* me)
{
  gzFile mefp;
  char mefilename[255];

  ensuredir("ME");

  snprintf(mefilename, 255, "ME/%s--%s%s--%s%s--%s.gz", 
	   Op->name, 
	   (Sa==0 ? "" : strjoin(SymmetrytoStr(Sa), ":")), 
	   filepart(mbfilea), 
	   (Sb==0 ? "" : strjoin(SymmetrytoStr(Sb), ":")), 
	   filepart(mbfileb), 
	   ProjectiontoStr(P));

  if (!(mefp = gzopen(mefilename, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", mefilename);
    return -1;
  }
    
  gzprintinfo(mefp);

  gzprintf(mefp, "<MBFile %s %s %d>\n", 
	   mbfilea, md5hash(mbfilea), Sa);
  gzprintf(mefp, "<MBFile %s %s %d>\n\n", 
	   mbfileb, md5hash(mbfileb), Sb);

  writeprojectedMBME(mefp, P, Op, Sa, Sb, me);

  gzclose(mefp);

  return 0;
} 


#define BUFSIZE 65536

static char* buf;

int readprojectedMBME(gzFile fp, 
		      const Projection* P, const ManyBodyOperator* Op,
		      Symmetry S, Symmetry Sp,
		      void* mbme)
{
  int jmax = P->jmax;
  int odd = P->odd;
  int size = Op->size;
  int rank = Op->rank;
  int dim = Op->dim;
  complex double (**me)[(rank+1)*size] = mbme;

  int p, j, m, k, l, r;

  // possibly initialize space for buffer
  if (!buf)
    buf = malloc(BUFSIZE);

  char *c;
  double ref, imf;

  // initialize matrix elements to zero
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      for (k=-j; k<=j; k=k+2)
	for (m=-j; m<=j; m=m+2)
	  for (l=0; l<dim; l++)
	    for (r=0 ; r<=rank; r++)
	      me[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = 0.0;
    }	

  do
    gzgets(fp, buf, BUFSIZE);
  while (strncmp(buf, "<ProjectedMBME ", 15) && !gzeof(fp));
  if (gzeof(fp)) {
    fprintf(stderr, "did't find <ProjectedMBME ...>\n");
    return -1;
  }
  char opname[80];
  int fileS, fileSp;
  sscanf(buf, "<ProjectedMBME %s %d %d>", opname, &fileS, &fileSp);
  if (strcmp(stripstr(opname, ">"), Op->name)) {
    fprintf(stderr, "not %s matrixelements\n", Op->name);
    return -1;
  }
  if (fileS != S || fileSp != Sp) {
    fprintf(stderr, "Symmetries don't match\n");
    return -1;
  }

  int filej, filep;
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      gzgets(fp, buf, BUFSIZE);
      sscanf(buf, "<ProjectedME %d %d>", &filej, &filep);
      if (filej != j || filep != p) {
	fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
	return -1;
      }

      for (m=-j; m<=j; m=m+2)
	for (k=-j; k<=j; k=k+2) {
	  if ((Op->rank != 0 || SymmetryAllowed(S, p, j, m)) && 
	      SymmetryAllowed(Sp, p, j, k)) {
	    gzgets(fp, buf, BUFSIZE);
	    c = strtok(buf, " ,()");
	    for (l=0; l<dim; l++)
	      for (r=0; r<=rank; r++) {
		ref=atof(c); c=strtok(NULL, " ,()");
		imf=atof(c); c=strtok(NULL, " ,()");
		me[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = ref+I*imf;
	    }
	  }
	}

      gzgets(fp, buf, BUFSIZE);
      if (strncmp(buf, "</ProjectedME>", 14)) {
	fprintf(stderr, "didn't find </ProjectedME>\n");
	return -1;
      }
    }
  gzgets(fp, buf, BUFSIZE);
  if (strncmp(buf, "</ProjectedMBME>", 16)) {
    fprintf(stderr, "didn't find </ProjectedMBME>\n");
    return -1;
  }

  return 0;
}




int readprojectedMBMEfromFile(const char* mbfilea, const char* mbfileb, 
			      const Projection* P,
			      const ManyBodyOperator* Op,
			      Symmetry Sa, Symmetry Sb,
			      void* mbme)
{
  gzFile mefp;
  char mefilename[255];
  char buf[BUFSIZE];
  char fnam[255], md5fnam[33];
  int fileS;

  snprintf(mefilename, 255, "ME/%s--%s%s--%s%s--%s.gz", 
	   Op->name, 
	   (Sa==0 ? "" : strjoin(SymmetrytoStr(Sa), ":")), 
	   filepart(mbfilea), 
	   (Sb==0 ? "" : strjoin(SymmetrytoStr(Sb), ":")), 
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
  do
    gzgets(mefp, buf, BUFSIZE);
  while (strncmp(buf, "<MBFile ", 7) && !gzeof(mefp));
  if (gzeof(mefp)) {
    fprintf(stderr, "...   did't find <MBFile ...>\n");
    gzclose(mefp);
    return -2;
  }
  sscanf(buf, "<MBFile %s %s %d>", fnam, md5fnam, &fileS);
  if (strncmp(md5fnam, md5hash(mbfilea), 32)) {
    fprintf(stderr, "...    MBFile does not match with existing ME\n");
    return -2;
  }
  if (fileS != Sa) {
    fprintf(stderr, "...    Symmetry does not match with existing ME\n");
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
  sscanf(buf, "<MBFile %s %s %d>", fnam, md5fnam, &fileS);
  if (strncmp(md5fnam, md5hash(mbfileb), 32)) {
    fprintf(stderr, "...    MBFile does not match with existing ME\n");
    gzclose(mefp);
    return -2;
  }
  if (fileS != Sb) {
    fprintf(stderr, "...    Symmetry does not match with existing ME\n");
    return -2;
  }

  // now read the matrix elements
  int res;
  res = readprojectedMBME(mefp, P, Op, Sa, Sb, mbme); 

  gzclose(mefp);
  return res;
}

// 07/13/09 important change: do not normalize Q and Qp anymore

void calcprojectedMBME(const Projection* P, const ManyBodyOperator* Op,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       Symmetry S, Symmetry Sp,
		       void* mbme)
{
  int size=Op->size;
  int rank=Op->rank;
  int dim=Op->dim;
  complex double (**val)[(rank+1)*size] = mbme;

  SlaterDet Qpp;
  SlaterDetAux X;

  int jmax = P->jmax;
  int odd = P->odd;

  // norms of SlaterDets
  initSlaterDetAux(Q, &X);

  // calcSlaterDetAuxod(Q, Q, &X);
  // double norm = sqrt(creal(X.ovlap));

  // initSlaterDetAux(Qp, &X);
  // calcSlaterDetAuxod(Qp, Qp, &X);
  // double normp = sqrt(creal(X.ovlap));  

  // set up cm integration
  cmintegrationpara cmpara;
  initcmintegration(P, Q, Qp, &cmpara);
  int ncm = cmpara.n;

  // set up angular momentum integration
  angintegrationpara angpara;
  initangintegration(P, Q, Qp, S, Sp, &angpara);
  int nang = angpara.n;

  int l, r;
  int p, j, m, k;

  // set matrix elements to zero
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      for (k=-j; k<=j; k=k+2)
	for (m=-j; m<=j; m=m+2)
	  for (l=0; l<dim; l++)
	    for (r=0; r<=rank; r++)
	      val[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = 0.0;
    }	

  initSlaterDet(Qp, &Qpp);
  
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
      
      copySlaterDet(Qp, &Qpp);
      moveSlaterDet(&Qpp, xcm);
      rotateSlaterDet(&Qpp, alpha, beta, gamma);

      for (ip=0; ip<=1; ip++) {
	  if (ip) invertSlaterDet(&Qpp);

	  // can only calculate Auxilliaries if Sldets are compatible
	  if (Q->A == Qp->A) {
            if (Q->Z == Qp->Z && Q->N == Qp->N)
              calcSlaterDetAuxod(Q, &Qpp, &X);
            else
              calcSlaterDetAuxodsingular(Q, &Qpp, &X);
          }
	  Op->me(Op->par, Q, &Qpp, &X, sval);

	  complex double w;
	  for (p=0; p<=1; p++)
	    for (j=odd; j<jmax; j=j+2)
	      for (k=-j; k<=j; k=k+2)
		for (m=-j; m<=j; m=m+2) {
		  if ((Op->rank != 0 || SymmetryAllowed(S, p, j, m)) &&
		      SymmetryAllowed(Sp, p, j, k)) {
		    w = weight * (p && ip%2 ? -1 : 1)*
		      (j+1)/(8*SQR(M_PI))*Djmkstar(j,m,k,alpha,beta,gamma);
		    for (l=0; l<dim; l++)
		      for (r=0; r<=rank; r++)
			val[idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] += 
			  w*sval[r+l*(rank+1)];
		  }	
		}	
      }

    }

  }
  freeSlaterDetAux(&X);
  freeSlaterDet(&Qpp);
}	


void calcprojectedMBMEs(const Projection* P, const ManyBodyOperators* Ops,
			const SlaterDet* Q, const SlaterDet* Qp,
			Symmetry S, Symmetry Sp,
			void* mbme)
{
  int size=Ops->size;
  int dim=Ops->dim;
  complex double ***val = mbme;

  SlaterDet Qpp;
  SlaterDetAux  X;

  int jmax = P->jmax;
  int odd = P->odd;

  // norms of SlaterDets
  initSlaterDetAux(Q, &X);
  // calcSlaterDetAuxod(Q, Q, &X);
  // double norm = sqrt(creal(X.ovlap));

  // initSlaterDetAux(Qp, &X);
  // calcSlaterDetAuxod(Qp, Qp, &X);
  // double normp = sqrt(creal(X.ovlap));  

  // set up cm integration
  cmintegrationpara cmpara;
  initcmintegration(P, Q, Qp, &cmpara);
  int ncm = cmpara.n;

  // set up angular momentum integration
  angintegrationpara angpara;
  initangintegration(P, Q, Qp, S, Sp, &angpara);
  int nang = angpara.n;

  int l, r;
  int o, p, j, m, k;

  // helpful for indexing matrix elements

  int ranko[Ops->n];
  for (o=0; o<Ops->n; o++)
    ranko[o] = Ops->Op[o].rank;

  int sizeo[Ops->n];
  for (o=0; o<Ops->n; o++)
    sizeo[o] = size*(ranko[o]+1);

  int io[Ops->n]; io[0] = 0;
  for (o=0; o<Ops->n-1; o++)
    io[o+1] = io[o]+sizeo[o];

  int no=0;
  for (o=0; o<Ops->n; o++)
    no += sizeo[o];


  // set matrix elements to zero
  for (o=0; o<Ops->n; o++)
    for (p=0; p<=1; p++)
      for (j=odd; j<jmax; j=j+2) {
	for (k=-j; k<=j; k=k+2)
	  for (m=-j; m<=j; m=m+2)
	    for (l=0; l<dim; l++)
	      for (r=0; r<=ranko[o]; r++)
		val[o][idxpij(jmax,p,j)][r+l*(ranko[o]+1)+idxjmk(j,m,k)*sizeo[o]] = 0.0;
    }	

  initSlaterDet(Qp, &Qpp);
  
  int icm; 
  double xcm[3]; double weightcm;

  int iang;
  double alpha, beta, gamma; double weightang;


  complex double sval[no];
  double weight;
  int ip;

  for (icm=0; icm<ncm; icm++) {
    getcmintegrationpoint(icm, &cmpara, xcm, &weightcm);

    for (iang=0; iang<nang; iang++) {
      getangintegrationpoint(iang, &angpara, &alpha, &beta, &gamma, &weightang);
      // weight = 1.0/(2*norm*normp)*weightcm*weightang;
      weight = 0.5*weightcm*weightang;
      
      copySlaterDet(Qp, &Qpp);
      moveSlaterDet(&Qpp, xcm);
      rotateSlaterDet(&Qpp, alpha, beta, gamma);

      for (ip=0; ip<=1; ip++) {
	  if (ip) invertSlaterDet(&Qpp);

	  // can only calculate Auxilliaries if Sldets are compatible
	  if (Q->A == Qp->A) {
            if (Q->Z == Qp->Z && Q->N == Qp->N)
              calcSlaterDetAuxod(Q, &Qpp, &X);
            else
              calcSlaterDetAuxodsingular(Q, &Qpp, &X);
          }
	  Ops->me(Ops->par, Q, &Qpp, &X, sval);

	  complex double w;
	  for (o=0; o<Ops->n; o++)
	    for (p=0; p<=1; p++)
	      for (j=odd; j<jmax; j=j+2)
		for (k=-j; k<=j; k=k+2)
		  for (m=-j; m<=j; m=m+2) {
		    if ((ranko[o] != 0 || SymmetryAllowed(S, p, j, m)) &&
			SymmetryAllowed(Sp, p, j, k)) {
		      w = weight * (p && ip%2 ? -1 : 1)*
			(j+1)/(8*SQR(M_PI))*Djmkstar(j,m,k,alpha,beta,gamma);
		      for (l=0; l<dim; l++)
			for (r=0; r<=ranko[o]; r++)
			  val[o][idxpij(jmax,p,j)][r+l*(ranko[o]+1)+idxjmk(j,m,k)*sizeo[o]] += 
			    w*sval[r+l*(ranko[o]+1)+io[o]];
		  }	
		}	
      }

    }

  }	
  freeSlaterDetAux(&X);
  freeSlaterDet(&Qpp);
}	


void hermitizeprojectedMBME(const Projection* P, const ManyBodyOperator* Op,
			    void* mbme, int n)
{
  int size=Op->size;
  int rank=Op->rank;
  int dim=Op->dim;
  complex double (***val)[(rank+1)*size] = mbme;

  int jmax = P->jmax;
  int odd = P->odd;

  int p, j, m, k;
  int a,b;
  int l, r;

  complex double valu, vall;

  // hermitize matrix elements
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2)
      for (b=0; b<n; b++)
	for (k=-j; k<=j; k=k+2)
	  for (a=0; a<n; a++)
	    for (m=-j; m<=j; m=m+2)
	      for (l=0; l<dim; l++) 
		for (r=0; r<=rank+1; r++) {
		  valu = val[a+b*n][idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)];
		  vall = val[b+a*n][idxpij(jmax,p,j)][idxjmk(j,k,m)][r+l*(rank+1)];
		  val[a+b*n][idxpij(jmax,p,j)][idxjmk(j,m,k)][r+l*(rank+1)] = 0.5*(valu+conj(vall));
		  val[b+a*n][idxpij(jmax,p,j)][idxjmk(j,k,m)][r+l*(rank+1)] = 0.5*(conj(valu)+vall);
	      }
}    

// calculates reduced matrix element divided by sqrt(2j+1)

void calcexpectprojectedMBME(const Projection* P,
			     const ManyBodyOperator* Op,
			     const void* mbme,
			     const Symmetry* S,
			     const Eigenstates* E,
			     void* expectmbme)
{
  int p,j,i;
  int k, k1, k2;
  int l, nu;
  int a, b;
  int n=E->n;
  int idx;
  int nj;
  int rank=Op->rank;
  int size=Op->size;
  int dim=Op->dim;
  int odd=P->odd;
  int jmax=P->jmax;
  complex double (***me)[(rank+1)*size] = mbme;
  complex double (**expme)[size] = expectmbme;
 
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {

      idx = idxpij(jmax,p,j);
      nj=n*(j+1);
      for (i=0; i<E->dim[idx]; i++) {

	for (l=0; l<dim; l++)
	  expme[idx][i][l] = 0.0;

	for (b=0; b<n; b++)
	  for (k2=-j; k2<=j; k2=k2+2)
	    for (a=0; a<n; a++)
	      for (k1=-j; k1<=j; k1=k1+2)	    
		for (k=max(-j,k1-rank); k<=min(j,k1+rank); k=k+2) {
		  nu=k1-k;
		  if (SymmetryAllowed(S[a], p, j, k1) &&
		      SymmetryAllowed(S[b], p, j, k2)) {
		    for (l=0; l<dim; l++)
		      expme[idx][i][l] +=
			clebsch(j, rank, j, k, nu, k1)*
			conj(E->V[idx][idxnjm(a,j,k1) + i*nj])*
			me[a+b*n][idx][idxjmk(j,k,k2)][(nu+rank)/2+l*(rank+1)]*
			E->V[idx][idxnjm(b,j,k2) + i*nj];
		  }
		}
      }		
    }		
}

// calculate the expectation value only for a selected state
void calcexpectprojectedMBMEipj(const Projection* P,
				const ManyBodyOperator* Op,
				const void* mbme,
				const Symmetry* S,
				const Eigenstates* E,
				int j, int p, int i,
				void* expectmbme)
{
  int k, k1, k2;
  int l, nu;
  int a, b;
  int n=E->n;
  int idx;
  int ii;
  int nj;
  int rank=Op->rank;
  int size=Op->size;
  int dim=Op->dim;
  int odd=P->odd;
  int jmax=P->jmax;
  complex double (***me)[(rank+1)*size] = mbme;
  complex double (**expme)[size] = expectmbme;
 
  idx = idxpij(jmax,p,j);
  ii = E->index[idx][i];
  nj=n*(j+1);

  // allocate space for expectation values
  expme[idx] = malloc(nj*size*sizeof(complex double));

  for (l=0; l<dim; l++)
    expme[idx][ii][l] = 0.0;

  for (b=0; b<n; b++)
    for (k2=-j; k2<=j; k2=k2+2)
      for (a=0; a<n; a++)
	for (k1=-j; k1<=j; k1=k1+2)	    
	  for (k=max(-j,k1-rank); k<=min(j,k1+rank); k=k+2) {
	    nu=k1-k;
	    if (SymmetryAllowed(S[a], p, j, k1) &&
		SymmetryAllowed(S[b], p, j, k2)) {
	      for (l=0; l<dim; l++)
		expme[idx][ii][l] +=
		  clebsch(j, rank, j, k, nu, k1)*
		  conj(E->V[idx][idxnjm(a,j,k1) + ii*nj])*
		  me[a+b*n][idx][idxjmk(j,k,k2)][(nu+rank)/2+l*(rank+1)]*
		  E->V[idx][idxnjm(b,j,k2) + ii*nj];
	    }
	  }
}


// calculate transition matrix elements for Eigenstates Efin <- Eini

// calculates reduced matrix element divided by sqrt(2jfin+1)

// slow but transparent
/*
void calctransitionprojectedMBME(const Projection* P,
				 const ManyBodyOperator* Op,
				 const void* mbme,
				 const Symmetry* Sfin,
				 const Symmetry* Sini,
				 const Eigenstates* Efin,
				 const Eigenstates* Eini,
				 void* transmbme)
{	
  int pini,jini,iini, pfin,jfin,ifin;
  int k, kfin, kini;
  int l, nu;
  int afin, aini;
  int nfin=Efin->n;
  int nini=Eini->n;
  int njini, njfin;
  int idxini, idxfin;
  int p=Op->pi;
  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  int oddini=P->odd;
  int jmax=P->jmax;
  complex double (***me)[(rank+1)*size] = mbme;
  complex double (****transme)[size] = transmbme;

  for (pini=0; pini<=1; pini++)
    for (jini=oddini; jini<jmax; jini=jini+2) {
      
      idxini = idxpij(jmax,pini,jini);
      njini=nini*(jini+1);

      pfin=(pini+p)%2;
      for (jfin=abs(jini-rank); jfin<=min(jmax-1,jini+rank); jfin=jfin+2) {

	idxfin = idxpij(jmax,pfin,jfin);
	njfin=nfin*(jfin+1);
	for (iini=0; iini<Eini->dim[idxini]; iini++)
	  for (ifin=0; ifin<Efin->dim[idxfin]; ifin++) {

	    for (l=0; l<dim; l++)
	      transme[idxfin][idxini][ifin][iini][l] = 0.0;

	    for (aini=0; aini<nini; aini++)
	      for (kini=-jini; kini<=jini; kini=kini+2)
		for (afin=0; afin<nfin; afin++)
		  for (kfin=-jfin; kfin<=jfin; kfin=kfin+2)
		    for (k=max(-jini,kfin-rank); k<=min(jini,kfin+rank); k=k+2) {
		      nu = kfin-k;
		      if (SymmetryAllowed(Sfin[afin], pfin, jfin, kfin) &&
			  SymmetryAllowed(Sini[aini], pini, jini, kini)) {
			for (l=0; l<dim; l++)
			  transme[idxfin][idxini][ifin][iini][l] +=
			    clebsch(jini, rank, jfin, k, nu, kfin)*
			    conj(Efin->V[idxfin][idxnjm(afin,jfin,kfin) + ifin*njfin])*
			    me[afin+aini*nfin][idxini][idxjmk(jini,k,kini)][(nu+rank)/2+l*(rank+1)]*

			    Eini->V[idxini][idxnjm(aini,jini,kini) + iini*njini];
		      }
		    }
	  }
      }			
    }
}	
*/


// much faster, calculate first product of matrix elements with initial vectors
// then product with final vectors in a second step
void calctransitionprojectedMBME(const Projection* P,
				 const ManyBodyOperator* Op,
				 const void* mbme,
				 const Symmetry* Sfin,
				 const Symmetry* Sini,
				 const Eigenstates* Efin,
				 const Eigenstates* Eini,
				 void* transmbme)
{	
  int pini,jini,iini, pfin,jfin,ifin;
  int k, kfin, kini;
  int l, nu;
  int afin, aini;
  int nfin=Efin->n;
  int nini=Eini->n;
  int njini, njfin;
  int idxini, idxfin;
  int p=Op->pi;
  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  int oddini=P->odd;
  int jmax=P->jmax;
  complex double (***me)[(rank+1)*size] = mbme;
  complex double (****transme)[size] = transmbme;

  complex double* meVwork = malloc(nfin*nini*(jmax+1)*(jmax+1)*size*(rank+1)*sizeof(complex double));

  for (pini=0; pini<=1; pini++)
    for (jini=oddini; jini<jmax; jini=jini+2) {
      
      idxini = idxpij(jmax,pini,jini);
      njini=nini*(jini+1);

      pfin=(pini+p)%2;
      for (jfin=abs(jini-rank); jfin<=min(jmax-1,jini+rank); jfin=jfin+2) {

	idxfin = idxpij(jmax,pfin,jfin);
	njfin=nfin*(jfin+1);

	complex double (*meV)[nini*(jini+1)][(jini+1)][size*(rank+1)] = meVwork;

	for (iini=0; iini<Eini->dim[idxini]; iini++)
	  for (afin=0; afin<nfin; afin++)
	    for (k=-jini; k<=jini; k=k+2)
	      for (l=0; l<dim; l++)
		for (nu=-rank; nu<=rank; nu=nu+2)
		  meV[afin][iini][idxjm(jini,k)][(nu+rank)/2+l*(rank+1)] = 0.0;

	for (iini=0; iini<Eini->dim[idxini]; iini++)
	  for (afin=0; afin<nfin; afin++)
	    for (k=-jini; k<=jini; k=k+2)
	      for (aini=0; aini<nini; aini++)
		for (kini=-jini; kini<=jini; kini=kini+2)
		  if (SymmetryAllowed(Sini[aini], pini, jini, kini)) {
	      
		    for (l=0; l<dim; l++)
		      for (nu=-rank; nu<=rank; nu=nu+2)
			meV[afin][iini][idxjm(jini,k)][(nu+rank)/2+l*(rank+1)] +=
			  me[afin+aini*nfin][idxini][idxjmk(jini,k,kini)][(nu+rank)/2+l*(rank+1)]*

			  Eini->V[idxini][idxnjm(aini,jini,kini) + iini*njini];
		  }


	for (iini=0; iini<Eini->dim[idxini]; iini++)
	  for (ifin=0; ifin<Efin->dim[idxfin]; ifin++) {

	    for (l=0; l<dim; l++)
	      transme[idxfin][idxini][ifin][iini][l] = 0.0;

	    for (afin=0; afin<nfin; afin++)
	      for (kfin=-jfin; kfin<=jfin; kfin=kfin+2)
		for (k=max(-jini,kfin-rank); k<=min(jini,kfin+rank); k=k+2) {
		  nu = kfin-k;
		  if (SymmetryAllowed(Sfin[afin], pfin, jfin, kfin)) {
		    for (l=0; l<dim; l++)
		      transme[idxfin][idxini][ifin][iini][l] +=
			clebsch(jini, rank, jfin, k, nu, kfin)*
			conj(Efin->V[idxfin][idxnjm(afin,jfin,kfin) + ifin*njfin])*
			meV[afin][iini][idxjm(jini,k)][(nu+rank)/2+l*(rank+1)];
		  }
		}
	  }	
      }			
    }
  
  free(meVwork);
}	


void calctransitionprojectedMBMEipj(const Projection* P,
				    const ManyBodyOperator* Op,
				    const void* mbme,
				    const Symmetry* Sfin,
				    const Symmetry* Sini,
				    const Eigenstates* Efin,
				    const Eigenstates* Eini,
				    int jfin, int pfin, int ifin,
				    int jini, int pini, int iini,
				    void* transmbme)
{	
  int k, kfin, kini;
  int l, nu;
  int afin, aini;
  int nfin=Efin->n;
  int nini=Eini->n;
  int njini, njfin;
  int idxini, idxfin;
  int iiini, iifin;
  int p=Op->pi;
  int size=Op->size;
  int dim=Op->dim;
  int rank=Op->rank;
  int oddini=P->odd;
  int jmax=P->jmax;
  complex double (***me)[(rank+1)*size] = mbme;
  complex double (****transme)[size] = transmbme;

  idxini = idxpij(jmax,pini,jini);
  iiini = Eini->index[idxini][iini];
  njini=nini*(jini+1);

  idxfin = idxpij(jmax,pfin,jfin);
  iifin = Efin->index[idxfin][ifin];
  njfin=nfin*(jfin+1);

  // allocate space for transition matrix elements
  transme[idxfin][idxini][iifin] = malloc(njini*size*sizeof(complex double));

  for (l=0; l<dim; l++)
    transme[idxfin][idxini][iifin][iiini][l] = 0.0;

  for (aini=0; aini<nini; aini++)
    for (kini=-jini; kini<=jini; kini=kini+2)
      for (afin=0; afin<nfin; afin++)
	for (kfin=-jfin; kfin<=jfin; kfin=kfin+2)
	  for (k=max(-jini,kfin-rank); k<=min(jini,kfin+rank); k=k+2) {
	    nu = kfin-k;
	    if (SymmetryAllowed(Sfin[afin], pfin, jfin, kfin) &&
		SymmetryAllowed(Sini[aini], pini, jini, kini)) {
	      for (l=0; l<dim; l++)
		transme[idxfin][idxini][iifin][iiini][l] +=
		  clebsch(jini, rank, jfin, k, nu, kfin)*
		  conj(Efin->V[idxfin][idxnjm(afin,jfin,kfin) + iifin*njfin])*
		  me[afin+aini*nfin][idxini][idxjmk(jini,k,kini)][(nu+rank)/2+l*(rank+1)]*
		  Eini->V[idxini][idxnjm(aini,jini,kini) + iiini*njini];
	    }
	  }
}	


void initEigenstates(const Projection* P, Eigenstates* E, int n)
{
  int jmax = P->jmax;

  E->n = n;

  E->v = initprojectedvector(P, sizeof(complex double), n);
  E->V = initprojectedmatrix(P, sizeof(complex double), n);
  E->norm = initprojectedvector(P, sizeof(complex double), n);
  E->dim = malloc((jmax+1)*sizeof(int));
  E->ngood = malloc((jmax+1)*sizeof(int));
  E->index = initprojectedvector(P, sizeof(int), n);
}


void initAmplitudes(const Projection* P, Amplitudes* A, int n)
{
  A->n = n;
  A->ngood = initprojectedvector(P, sizeof(int), n);
  A->amp = initprojectedmatrix(P, sizeof(complex double), n);
}


void calcEigenstates(const Projection* P,
		     const Interaction* Int,
		     const Observablesod ***obsme, 
		     Eigenstates* E, double thresh)
{
  int n = E->n;
  int odd = P->odd;
  int jmax = P->jmax;

  complex double* H = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* N = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* N2 = malloc(SQR(n*(jmax+1))*sizeof(complex double));

  int a, b;
  int p, j, m, k, i;
  int nj, ipj;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {

      nj = n*(j+1);
      ipj = idxpij(jmax,p,j);

      for (b=0; b<n; b++)
	for (k=-j; k<=j; k=k+2)
	  for (a=0; a<n; a++)
	    for (m=-j; m<=j; m=m+2) {
	      N[idxnjm(a,j,m) + idxnjm(b,j,k)*nj] = 
		obsme[a+b*n][ipj][idxjmk(j,m,k)].n;

	      H[idxnjm(a,j,m) + idxnjm(b,j,k)*nj] = 
		obsme[a+b*n][ipj][idxjmk(j,m,k)].h;
	    }

      generalizedeigensystem(H, N, nj, thresh, 
			     E->v[ipj], E->V[ipj], 
			     &E->dim[ipj]);

      // normalize eigenvectors such that norm = \sum_M |<Q|Q;JMalpha>|^2

      multcmat(N, N, N2, nj);

      for (i=0; i<E->dim[ipj]; i++) {
        
        double norm2=0.0, ovl2 = 0.0;
	for (b=0; b<nj; b++)
	  for (a=0; a<nj; a++) {
            norm2 += conj(E->V[ipj][a+i*nj])*N[a+b*nj]*E->V[ipj][b+i*nj];
            ovl2 += conj(E->V[ipj][a+i*nj])*N2[a+b*nj]*E->V[ipj][b+i*nj];
          }
        
        // sort out unphysical (due to numerics) states
        double scale;
        if (norm2 < 0.0 || ovl2 < 0.0 || isnan(norm2) || isnan(ovl2))
          scale = 0.0;
        else
          scale = sqrt(ovl2)/norm2;

        for (a=0; a<nj; a++)
          E->V[ipj][a+i*nj] *= scale;

	E->norm[ipj][i] = 0.0;
	for (b=0; b<nj; b++)
	  for (a=0; a<nj; a++)
	    E->norm[ipj][i] += 
	      conj(E->V[ipj][a+i*nj])*N[a+b*nj]*E->V[ipj][b+i*nj];
      }

   }

  free(H); free(N); free(N2);
}


void calcEigenstatesK(const Projection* P,
		      const Interaction* Int,
		      const Observablesod ***obsme, 
		      Eigenstates* E, int K, double thresh)
{
  int n = E->n;
  int odd = P->odd;
  int jmax = P->jmax;

  complex double* H = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* N = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  
  int a, b;
  int p, j, m, k, i;
  int nj, ipj;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      nj = n*(j+1);
      ipj = idxpij(jmax,p,j);

      for (b=0; b<n; b++)
	for (k=-j; k<=j; k=k+2)
	  for (a=0; a<n; a++)
	    for (m=-j; m<=j; m=m+2) {
	      N[idxnjm(a,j,m) + idxnjm(b,j,k)*nj] =  
		(m == K && k == K) ? obsme[a+b*n][ipj][idxjmk(j,m,k)].n : 0.0;

	      H[idxnjm(a,j,m) + idxnjm(b,j,k)*nj] = 
		(m == K && k == K) ? obsme[a+b*n][ipj][idxjmk(j,m,k)].h : 0.0;
	    }

      generalizedeigensystem(H, N, nj, thresh, 
			     E->v[ipj], E->V[ipj], 
			     &E->dim[ipj]);

      for (i=0; i<E->dim[ipj]; i++) {
	E->norm[ipj][i] = 0.0;
	for (b=0; b<nj; b++)
	  for (a=0; a<nj; a++)
	    E->norm[ipj][i] += 
	      conj(E->V[ipj][a+i*nj])*N[a+b*nj]*E->V[ipj][b+i*nj];
      }

   }

  free(H); free(N);
}

/*
void calcMultiEigenstates(const Projection* P,
			  const Interaction* Int,
			  const Observablesod ***obsme,
			  const Eigenstates* Ep,
			  Eigenstates* multiE, Amplitudes* multiA, 
			  double thresh)
{
  int n=multiE->n;
  int odd=P->odd;
  int jmax=P->jmax;

  complex double* H = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* N = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* v = malloc(n*(jmax+1)*sizeof(complex double));
  complex double* V = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  
  int a, b;
  int p, j, ipj, i, m, k;
  int dim, d;
  int ai, bi, iai, ibi, idxa, idxb;
  complex double normi2, norma2;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      ipj=idxpij(jmax,p,j);
      
      dim=0;
      for (a=0; a<n; a++)
	dim += Ep[a].ngood[ipj];

      // do we have at least a one-dimensional space
      if (dim ==0) {
	multiE->dim[ipj] = 0;
      } else {

      idxb=-1;
      for (b=0; b<n; b++)
	for (bi=0; bi<Ep[b].ngood[ipj]; bi++) {
	  ibi=Ep[b].index[ipj][bi];
	  idxb++;
	  idxa=-1;
	  for (a=0; a<n; a++)	
	    for (ai=0; ai<Ep[a].ngood[ipj]; ai++) {
	      iai=Ep[a].index[ipj][ai];
	      idxa++;
	      N[idxa+idxb*dim] = 0.0;
	      H[idxa+idxb*dim] = 0.0;

	      for (k=-j; k<=j; k=k+2)
		for (m=-j; m<=j; m=m+2) {
		  N[idxa+idxb*dim] +=
		    conj(Ep[a].V[ipj][idxjm(j,m)+iai*(j+1)])*
		    obsme[a+b*n][ipj][idxjmk(j,m,k)].n*
		    Ep[b].V[ipj][idxjm(j,k)+ibi*(j+1)];
			 
		  H[idxa+idxb*dim] +=
		    conj(Ep[a].V[ipj][idxjm(j,m)+iai*(j+1)])*
		    obsme[a+b*n][ipj][idxjmk(j,m,k)].h*
		    Ep[b].V[ipj][idxjm(j,k)+ibi*(j+1)];

		}
	    }	
	}      

      generalizedeigensystem(H, N, dim, thresh, 
			     v, V, 
			     &d);

      // embed solution into full space

      multiE->dim[ipj] = d;
      int nj=n*(j+1);

      for (i=0; i<d; i++) {
	multiE->v[ipj][i] = v[i];

	for (k=-j; k<=j; k=k+2) {
	  idxa = -1;
	  for (a=0; a<n; a++) {
	    multiE->V[ipj][idxnjm(a,j,k)+i*nj] = 0.0;

	    for (ai=0; ai<Ep[a].ngood[ipj]; ai++) {
	      idxa++;
	      iai = Ep[a].index[ipj][ai];

	      multiE->V[ipj][idxnjm(a,j,k)+i*nj] += 
                Ep[a].V[ipj][idxjm(j,k)+iai*(j+1)]*
		  V[idxa+i*dim];
	    }
	  }				      
	}
      }

      // calculate Amplitudes

      for (i=0; i<d; i++) {

	normi2 = 0.0;
	for (idxb=0; idxb<dim; idxb++)
	  for (idxa=0; idxa<dim; idxa++)
	    normi2 += conj(V[idxa+i*dim])*N[idxa+idxb*dim]*V[idxb+i*dim];

	multiE->norm[ipj][i] = normi2;

	idxa = -1;
	for (a=0; a<n; a++) {
	  multiA->ngood[ipj][a] = Ep[a].ngood[ipj];
	  
	  for (ai=0; ai<Ep[a].ngood[ipj]; ai++) {
	    idxa++;
	    norma2 = N[idxa+idxa*dim];
	    multiA->amp[ipj][ai+a*(j+1)+i*n*(j+1)] = 0.0;
	    for (idxb=0; idxb<dim; idxb++)
	      multiA->amp[ipj][ai+a*(j+1)+i*n*(j+1)] +=
		N[idxa+idxb*dim]*V[idxb+i*dim]/csqrt(norma2*normi2);
	  }
	}
      }
      }

   }

  free(H); free(N);
  free(v); free(V);
}
*/

// use basis states |Q^i;JMalpha> normalized to 1

void calcMultiEigenstates(const Projection* P,
			  const Interaction* Int,
			  const Observablesod ***obsme,
			  const Eigenstates* Ep,
			  Eigenstates* multiE, Amplitudes* multiA, 
			  double thresh)
{
  int n=multiE->n;
  int odd=P->odd;
  int jmax=P->jmax;

  complex double* H = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* N = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  complex double* v = malloc(n*(jmax+1)*sizeof(complex double));
  complex double* V = malloc(SQR(n*(jmax+1))*sizeof(complex double));
  
  int a, b;
  int p, j, ipj, i, m, k;
  int dim, d;
  int ai, bi, iai, ibi, idxa, idxb;
  double norma2, normb2, normi2;        

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      ipj=idxpij(jmax,p,j);

      dim=0;
      for (a=0; a<n; a++)
	dim += Ep[a].ngood[ipj];

      // do we have at least a one-dimensional space
      if (dim ==0) {
	multiE->dim[ipj] = 0;
      } else {

      idxb=-1;
      for (b=0; b<n; b++)
	for (bi=0; bi<Ep[b].ngood[ipj]; bi++) {
	  ibi=Ep[b].index[ipj][bi];
          normb2=Ep[b].norm[ipj][ibi];
	  idxb++;
	  idxa=-1;
	  for (a=0; a<n; a++)	
	    for (ai=0; ai<Ep[a].ngood[ipj]; ai++) {
	      iai=Ep[a].index[ipj][ai];
              norma2=Ep[a].norm[ipj][iai];
	      idxa++;
	      N[idxa+idxb*dim] = 0.0;
	      H[idxa+idxb*dim] = 0.0;

	      for (k=-j; k<=j; k=k+2)
		for (m=-j; m<=j; m=m+2) {
		  N[idxa+idxb*dim] +=
		    conj(Ep[a].V[ipj][idxjm(j,m)+iai*(j+1)])*
		    obsme[a+b*n][ipj][idxjmk(j,m,k)].n*
		    Ep[b].V[ipj][idxjm(j,k)+ibi*(j+1)]/
                    sqrt(norma2*normb2);
			 
		  H[idxa+idxb*dim] +=
		    conj(Ep[a].V[ipj][idxjm(j,m)+iai*(j+1)])*
		    obsme[a+b*n][ipj][idxjmk(j,m,k)].h*
		    Ep[b].V[ipj][idxjm(j,k)+ibi*(j+1)]/
                    sqrt(norma2*normb2);

		}
	    }	
	}      

      /* debugging

      fprintf(stderr, "\nJ^pi = %s\n\n", AngmomtoStr(j, p));
      
      fprintf(stderr, "N matrix\n");
      fprintcmat(stderr, dim, N);

      fprintf(stderr, "H matrix\n");
      fprintcmat(stderr, dim, H);

      */

      generalizedeigensystem(H, N, dim, thresh, 
			     v, V, 
			     &d);

      // embed solution into full space

      multiE->dim[ipj] = d;
      int nj=n*(j+1);

      for (i=0; i<d; i++) {
	multiE->v[ipj][i] = v[i];

	for (k=-j; k<=j; k=k+2) {
	  idxa = -1;
	  for (a=0; a<n; a++) {
	    multiE->V[ipj][idxnjm(a,j,k)+i*nj] = 0.0;

	    for (ai=0; ai<Ep[a].ngood[ipj]; ai++) {
	      idxa++;
	      iai = Ep[a].index[ipj][ai];
              norma2 = Ep[a].norm[ipj][iai];

	      multiE->V[ipj][idxnjm(a,j,k)+i*nj] += 
                Ep[a].V[ipj][idxjm(j,k)+iai*(j+1)]/sqrt(norma2)*
		  V[idxa+i*dim];
	    }
	  }				      
	}
      }

      // calculate Amplitudes

      for (i=0; i<d; i++) {

	normi2 = 0.0;
	for (idxb=0; idxb<dim; idxb++)
	  for (idxa=0; idxa<dim; idxa++)
	    normi2 += conj(V[idxa+i*dim])*N[idxa+idxb*dim]*V[idxb+i*dim];

	multiE->norm[ipj][i] = normi2;

	idxa = -1;
	for (a=0; a<n; a++) {
	  multiA->ngood[ipj][a] = Ep[a].ngood[ipj];
	  
	  for (ai=0; ai<Ep[a].ngood[ipj]; ai++) {
	    idxa++;
	    norma2 = N[idxa+idxa*dim];

	    multiA->amp[ipj][ai+a*(j+1)+i*n*(j+1)] = 0.0;
	    for (idxb=0; idxb<dim; idxb++)
	      multiA->amp[ipj][ai+a*(j+1)+i*n*(j+1)] +=
		N[idxa+idxb*dim]*V[idxb+i*dim]/sqrt(norma2*normi2);
	  }
	}
      }

      }
   }

  free(H); free(N);
  free(v); free(V);
}


int writeEigenstates(FILE* fp, 
		     const Projection* P, const Eigenstates* E)
{
  int p, j;
  int i;
  int odd=P->odd;
  int jmax=P->jmax;
  int n=E->n;
  int idx;
  int nj, dimj;

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      idx = idxpij(jmax,p,j);
      nj = n*(j+1);
      dimj = E->dim[idx];

      fprintf(fp, "<Eigenstates %d %d>\n", j, p);
      fprintf(fp, "<Indices %2d %2d>\n", dimj, E->ngood[idx]);
      if (dimj > 0) {
	for (i=0; i<dimj; i++)
	  fprintf(fp, "%2d   ", E->index[idx][i]);
	fprintf(fp, "\n</Indices>\n");
	fprintf(fp, "<Norms>\n");
	for (i=0; i<dimj; i++)
	  fprintf(fp, "(%10.8f,%10.8f)   ",
		  creal(E->norm[idx][i]),
		  cimag(E->norm[idx][i]));
	fprintf(fp, "\n</Norms>\n");
	fprintf(fp, "<Energies>\n");
	for (i=0; i<dimj; i++)
	  fprintf(fp, "(%8.5f,%8.5f)   ", 
		  creal(E->v[idx][i]),
		  cimag(E->v[idx][i]));
	fprintf(fp, "\n</Energies>\n");
	fprintf(fp, "<Eigenvectors>\n");
	fprintcmatcols(fp, nj, dimj, E->V[idx]);
	fprintf(fp, "</Eigenvectors>\n");
      }
      fprintf(fp, "</Eigenstates>\n");
    }
  
  return 0;
}


#define BUFSIZE 65536

int readEigenstates(FILE* fp,
		    const Projection* P, Eigenstates* E, int n)
{
  int odd=P->odd;
  int jmax=P->jmax;
  char buf[BUFSIZE];
  int p, j;
  int filep, filej, filedim, filegood;
  int idx;
  int nj;
  int dimj;

  initEigenstates(P, E, n);

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      idx = idxpij(jmax,p,j);
      nj = n*(j+1);

      do 
	fgets(buf, BUFSIZE, fp);
      while (strncmp(buf, "<Eigenstates ", 13) && !feof(fp));
      if (feof(fp)) {
	fprintf(stderr, "didn't find <Eigenstates ...>\n");
	return -1;
      }
      sscanf(buf, "<Eigenstates %d %d>", &filej, &filep);
      if (j != filej || p != filep) {
	fprintf(stderr, "expected <Eigenstates %d %d>\n", j, p);
	return -1;
      }
      fgets(buf, BUFSIZE, fp);
      sscanf(buf, "<Indices %d %d>", &filedim, &filegood);
      E->dim[idx] = filedim;
      E->ngood[idx] = filegood;
      if (filedim > 0) {
	dimj = filedim;
	freadivec(fp, dimj, E->index[idx]);
	fgets(buf, BUFSIZE, fp); // </Indices>
	fgets(buf, BUFSIZE, fp); // <Norms>
	freadcvec(fp, dimj, E->norm[idx]);
	fgets(buf, BUFSIZE, fp); // </Norms>
	fgets(buf, BUFSIZE, fp); // <Energies>
	freadcvec(fp, dimj, E->v[idx]);
	fgets(buf, BUFSIZE, fp); // </Energies>
	fgets(buf, BUFSIZE, fp); // <Eigenvectors>
	freadcmatcols(fp, nj, dimj, E->V[idx]);
	fgets(buf, BUFSIZE, fp); // </Eigenvectors>
      }
      fgets(buf, BUFSIZE, fp);
      if (strncmp(buf, "</Eigenstates>", 14)) {
	fprintf(stderr, "didn't find </Eigenstates>\n");
	return -1;
      }
    }
  return 0;
}



static int cmpmerit(void* ap, void* bp) 
{
  double a=*(double*) ap;
  double b=*(double*) bp;
  
  return (a<=b ? (a<b ? 1 : 0) : -1);
}

#define MAXSTATES 15
#define ERRJMAX 0.25            // |<J^2> - j(j+1)| should be smaller than ERRJMAX

void sortEigenstates(const Projection* P,
		     const Interaction* Int,
		     const Observablesod** obs,
		     Eigenstates* E, 
		     double minnorm, int all)
{
  int odd=P->odd;
  int jmax=P->jmax;
  int n=E->n;

  int p,j,i;
  int dimj;
  int nj;
  int ngood;

  double maxnorm;

  struct merit {
    double val;
    int idx;
  } merits[n*(jmax+1)];

  maxnorm = 0.0;
  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2)
      for (i=0; i<E->dim[idxpij(jmax,p,j)]; i++)
	maxnorm = fmax(maxnorm, creal(obs[idxpij(jmax,p,j)][i].n));

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      nj = n*(j+1);
      dimj = E->dim[idxpij(jmax,p,j)];

      // determine merits of eigenstates
      for (i=0; i<dimj; i++) {
	merits[i].idx = i;

        // norm too small ?
	if (creal(obs[idxpij(jmax,p,j)][i].n) < minnorm*maxnorm)
	  merits[i].val = -100000.0;	
        // J2 looks suspicious ?
        else if (cabs(obs[idxpij(jmax,p,j)][i].j2/obs[idxpij(jmax,p,j)][i].n - 0.25*j*(j+2)) > ERRJMAX)
          merits[i].val = -100000.0;
        // else use the energy
	else
	  merits[i].val = -creal(obs[idxpij(jmax,p,j)][i].h/obs[idxpij(jmax,p,j)][i].n);
      }
	
      // sort according to merits value
      qsort(merits, dimj, sizeof(struct merit), cmpmerit);

      i=0; while (i<dimj && ( all || merits[i].val > -100000.0)) i++;
      ngood = min(i,MAXSTATES);
      E->ngood[idxpij(jmax,p,j)] = ngood;
      for (i=0; i<dimj; i++)
	E->index[idxpij(jmax,p,j)][i] = merits[i].idx;

    }
}


// don't check whether MBFile still matches 
int readEigenstatesfromFile(const char* fname,
			    const Projection* P, Eigenstates* E, int n)
{
  FILE* fp;
  char fullname[255];

  snprintf(fullname, 255, "%s.%s.states", fname, ProjectiontoStr(P));

  if (!(fp = fopen(fullname, "r"))) {  
    fprintf(stderr, "couldn't open %s for reading\n", fullname);
    return -1;
  }

  fprintf(stderr, "... reading Eigenstates from file %s\n", fullname);

  int err=readEigenstates(fp, P, E, n);
  fclose(fp);

  return err;
}


int writeMulticonfigfile(const char* fname,
			 const Projection* P,
			 Symmetry* S,
			 const char** mbfile,
			 const Eigenstates* E, int n)
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
  fprintf(fp, "<Multiconfprojected %d %s>\n", n, ProjectiontoStr(P));

  int i;
  for (i=0; i<n; i++)
    fprintf(fp, "<MBFile %s %s %d>\n",
	    mbfile[i], md5hash(mbfile[i]), S[i]);
  fprintf(fp, "\n");
    
  writeEigenstates(fp, P, E);

  fprintf(fp, "</Multiconfprojected>\n");
  fclose(fp);

  return 0;
}


#define BUFSIZE 1024
int readMulticonfigfile(const char* fname,
			char*** slaterdetfilep,
			Projection* P,
			SlaterDet** Qp,
			Symmetry** Sp,
			Eigenstates* E,
			int* nstates)
{
  FILE* fp;
  char buf[BUFSIZE];
  int i, n;
  char sldetfname[255], md5sldetfname[33];

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
  char projpar[40];
  sscanf(buf, "<Multiconfprojected %d %s>", &n, projpar);
  *nstates = n;

  // read the Slater determinants
  *slaterdetfilep = (char**) malloc(n*sizeof(char *));
  for (i=0; i<n; i++)
    (*slaterdetfilep)[i] = (char*) malloc(255*sizeof(char));
  *Qp = malloc(n*sizeof(SlaterDet));
  *Sp = malloc(n*sizeof(Symmetry));
  int fileS;
  for (i=0; i<n; i++) {
    do	
      fgets(buf, BUFSIZE, fp);
    while (strncmp(buf, "<MBFile ", 7) && !feof(fp));
    if (feof(fp)) {
      fprintf(stderr, "...   did't find <SlaterDetFile ...>\n");
      fclose(fp);
      return -1;
    }
    sscanf(buf, "<MBFile %s %s %d>", sldetfname, md5sldetfname, &fileS);
    if (strncmp(md5sldetfname, md5hash(sldetfname), 32)) {
      fprintf(stderr, "...    SlaterDetFile %s changed\n", sldetfname);
      return -1;
    }
    strcpy((*slaterdetfilep)[i], sldetfname);
    if (readSlaterDetfromFile(&(*Qp)[i], sldetfname))
      return -1;
    (*Sp)[i] = fileS;
  }
  
  // odd or even ?
  int odd = Qp[0]->A % 2;

  // initialize Projection
  initProjection(P, odd, projpar);

  // read the Eigenstates
  if (readEigenstates(fp, P, E, n))
    return -1;

  return 0;
}

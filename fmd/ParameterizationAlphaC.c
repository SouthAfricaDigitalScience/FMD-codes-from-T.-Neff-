/**

  \file ParameterizationAlphaC.c

  Parametrization of SlaterDet.

  alpha-cluster nuclei with fixed width a and vanishing momenta


  (c) 2003 Thomas Neff

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>


#include "Gaussian.h"
#include "SlaterDet.h"
#include "gradGaussian.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"
#include "ParameterizationAlphaC.h"

#include "misc/utils.h"


Parameterization ParameterizationAlphaC = {
  name : "AlphaC",
  Pararead : AlphaCread,
  Parawrite : AlphaCwrite,
  Paraclone : AlphaCclone,
  ParainitSlaterDet : AlphaCinitSlaterDet,
  ParatoSlaterDet : AlphaCtoSlaterDet,
  ParaprojectgradSlaterDet : AlphaCprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 3		


void AlphaCclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = (double*) malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;

  // internals stores width a
  qp->internals = malloc(sizeof(double));
  double* qpa = qp->internals;
  double* qa = q->internals;
  *qpa = *qa;
}


#define BUFSIZE 255

int AlphaCread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;


  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "AlphaC")) {
    fprintf(stderr, "not a AlphaC parameter set\n");
    return -1;
  }
  strcpy(q->name, stripstr(qname, ">"));

  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d %d %d", &q->A, &q->Z, &q->N) != 3) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
  }
  if (q->A != 2*q->Z || q->A != 2*q->N) {
    fprintf(stderr, "not an alpha nucleus in line\n>>%s<<\n",
	    stripstr(buf, "\n"));
    return -1;
  }
  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d", &q->ngauss) != 1) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
  }

  q->n = q->A/4*NGAUSS;

  q->x = (double*) malloc(q->n*sizeof(double));

  double a, b0re, b1re, b2re; 

  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%lf", &a) != 1) {      
    fprintf(stderr, "malformed line\n>>%s<<\n", 
	    stripstr(buf,"\n"));
    return -1;
  }
  q->internals = malloc(sizeof(double));
  double* qa = q->internals;
  *qa = a;

  for (i=0; i<q->A/4; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "%lf %lf %lf", 
	       &b0re, &b1re, &b2re) != 3) {
      fprintf(stderr, "malformed line\n>>%s<<\n", 
	      stripstr(buf,"\n"));
      return -1;
    }
    q->x[i*NGAUSS+ 0] = b0re;
    q->x[i*NGAUSS+ 1] = b1re;
    q->x[i*NGAUSS+ 2] = b2re;
  }     

  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</Para>", 7)) {
    fprintf(stderr, "didn't find </Para>\n");
    return -1;
  }	

  return 0;
}


int AlphaCwrite(FILE* fp, const Para* q)
{
  int i;

  fprintf(fp, "<Para AlphaC %s>\n", q->name);
  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);
  double *qa = q->internals;
  fprintf(fp, "%12.8f\n", *qa);
  for (i=0; i<q->A/4; i++)
    fprintf(fp,
	    "%12.8f  %12.8f  %12.8f\n",
	    q->x[i*NGAUSS+ 0], q->x[i*NGAUSS+ 1], q->x[i*NGAUSS+ 2]);
  fprintf(fp, "</Para>\n");

  return 0;
}


void AlphaCinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k;

  Q->A = q->A; Q->Z = q->Z; Q->N = q->N;
  Q->ngauss = q->ngauss;

  Q->idx = (int*) malloc(q->A*sizeof(int));
  for (k=0; k<q->A; k++)
    Q->idx[k] = k;
    
  Q->ng = (int*) malloc(q->A*sizeof(int));
  for (k=0; k<q->A; k++)
    Q->ng[k] = 1;

  Q->G = (Gaussian*) malloc(q->ngauss*sizeof(Gaussian));
  for (k=0; k<q->A/4; k++) {
    Q->G[4*k+0].xi = 1; 
    Q->G[4*k+1].xi = 1; 
    Q->G[4*k+2].xi = -1; 
    Q->G[4*k+3].xi = -1; 
  }
}


void AlphaCtoSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  double* qa = q->internals;

  for (k=0; k<q->A/4; k++)
    for (ki=0; ki<4; ki++) {
    Q->G[4*k+ki].chi[0] = ki % 2 ? 1 : 0;
    Q->G[4*k+ki].chi[1] = ki % 2 ? 0 : 1;
    Q->G[4*k+ki].a      = *qa;
    Q->G[4*k+ki].b[0]   = q->x[k*NGAUSS+ 0];
    Q->G[4*k+ki].b[1]   = q->x[k*NGAUSS+ 1];
    Q->G[4*k+ki].b[2]   = q->x[k*NGAUSS+ 2];
  }
}


void AlphaCprojectgradSlaterDet(const Para* q,
				const gradSlaterDet* dQ, double* dq)
{
  int k, ki, i;

  // only converting a complex in a real gradient

  for (k=0; k<q->A/4; k++) {
    for (i=0; i<NGAUSS; i++)
      dq[k*NGAUSS+ i] = 0.0;
    for (ki=0; ki<4; ki++) {
      dq[k*NGAUSS+ 0] += 2.0*creal(dQ->gradval[4*k+ki].b[0]);
      dq[k*NGAUSS+ 1] += 2.0*creal(dQ->gradval[4*k+ki].b[1]);
      dq[k*NGAUSS+ 2] += 2.0*creal(dQ->gradval[4*k+ki].b[2]);
    }
  }
}

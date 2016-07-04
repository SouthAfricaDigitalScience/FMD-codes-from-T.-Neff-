/**

  \file ParameterizationAlphaA.c

  Parametrization of SlaterDet.

  alpha-cluster nuclei with fixed width a


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
#include "ParameterizationAlphaA.h"

#include "misc/utils.h"


Parameterization ParameterizationAlphaA = {
  name : "AlphaA",
  Pararead : AlphaAread,
  Parawrite : AlphaAwrite,
  Paraclone : AlphaAclone,
  ParainitSlaterDet : AlphaAinitSlaterDet,
  ParatoSlaterDet : AlphaAtoSlaterDet,
  ParaprojectgradSlaterDet : AlphaAprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 6		


void AlphaAclone(const Para* q, Para* qp)
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

int AlphaAread(FILE* fp, Para* q)
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
  if (strcmp(Pname, "AlphaA")) {
    fprintf(stderr, "not a AlphaA parameter set\n");
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

  double a, b0re, b0im, b1re, b1im, b2re, b2im; 

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
    if (sscanf(buf, "(%lf,%lf) (%lf,%lf) (%lf,%lf)", 
	       &b0re, &b0im, &b1re, &b1im, &b2re, &b2im) != 6) {
      fprintf(stderr, "malformed line\n>>%s<<\n", 
	      stripstr(buf,"\n"));
      return -1;
    }
    q->x[i*NGAUSS+ 0] = b0re; q->x[i*NGAUSS+ 1] = b0im;
    q->x[i*NGAUSS+ 2] = b1re; q->x[i*NGAUSS+ 3] = b1im;
    q->x[i*NGAUSS+ 4] = b2re; q->x[i*NGAUSS+ 5] = b2im;
  }     

  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</Para>", 7)) {
    fprintf(stderr, "didn't find </Para>\n");
    return -1;
  }	

  return 0;
}


int AlphaAwrite(FILE* fp, const Para* q)
{
  int i;

  fprintf(fp, "<Para AlphaA %s>\n", q->name);
  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);
  double *qa = q->internals;
  fprintf(fp, "%12.8f\n", *qa);
  for (i=0; i<q->A/4; i++)
    fprintf(fp,
	    "(%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	    q->x[i*NGAUSS+ 0], q->x[i*NGAUSS+ 1], 
	    q->x[i*NGAUSS+ 2], q->x[i*NGAUSS+ 3],
	    q->x[i*NGAUSS+ 4], q->x[i*NGAUSS+ 5]);
  fprintf(fp, "</Para>\n");

  return 0;
}


void AlphaAinitSlaterDet(const Para* q, SlaterDet* Q)
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


void AlphaAtoSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  double* qa = q->internals;

  for (k=0; k<q->A/4; k++)
    for (ki=0; ki<4; ki++) {
    Q->G[4*k+ki].chi[0] = ki % 2 ? 1 : 0;
    Q->G[4*k+ki].chi[1] = ki % 2 ? 0 : 1;
    Q->G[4*k+ki].a      = *qa;
    Q->G[4*k+ki].b[0]   = q->x[k*NGAUSS+ 0]+I*q->x[k*NGAUSS+ 1];
    Q->G[4*k+ki].b[1]   = q->x[k*NGAUSS+ 2]+I*q->x[k*NGAUSS+ 3];
    Q->G[4*k+ki].b[2]   = q->x[k*NGAUSS+ 4]+I*q->x[k*NGAUSS+ 5];
  }
}


void AlphaAprojectgradSlaterDet(const Para* q,
				const gradSlaterDet* dQ, double* dq)
{
  int k, ki, i;

  // only converting a complex in a real gradient

  for (k=0; k<q->A/4; k++) {
    for (i=0; i<NGAUSS; i++)
      dq[k*NGAUSS+ i] = 0.0;
    for (ki=0; ki<4; ki++) {
      dq[k*NGAUSS+ 0] += 2.0*creal(dQ->gradval[4*k+ki].b[0]);
      dq[k*NGAUSS+ 1] += 2.0*cimag(dQ->gradval[4*k+ki].b[0]);
      dq[k*NGAUSS+ 2] += 2.0*creal(dQ->gradval[4*k+ki].b[1]);
      dq[k*NGAUSS+ 3] += 2.0*cimag(dQ->gradval[4*k+ki].b[1]);
      dq[k*NGAUSS+ 4] += 2.0*creal(dQ->gradval[4*k+ki].b[2]);
      dq[k*NGAUSS+ 5] += 2.0*cimag(dQ->gradval[4*k+ki].b[2]);
    }
  }
}

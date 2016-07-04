/**

  \file ParameterizationFMD.c

  Parametrization of SlaterDet.

  One-to-one mapping between Parameterization and SlaterDet


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
#include "ParameterizationFMD.h"

#include "misc/utils.h"
#include "misc/physics.h"


typedef struct {
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
} FMDInternals;


Parameterization ParameterizationFMD = {
  name : "FMD",
  Pararead : FMDread,
  Parawrite : FMDwrite,
  Paraclone : FMDclone,
  ParainitSlaterDet : FMDinitSlaterDet,
  ParatoSlaterDet : FMDtoSlaterDet,
  ParaprojectgradSlaterDet : FMDprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 12		


void FMDclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(FMDInternals));
  FMDInternals *qpint = qp->internals;
  FMDInternals *qint = q->internals;

  qpint->ng = malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->ng[i] = qint->ng[i];
  qpint->idx = malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->idx[i] = qint->idx[i];
  qpint->xi = malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->xi[i] = qint->xi[i];

}


#define BUFSIZE 255

int FMDread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;
  FMDInternals* qint;

  qint = q->internals = (FMDInternals*) malloc(sizeof(FMDInternals));
					

  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "FMD")) {
    fprintf(stderr, "not a FMD parameter set\n");
    return -1;
  }
  strcpy(q->name, stripstr(qname, ">"));

  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d %d %d", &q->A, &q->Z, &q->N) != 3) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    exit(-1);
  }
  if (q->A != q->Z+q->N) {
    fprintf(stderr, "proton and neutron number don't sum up in line\n>>%s<<\n",
	    stripstr(buf, "\n"));
    return -1;
  }
  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d", &q->ngauss) != 1) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
  }

  qint->idx = malloc(q->A*sizeof(int));
  qint->ng = malloc(q->A*sizeof(int));
  qint->xi = malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++) {
    qint->idx[i] = 0;
    qint->ng[i] = 0;
    qint->xi[i] = 0;
  }
  
  q->n = q->ngauss*NGAUSS;

  q->x = malloc(q->n*sizeof(double));

  int k, xi;
  double chi0re, chi0im, chi1re, chi1im;
  double are, aim, b0re, b0im, b1re, b1im, b2re, b2im; 
  for (i=0; i<q->ngauss; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "%d %d "
	       "(%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf)", 
	       &k, &xi, &chi0re, &chi0im, &chi1re, &chi1im, 
	       &are, &aim, &b0re, &b0im, &b1re, &b1im, &b2re, &b2im) != 14) {
      fprintf(stderr, "malformed line\n>>%s<<\n", 
	      stripstr(buf,"\n"));
      return -1;
    }
    qint->ng[k]++;
    qint->xi[k] = xi;
    q->x[i*NGAUSS+ 0] = chi0re; q->x[i*NGAUSS+ 1] = chi0im;
    q->x[i*NGAUSS+ 2] = chi1re; q->x[i*NGAUSS+ 3] = chi1im;
    q->x[i*NGAUSS+ 4] = are;    q->x[i*NGAUSS+ 5] = aim;
    q->x[i*NGAUSS+ 6] = b0re;   q->x[i*NGAUSS+ 7] = b0im;
    q->x[i*NGAUSS+ 8] = b1re;   q->x[i*NGAUSS+ 9] = b1im;
    q->x[i*NGAUSS+10] = b2re;   q->x[i*NGAUSS+11] = b2im;
  }     
  for (k=1; k<q->A; k++) 
    qint->idx[k] = qint->idx[k-1]+qint->ng[k-1];

  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</Para>", 7)) {
    fprintf(stderr, "didn't find </Para>\n");
    return -1;
  }

  return 0;
}


int FMDwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  FMDInternals* qint = q->internals;

  fprintf(fp, "<Para FMD %s>\n", q->name);
  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);
  for (k=0; k<q->A; k++)
    for (ki=0; ki<qint->ng[k]; ki++) {
      i = qint->idx[k]+ki;
      fprintf(fp,
	      "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	      k, qint->xi[k],
	      q->x[i*NGAUSS+ 0], q->x[i*NGAUSS+ 1],
	      q->x[i*NGAUSS+ 2], q->x[i*NGAUSS+ 3],
	      q->x[i*NGAUSS+ 4], q->x[i*NGAUSS+ 5],
	      q->x[i*NGAUSS+ 6], q->x[i*NGAUSS+ 7],
	      q->x[i*NGAUSS+ 8], q->x[i*NGAUSS+ 9],
	      q->x[i*NGAUSS+10], q->x[i*NGAUSS+11]);
  }
  fprintf(fp, "</Para>\n");

  return 0;
}


void FMDinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  FMDInternals* qint = q->internals;

  Q->A = q->A; Q->Z = q->Z; Q->N = q->N;
  Q->ngauss = q->ngauss;

  Q->idx = (int*) malloc(q->A*sizeof(int));
  for (k=0; k<q->A; k++)
    Q->idx[k] = qint->idx[k];

  Q->ng = (int*) malloc(q->A*sizeof(int));
  for (k=0; k<q->A; k++)
    Q->ng[k] = qint->ng[k];

  Q->G = (Gaussian*) malloc(q->ngauss*sizeof(Gaussian));
  for (k=0; k<q->A; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->idx[k]+ki].xi = qint->xi[k]; 
}


void FMDtoSlaterDet(const Para* q, SlaterDet* Q)
{
  int i;

  // this is a one-to-one mapping

  for (i=0; i<q->ngauss; i++) {
    Q->G[i].chi[0] = q->x[i*NGAUSS+ 0] + I* q->x[i*NGAUSS+ 1];
    Q->G[i].chi[1] = q->x[i*NGAUSS+ 2] + I* q->x[i*NGAUSS+ 3];
    Q->G[i].a      = q->x[i*NGAUSS+ 4] + I* q->x[i*NGAUSS+ 5];
    Q->G[i].b[0]   = q->x[i*NGAUSS+ 6] + I* q->x[i*NGAUSS+ 7];
    Q->G[i].b[1]   = q->x[i*NGAUSS+ 8] + I* q->x[i*NGAUSS+ 9];
    Q->G[i].b[2]   = q->x[i*NGAUSS+10] + I* q->x[i*NGAUSS+11];
  }
}


void FMDprojectgradSlaterDet(const Para* q,
			     const gradSlaterDet* dQ, double* dq)
{
  int i;

  // only converting a complex in a real gradient

  for (i=0; i<q->ngauss; i++) {
    dq[i*NGAUSS+ 0] = 2.0*creal(dQ->gradval[i].chi[0]);
    dq[i*NGAUSS+ 1] = 2.0*cimag(dQ->gradval[i].chi[0]);
    dq[i*NGAUSS+ 2] = 2.0*creal(dQ->gradval[i].chi[1]);
    dq[i*NGAUSS+ 3] = 2.0*cimag(dQ->gradval[i].chi[1]);
    dq[i*NGAUSS+ 4] = 2.0*creal(dQ->gradval[i].a);
    dq[i*NGAUSS+ 5] = 2.0*cimag(dQ->gradval[i].a);
    dq[i*NGAUSS+ 6] = 2.0*creal(dQ->gradval[i].b[0]);
    dq[i*NGAUSS+ 7] = 2.0*cimag(dQ->gradval[i].b[0]);
    dq[i*NGAUSS+ 8] = 2.0*creal(dQ->gradval[i].b[1]);
    dq[i*NGAUSS+ 9] = 2.0*cimag(dQ->gradval[i].b[1]);
    dq[i*NGAUSS+10] = 2.0*creal(dQ->gradval[i].b[2]);
    dq[i*NGAUSS+11] = 2.0*cimag(dQ->gradval[i].b[2]);
  }
}


void SlaterDetinitFMD(const SlaterDet* Q, Para* q)
{
  int i;
  FMDInternals* qint;

  strcpy(q->name, nucleusname(Q->A, Q->Z));

  q->A = Q->A;
  q->Z = Q->Z;
  q->N = Q->N;

  qint = q->internals = malloc(sizeof(FMDInternals));
					
  qint->idx = (int*) malloc(Q->A*sizeof(int));
  qint->ng = (int*) malloc(Q->A*sizeof(int));
  qint->xi = (int*) malloc(Q->A*sizeof(int));
  for (i=0; i<Q->A; i++) {
    qint->idx[i] = Q->idx[i];
    qint->ng[i] = Q->ng[i];
    qint->xi[i] = Q->G[Q->idx[i]].xi;
  }

  q->ngauss = Q->ngauss;
  q->n = Q->ngauss*NGAUSS;
  q->x = (double*) malloc(q->n*sizeof(double));

  for (i=0; i<Q->ngauss; i++) {
    q->x[i*NGAUSS+ 0] = creal(Q->G[i].chi[0]);
    q->x[i*NGAUSS+ 1] = cimag(Q->G[i].chi[0]);
    q->x[i*NGAUSS+ 2] = creal(Q->G[i].chi[1]);
    q->x[i*NGAUSS+ 3] = cimag(Q->G[i].chi[1]);
    q->x[i*NGAUSS+ 4] = creal(Q->G[i].a);
    q->x[i*NGAUSS+ 5] = cimag(Q->G[i].a);
    q->x[i*NGAUSS+ 6] = creal(Q->G[i].b[0]);
    q->x[i*NGAUSS+ 7] = cimag(Q->G[i].b[0]);
    q->x[i*NGAUSS+ 8] = creal(Q->G[i].b[1]);
    q->x[i*NGAUSS+ 9] = cimag(Q->G[i].b[1]);
    q->x[i*NGAUSS+10] = creal(Q->G[i].b[2]);
    q->x[i*NGAUSS+11] = cimag(Q->G[i].b[2]);
  }
}

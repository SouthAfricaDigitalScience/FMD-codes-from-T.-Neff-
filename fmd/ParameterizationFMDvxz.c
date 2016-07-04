/**

  \file ParameterizationFMDvxz.c

  Parametrization of SlaterDet symmetric under reflection on xz-plane


  (c) 2009 Thomas Neff

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>


#include "Gaussian.h"
#include "SlaterDet.h"
#include "gradGaussian.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"
#include "ParameterizationFMDvxz.h"

#include "misc/utils.h"
#include "misc/physics.h"


typedef struct {
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
} FMDvxzInternals;


Parameterization ParameterizationFMDvxz = {
  name : "FMDvxz",
  Pararead : FMDvxzread,
  Parawrite : FMDvxzwrite,
  Paraclone : FMDvxzclone,
  ParainitSlaterDet : FMDvxzinitSlaterDet,
  ParatoSlaterDet : FMDvxztoSlaterDet,
  ParaprojectgradSlaterDet : FMDvxzprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 12		


void FMDvxzclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = (double*) malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(FMDvxzInternals));
  FMDvxzInternals *qpint = qp->internals;
  FMDvxzInternals *qint = q->internals;

  qpint->ng = malloc(q->A/2*sizeof(int));
  for (i=0; i<q->A/2; i++)
    qpint->ng[i] = qint->ng[i];
  qpint->idx = malloc(q->A/2*sizeof(int));
  for (i=0; i<q->A/2; i++)
    qpint->idx[i] = qint->idx[i];
  qpint->xi = malloc(q->A/2*sizeof(int));
  for (i=0; i<q->A/2; i++)
    qpint->xi[i] = qint->xi[i];

}


#define BUFSIZE 255

int FMDvxzread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;
  FMDvxzInternals* qint;

  qint = q->internals = malloc(sizeof(FMDvxzInternals));
					
  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "FMDvxz")) {
    fprintf(stderr, "not a FMDvxz parameter set\n");
    return -1;
  }
  strcpy(q->name, stripstr(qname, ">"));

  // read nucleon parameters
  fgets(buf, BUFSIZE, fp);
  if (sscanf(buf, "%d %d %d", &q->A, &q->Z, &q->N) != 3) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
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

  qint->idx = (int*) malloc(q->A/2*sizeof(int));
  qint->ng = (int*) malloc(q->A/2*sizeof(int));
  qint->xi = (int*) malloc(q->A/2*sizeof(int));
  for (i=0; i<q->A/2; i++) {
    qint->idx[i] = 0;
    qint->ng[i] = 0;
    qint->xi[i] = 0;
  }
  
  q->n = q->ngauss/2*NGAUSS;

  q->x = (double*) malloc(q->n*sizeof(double));

  int k, xi;
  double chi0re, chi0im, chi1re, chi1im;
  double are, aim, b0re, b0im, b1re, b1im, b2re, b2im; 
  for (i=0; i<q->ngauss/2; i++) {
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
  for (k=1; k<q->A/2; k++) 
    qint->idx[k] = qint->idx[k-1]+qint->ng[k-1];

  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</Para>", 7)) {
    fprintf(stderr, "didn't find </Para>\n");
    return -1;
  }	

  return 0;
}


int FMDvxzwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  FMDvxzInternals* qint = q->internals;

  fprintf(fp, "<Para FMDvxz %s>\n", q->name);

  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);

  for (k=0; k<q->A/2; k++)
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


void FMDvxzinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  FMDvxzInternals* qint = q->internals;

  Q->A = q->A; 
  Q->Z = q->Z; 
  Q->N = q->N;
  Q->ngauss = q->ngauss;

  Q->idx = (int*) malloc(Q->A*sizeof(int));
  for (k=0; k<q->A/2; k++)
    Q->idx[k] = qint->idx[k];
  for (k=0; k<q->A/2; k++)
    Q->idx[k + q->A/2] = qint->idx[k] + q->ngauss/2;

  Q->ng = (int*) malloc(Q->A*sizeof(int));
  for (k=0; k<q->A/2; k++)
    Q->ng[k] = qint->ng[k];
  for (k=0; k<q->A/2; k++)
    Q->ng[k + q->A/2] = qint->ng[k];

  Q->G = (Gaussian*) malloc(Q->ngauss*sizeof(Gaussian));
  for (k=0; k<q->A/2; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->idx[k]+ki].xi = qint->xi[k]; 
  for (k=0; k<q->A/2; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->idx[k]+ki + q->ngauss/2].xi = qint->xi[k]; 
}


void FMDvxztoSlaterDet(const Para* q, SlaterDet* Q)
{
  FMDvxzInternals* qint = q->internals;

  int i;

  complex double chiu, chid, a, bx, by, bz;

  // the nucleons
  for (i=0; i<q->ngauss/2; i++) {
    chiu = q->x[i*NGAUSS+ 0] + I* q->x[i*NGAUSS+ 1];
    chid = q->x[i*NGAUSS+ 2] + I* q->x[i*NGAUSS+ 3];
    a =    q->x[i*NGAUSS+ 4] + I* q->x[i*NGAUSS+ 5];
    bx =   q->x[i*NGAUSS+ 6] + I* q->x[i*NGAUSS+ 7];
    by =   q->x[i*NGAUSS+ 8] + I* q->x[i*NGAUSS+ 9];
    bz =   q->x[i*NGAUSS+10] + I* q->x[i*NGAUSS+11];

    Q->G[i].chi[0] = chiu;
    Q->G[i].chi[1] = chid;
    Q->G[i].a      = a;
    Q->G[i].b[0]   = bx;
    Q->G[i].b[1]   = by;
    Q->G[i].b[2]   = bz;

    Q->G[i + q->ngauss/2].chi[0] = -chid;
    Q->G[i + q->ngauss/2].chi[1] = chiu;
    Q->G[i + q->ngauss/2].a      = a;
    Q->G[i + q->ngauss/2].b[0]   = bx;
    Q->G[i + q->ngauss/2].b[1]   = -by;
    Q->G[i + q->ngauss/2].b[2]   = bz;

  }
}


void FMDvxzprojectgradSlaterDet(const Para* q,
                                const gradSlaterDet* dQ, double* dq)
{
  FMDvxzInternals* qint = q->internals;
  int i;

  // gradients with respect to nucleons

  for (i=0; i<q->ngauss/2; i++) {
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

  for (i=0; i<q->ngauss/2; i++) {
    dq[i*NGAUSS+ 0] += 2.0*creal(dQ->gradval[i + q->ngauss/2].chi[1]);
    dq[i*NGAUSS+ 1] += 2.0*cimag(dQ->gradval[i + q->ngauss/2].chi[1]);
    dq[i*NGAUSS+ 2] += -2.0*creal(dQ->gradval[i + q->ngauss/2].chi[0]);
    dq[i*NGAUSS+ 3] += -2.0*cimag(dQ->gradval[i + q->ngauss/2].chi[0]);
    dq[i*NGAUSS+ 4] += 2.0*creal(dQ->gradval[i + q->ngauss/2].a);
    dq[i*NGAUSS+ 5] += 2.0*cimag(dQ->gradval[i + q->ngauss/2].a);
    dq[i*NGAUSS+ 6] += 2.0*creal(dQ->gradval[i + q->ngauss/2].b[0]);
    dq[i*NGAUSS+ 7] += 2.0*cimag(dQ->gradval[i + q->ngauss/2].b[0]);
    dq[i*NGAUSS+ 8] += -2.0*creal(dQ->gradval[i + q->ngauss/2].b[1]);
    dq[i*NGAUSS+ 9] += -2.0*cimag(dQ->gradval[i + q->ngauss/2].b[1]);
    dq[i*NGAUSS+10] += 2.0*creal(dQ->gradval[i + q->ngauss/2].b[2]);
    dq[i*NGAUSS+11] += 2.0*cimag(dQ->gradval[i + q->ngauss/2].b[2]);
  }

}

/**

  \file ParameterizationFMDd3h.c

  Parametrization of SlaterDet with d3h symmetry


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
#include "ParameterizationFMDd3h.h"

#include "misc/utils.h"
#include "misc/physics.h"


// rotation by 120 deg
#define R1xx -0.5
#define R1xy -0.5*sqrt(3.0)
#define R1yx 0.5*sqrt(3.0)
#define R1yy -0.5

#define R1uu 0.5
#define R1ud -0.5*sqrt(3.0)
#define R1du 0.5*sqrt(3.0)
#define R1dd 0.5

// rotation by 240 deg
#define R2xx -0.5
#define R2xy 0.5*sqrt(3.0)
#define R2yx -0.5*sqrt(3.0)
#define R2yy -0.5

#define R2uu -0.5
#define R2ud -0.5*sqrt(3.0)
#define R2du 0.5*sqrt(3.0)
#define R2dd -0.5


typedef struct {
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
} FMDd3hInternals;


Parameterization ParameterizationFMDd3h = {
  name : "FMDd3h",
  Pararead : FMDd3hread,
  Parawrite : FMDd3hwrite,
  Paraclone : FMDd3hclone,
  ParainitSlaterDet : FMDd3hinitSlaterDet,
  ParatoSlaterDet : FMDd3htoSlaterDet,
  ParaprojectgradSlaterDet : FMDd3hprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 12		


void FMDd3hclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = (double*) malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(FMDd3hInternals));
  FMDd3hInternals *qpint = qp->internals;
  FMDd3hInternals *qint = q->internals;

  qpint->ng = malloc(q->A/3*sizeof(int));
  for (i=0; i<q->A/3; i++)
    qpint->ng[i] = qint->ng[i];
  qpint->idx = malloc(q->A/3*sizeof(int));
  for (i=0; i<q->A/3; i++)
    qpint->idx[i] = qint->idx[i];
  qpint->xi = malloc(q->A/3*sizeof(int));
  for (i=0; i<q->A/3; i++)
    qpint->xi[i] = qint->xi[i];

}


#define BUFSIZE 255

int FMDd3hread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;
  FMDd3hInternals* qint;

  qint = q->internals = malloc(sizeof(FMDd3hInternals));
					
  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "FMDd3h")) {
    fprintf(stderr, "not a FMDd3h parameter set\n");
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

  qint->idx = (int*) malloc(q->A/3*sizeof(int));
  qint->ng = (int*) malloc(q->A/3*sizeof(int));
  qint->xi = (int*) malloc(q->A/3*sizeof(int));
  for (i=0; i<q->A/3; i++) {
    qint->idx[i] = 0;
    qint->ng[i] = 0;
    qint->xi[i] = 0;
  }
  
  q->n = q->ngauss/3*NGAUSS;

  q->x = (double*) malloc(q->n*sizeof(double));

  int k, xi;
  double chi0re, chi0im, chi1re, chi1im;
  double are, aim, b0re, b0im, b1re, b1im, b2re, b2im; 
  for (i=0; i<q->ngauss/3; i++) {
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
  for (k=1; k<q->A/3; k++) 
    qint->idx[k] = qint->idx[k-1]+qint->ng[k-1];

  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</Para>", 7)) {
    fprintf(stderr, "didn't find </Para>\n");
    return -1;
  }	

  return 0;
}


int FMDd3hwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  FMDd3hInternals* qint = q->internals;

  fprintf(fp, "<Para FMDd3h %s>\n", q->name);

  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);

  for (k=0; k<q->A/3; k++)
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


void FMDd3hinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  FMDd3hInternals* qint = q->internals;

  Q->A = q->A; 
  Q->Z = q->Z; 
  Q->N = q->N;
  Q->ngauss = q->ngauss;

  Q->idx = (int*) malloc(Q->A*sizeof(int));
  for (k=0; k<q->A/3; k++)
    Q->idx[k] = qint->idx[k];
  for (k=0; k<q->A/3; k++)
    Q->idx[k + q->A/3] = qint->idx[k] + q->ngauss/3;
  for (k=0; k<q->A/3; k++)
    Q->idx[k + 2*q->A/3] = qint->idx[k] + 2*q->ngauss/3;

  Q->ng = (int*) malloc(Q->A*sizeof(int));
  for (k=0; k<q->A/3; k++)
    Q->ng[k] = qint->ng[k];
  for (k=0; k<q->A/3; k++)
    Q->ng[k + q->A/3] = qint->ng[k];
  for (k=0; k<q->A/3; k++)
    Q->ng[k + 2*q->A/3] = qint->ng[k];

  Q->G = (Gaussian*) malloc(Q->ngauss*sizeof(Gaussian));
  for (k=0; k<q->A/3; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->idx[k]+ki].xi = qint->xi[k]; 
  for (k=0; k<q->A/3; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->idx[k]+ki + q->ngauss/3].xi = qint->xi[k]; 
  for (k=0; k<q->A/3; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->idx[k]+ki + 2*q->ngauss/3].xi = qint->xi[k]; 
}


void FMDd3htoSlaterDet(const Para* q, SlaterDet* Q)
{
  FMDd3hInternals* qint = q->internals;

  int i;

  complex double chiu, chid, a, bx, by, bz;

  // the nucleons
  for (i=0; i<q->ngauss/3; i++) {
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

    Q->G[i + q->ngauss/3].chi[0] = R1uu*chiu + R1ud*chid;
    Q->G[i + q->ngauss/3].chi[1] = R1du*chiu + R1dd*chid;
    Q->G[i + q->ngauss/3].a      = a;
    Q->G[i + q->ngauss/3].b[0]   = R1xx*bx + R1xy*by;
    Q->G[i + q->ngauss/3].b[1]   = R1yx*bx + R1yy*by;
    Q->G[i + q->ngauss/3].b[2]   = bz;

    Q->G[i + 2*q->ngauss/3].chi[0] = R2uu*chiu + R2ud*chid;
    Q->G[i + 2*q->ngauss/3].chi[1] = R2du*chiu + R2dd*chid;
    Q->G[i + 2*q->ngauss/3].a      = a;
    Q->G[i + 2*q->ngauss/3].b[0]   = R2xx*bx + R2xy*by;
    Q->G[i + 2*q->ngauss/3].b[1]   = R2yx*bx + R2yy*by;
    Q->G[i + 2*q->ngauss/3].b[2]   = bz;
  }
}


void FMDd3hprojectgradSlaterDet(const Para* q,
                                const gradSlaterDet* dQ, double* dq)
{
  FMDd3hInternals* qint = q->internals;
  int i;

  // gradients with respect to nucleons

  for (i=0; i<q->ngauss/3; i++) {
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

  for (i=0; i<q->ngauss/3; i++) {
    dq[i*NGAUSS+ 0] += R1uu* 2.0*creal(dQ->gradval[i + q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 1] += R1uu* 2.0*cimag(dQ->gradval[i + q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 0] += R1du* 2.0*creal(dQ->gradval[i + q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 1] += R1du* 2.0*cimag(dQ->gradval[i + q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 2] += R1ud* 2.0*creal(dQ->gradval[i + q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 3] += R1ud* 2.0*cimag(dQ->gradval[i + q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 2] += R1dd* 2.0*creal(dQ->gradval[i + q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 3] += R1dd* 2.0*cimag(dQ->gradval[i + q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 4] += 2.0*creal(dQ->gradval[i + q->ngauss/3].a);
    dq[i*NGAUSS+ 5] += 2.0*cimag(dQ->gradval[i + q->ngauss/3].a);
    dq[i*NGAUSS+ 6] += R1xx* 2.0*creal(dQ->gradval[i + q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 7] += R1xx* 2.0*cimag(dQ->gradval[i + q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 6] += R1yx* 2.0*creal(dQ->gradval[i + q->ngauss/3].b[1]);
    dq[i*NGAUSS+ 7] += R1yx* 2.0*cimag(dQ->gradval[i + q->ngauss/3].b[1]);
    dq[i*NGAUSS+ 8] += R1xy* 2.0*creal(dQ->gradval[i + q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 9] += R1xy* 2.0*cimag(dQ->gradval[i + q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 8] += R1yy* 2.0*creal(dQ->gradval[i + q->ngauss/3].b[1]);
    dq[i*NGAUSS+ 9] += R1yy* 2.0*cimag(dQ->gradval[i + q->ngauss/3].b[1]);
    dq[i*NGAUSS+10] += 2.0*creal(dQ->gradval[i + q->ngauss/3].b[2]);
    dq[i*NGAUSS+11] += 2.0*cimag(dQ->gradval[i + q->ngauss/3].b[2]);
  }

  for (i=0; i<q->ngauss/3; i++) {
    dq[i*NGAUSS+ 0] += R2uu* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 1] += R2uu* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 0] += R2du* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 1] += R2du* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 2] += R2ud* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 3] += R2ud* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].chi[0]);
    dq[i*NGAUSS+ 2] += R2dd* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 3] += R2dd* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].chi[1]);
    dq[i*NGAUSS+ 4] += 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].a);
    dq[i*NGAUSS+ 5] += 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].a);
    dq[i*NGAUSS+ 6] += R2xx* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 7] += R2xx* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 6] += R2yx* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].b[1]);
    dq[i*NGAUSS+ 7] += R2yx* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].b[1]);
    dq[i*NGAUSS+ 8] += R2xy* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 9] += R2xy* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].b[0]);
    dq[i*NGAUSS+ 8] += R2yy* 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].b[1]);
    dq[i*NGAUSS+ 9] += R2yy* 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].b[1]);
    dq[i*NGAUSS+10] += 2.0*creal(dQ->gradval[i + 2*q->ngauss/3].b[2]);
    dq[i*NGAUSS+11] += 2.0*cimag(dQ->gradval[i + 2*q->ngauss/3].b[2]);
  }
}

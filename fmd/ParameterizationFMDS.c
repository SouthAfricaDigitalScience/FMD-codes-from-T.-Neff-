/**

  \file ParameterizationFMDS.h

  Parametrization of SlaterDet.

  FMD Parameterization of Slater determinant with fixed spins


  (c) 2004 Thomas Neff

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
#include "ParameterizationFMDS.h"

#include "misc/utils.h"
#include "misc/physics.h"


typedef struct {
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
  int* chi;             ///< spin of nucleons
} FMDSInternals;


Parameterization ParameterizationFMDS = {
  name : "FMDS",
  Pararead : FMDSread,
  Parawrite : FMDSwrite,
  Paraclone : FMDSclone,
  ParainitSlaterDet : FMDSinitSlaterDet,
  ParatoSlaterDet : FMDStoSlaterDet,
  ParaprojectgradSlaterDet : FMDSprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 8		


void FMDSclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = (double*) malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(FMDSInternals));
  FMDSInternals *qpint = qp->internals;
  FMDSInternals *qint = q->internals;

  qpint->ng = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->ng[i] = qint->ng[i];
  qpint->idx = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->idx[i] = qint->idx[i];
  qpint->xi = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->xi[i] = qint->xi[i];
  qpint->chi = (int*) malloc(q->ngauss*sizeof(int));
  for (i=0; i<q->ngauss; i++)
    qpint->chi[i] = qint->chi[i];
}


#define BUFSIZE 255

int FMDSread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;
  FMDSInternals* qint;

  qint = q->internals = (FMDSInternals*) malloc(sizeof(FMDSInternals));
					

  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "FMDS")) {
    fprintf(stderr, "not a FMDS parameter set\n");
    return -1;
  }
  strcpy(q->name, stripstr(qname, ">"));

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

  qint->idx = (int*) malloc(q->A*sizeof(int));
  qint->ng = (int*) malloc(q->A*sizeof(int));
  qint->xi = (int*) malloc(q->A*sizeof(int));
  qint->chi = (int*) malloc(q->ngauss*sizeof(int));
  for (i=0; i<q->A; i++) {
    qint->idx[i] = 0;
    qint->ng[i] = 0;
    qint->xi[i] = 0;
  }
  for (i=0; i<q->ngauss; i++)
    qint->chi[i] = 0;

  q->n = q->ngauss*NGAUSS;

  q->x = (double*) malloc(q->n*sizeof(double));

  int k, xi, chi;
  double are, aim, b0re, b0im, b1re, b1im, b2re, b2im; 
  for (i=0; i<q->ngauss; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "%d %d %d "
	       "(%lf,%lf) (%lf,%lf) (%lf,%lf) (%lf,%lf)", 
	       &k, &xi, &chi,
	       &are, &aim, &b0re, &b0im, &b1re, &b1im, &b2re, &b2im) != 11) {
      fprintf(stderr, "malformed line\n>>%s<<\n", 
	      stripstr(buf,"\n"));
      return -1;
    }
    qint->ng[k]++;
    qint->xi[k] = xi;
    qint->chi[i] = chi;
    q->x[i*NGAUSS+ 0] = are;    q->x[i*NGAUSS+ 1] = aim;
    q->x[i*NGAUSS+ 2] = b0re;   q->x[i*NGAUSS+ 3] = b0im;
    q->x[i*NGAUSS+ 4] = b1re;   q->x[i*NGAUSS+ 5] = b1im;
    q->x[i*NGAUSS+ 6] = b2re;   q->x[i*NGAUSS+ 7] = b2im;
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


int FMDSwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  FMDSInternals* qint = q->internals;

  fprintf(fp, "<Para FMDS %s>\n", q->name);
  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);
  for (k=0; k<q->A; k++)
    for (ki=0; ki<qint->ng[k]; ki++) {
      i = qint->idx[k]+ki;
      fprintf(fp,
	      "%2d  %+2d  %+2d  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	      k, qint->xi[k], qint->chi[i],
	      q->x[i*NGAUSS+ 0], q->x[i*NGAUSS+ 1],
	      q->x[i*NGAUSS+ 2], q->x[i*NGAUSS+ 3],
	      q->x[i*NGAUSS+ 4], q->x[i*NGAUSS+ 5],
	      q->x[i*NGAUSS+ 6], q->x[i*NGAUSS+ 7]);
  }
  fprintf(fp, "</Para>\n");

  return 0;
}


void FMDSinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  FMDSInternals* qint = q->internals;

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


void FMDStoSlaterDet(const Para* q, SlaterDet* Q)
{
  FMDSInternals* qint = q->internals;
  int i;

  // this is a one-to-one mapping

  for (i=0; i<q->ngauss; i++) {
    Q->G[i].chi[0] = qint->chi[i] == +1 ? 1.0 : 0.0;
    Q->G[i].chi[1] = qint->chi[i] == -1 ? 1.0 : 0.0;
    Q->G[i].a      = q->x[i*NGAUSS+ 0] + I* q->x[i*NGAUSS+ 1];
    Q->G[i].b[0]   = q->x[i*NGAUSS+ 2] + I* q->x[i*NGAUSS+ 3];
    Q->G[i].b[1]   = q->x[i*NGAUSS+ 4] + I* q->x[i*NGAUSS+ 5];
    Q->G[i].b[2]   = q->x[i*NGAUSS+ 6] + I* q->x[i*NGAUSS+ 7];
  }
}


void FMDSprojectgradSlaterDet(const Para* q,
			     const gradSlaterDet* dQ, double* dq)
{
  int i;

  // only converting a complex in a real gradient

  for (i=0; i<q->ngauss; i++) {
    dq[i*NGAUSS+ 0] = 2.0*creal(dQ->gradval[i].a);
    dq[i*NGAUSS+ 1] = 2.0*cimag(dQ->gradval[i].a);
    dq[i*NGAUSS+ 2] = 2.0*creal(dQ->gradval[i].b[0]);
    dq[i*NGAUSS+ 3] = 2.0*cimag(dQ->gradval[i].b[0]);
    dq[i*NGAUSS+ 4] = 2.0*creal(dQ->gradval[i].b[1]);
    dq[i*NGAUSS+ 5] = 2.0*cimag(dQ->gradval[i].b[1]);
    dq[i*NGAUSS+ 6] = 2.0*creal(dQ->gradval[i].b[2]);
    dq[i*NGAUSS+ 7] = 2.0*cimag(dQ->gradval[i].b[2]);
  }
}

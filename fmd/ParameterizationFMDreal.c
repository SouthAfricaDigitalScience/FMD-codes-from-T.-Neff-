/**

  \file ParameterizationFMDreal.c

  Parametrization of SlaterDet.

  Parameterization uses only real parts of a,b for SlaterDet


  (c) 2005 Thomas Neff

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
#include "ParameterizationFMDreal.h"

#include "misc/utils.h"
#include "misc/physics.h"


typedef struct {
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
} FMDrealInternals;


Parameterization ParameterizationFMDreal = {
  name : "FMDreal",
  Pararead : FMDrealread,
  Parawrite : FMDrealwrite,
  Paraclone : FMDrealclone,
  ParainitSlaterDet : FMDrealinitSlaterDet,
  ParatoSlaterDet : FMDrealtoSlaterDet,
  ParaprojectgradSlaterDet : FMDrealprojectgradSlaterDet
};


// number of real parameter of a Gaussian
#define NGAUSS 8		


void FMDrealclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = (double*) malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(FMDrealInternals));
  FMDrealInternals *qpint = qp->internals;
  FMDrealInternals *qint = q->internals;

  qpint->ng = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->ng[i] = qint->ng[i];
  qpint->idx = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->idx[i] = qint->idx[i];
  qpint->xi = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++)
    qpint->xi[i] = qint->xi[i];

}


#define BUFSIZE 255

int FMDrealread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;
  FMDrealInternals* qint;

  qint = q->internals = (FMDrealInternals*) malloc(sizeof(FMDrealInternals));
					

  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "FMDreal")) {
    fprintf(stderr, "not a FMDreal parameter set\n");
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

  qint->idx = (int*) malloc(q->A*sizeof(int));
  qint->ng = (int*) malloc(q->A*sizeof(int));
  qint->xi = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++) {
    qint->idx[i] = 0;
    qint->ng[i] = 0;
    qint->xi[i] = 0;
  }
  
  q->n = q->ngauss*NGAUSS;

  q->x = (double*) malloc(q->n*sizeof(double));

  int k, xi;
  double chi0re, chi0im, chi1re, chi1im;
  double are, b0re, b1re, b2re; 
  for (i=0; i<q->ngauss; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "%d %d "
	       "(%lf,%lf) (%lf,%lf)  %lf  %lf %lf %lf", 
	       &k, &xi, &chi0re, &chi0im, &chi1re, &chi1im, 
	       &are, &b0re, &b1re, &b2re) != 10) {
      fprintf(stderr, "malformed line\n>>%s<<\n", 
	      stripstr(buf,"\n"));
      return -1;
    }
    qint->ng[k]++;
    qint->xi[k] = xi;
    q->x[i*NGAUSS+ 0] = chi0re; q->x[i*NGAUSS+ 1] = chi0im;
    q->x[i*NGAUSS+ 2] = chi1re; q->x[i*NGAUSS+ 3] = chi1im;
    q->x[i*NGAUSS+ 4] = are;
    q->x[i*NGAUSS+ 5] = b0re;
    q->x[i*NGAUSS+ 6] = b1re;
    q->x[i*NGAUSS+ 7] = b2re;
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


int FMDrealwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  FMDrealInternals* qint = q->internals;

  fprintf(fp, "<Para FMDreal %s>\n", q->name);
  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);
  for (k=0; k<q->A; k++)
    for (ki=0; ki<qint->ng[k]; ki++) {
      i = qint->idx[k]+ki;
      fprintf(fp,
	      "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  %12.8f  %12.8f %12.8f %12.8f\n",
	      k, qint->xi[k],
	      q->x[i*NGAUSS+ 0], q->x[i*NGAUSS+ 1],
	      q->x[i*NGAUSS+ 2], q->x[i*NGAUSS+ 3],
	      q->x[i*NGAUSS+ 4],
	      q->x[i*NGAUSS+ 5],
	      q->x[i*NGAUSS+ 6],
	      q->x[i*NGAUSS+ 7]);
  }
  fprintf(fp, "</Para>\n");

  return 0;
}


void FMDrealinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  FMDrealInternals* qint = q->internals;

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


void FMDrealtoSlaterDet(const Para* q, SlaterDet* Q)
{
  int i;

  // this is a one-to-one mapping

  for (i=0; i<q->ngauss; i++) {
    Q->G[i].chi[0] = q->x[i*NGAUSS+ 0] + I* q->x[i*NGAUSS+ 1];
    Q->G[i].chi[1] = q->x[i*NGAUSS+ 2] + I* q->x[i*NGAUSS+ 3];
    Q->G[i].a      = q->x[i*NGAUSS+ 4];
    Q->G[i].b[0]   = q->x[i*NGAUSS+ 5];
    Q->G[i].b[1]   = q->x[i*NGAUSS+ 6];
    Q->G[i].b[2]   = q->x[i*NGAUSS+ 7];
  }
}


void FMDrealprojectgradSlaterDet(const Para* q,
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
    dq[i*NGAUSS+ 5] = 2.0*creal(dQ->gradval[i].b[0]);
    dq[i*NGAUSS+ 6] = 2.0*creal(dQ->gradval[i].b[1]);
    dq[i*NGAUSS+ 7] = 2.0*creal(dQ->gradval[i].b[2]);
  }
}


void SlaterDetinitFMDreal(const SlaterDet* Q, Para* q)
{
  int i;
  FMDrealInternals* qint;

  strcpy(q->name, nucleusname(Q->A, Q->Z));

  q->A = Q->A;
  q->Z = Q->Z;
  q->N = Q->N;

  qint = q->internals = malloc(sizeof(FMDrealInternals));
					
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
    q->x[i*NGAUSS+ 5] = creal(Q->G[i].b[0]);
    q->x[i*NGAUSS+ 6] = creal(Q->G[i].b[1]);
    q->x[i*NGAUSS+ 7] = creal(Q->G[i].b[2]);
  }
}

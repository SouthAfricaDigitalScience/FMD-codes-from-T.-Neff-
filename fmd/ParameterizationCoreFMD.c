/**

  \file ParameterizationCoreFMD.h

  Parametrization of SlaterDet.

  Core plus FMD nucleons


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
#include "ParameterizationCoreFMD.h"

#include "misc/utils.h"
#include "misc/physics.h"


typedef struct {
  char corefname[255];
  char corefmd5hash[255];
  SlaterDet* core;
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
} CoreFMDInternals;


Parameterization ParameterizationCoreFMD = {
  name : "CoreFMD",
  Pararead : CoreFMDread,
  Parawrite : CoreFMDwrite,
  Paraclone : CoreFMDclone,
  ParainitSlaterDet : CoreFMDinitSlaterDet,
  ParatoSlaterDet : CoreFMDtoSlaterDet,
  ParaprojectgradSlaterDet : CoreFMDprojectgradSlaterDet
};


// number of parameters for Core
#define NCORE 3

// number of real parameter of a Gaussian
#define NGAUSS 12		


void CoreFMDclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = (double*) malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(CoreFMDInternals));
  CoreFMDInternals *qpint = qp->internals;
  CoreFMDInternals *qint = q->internals;

  strcpy(qpint->corefname, qint->corefname);
  strcpy(qpint->corefmd5hash, qint->corefmd5hash);

  qpint->core = malloc(sizeof(SlaterDet));
  cloneSlaterDet(qint->core, qpint->core);

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

int CoreFMDread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i;
  CoreFMDInternals* qint;

  qint = q->internals = malloc(sizeof(CoreFMDInternals));
  qint->core = malloc(sizeof(SlaterDet));
					
  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "CoreFMD")) {
    fprintf(stderr, "not a CoreFMD parameter set\n");
    return -1;
  }
  strcpy(q->name, stripstr(qname, ">"));

  // read core
  double corepos[3];

  fgets(buf, BUFSIZE, fp);
  sscanf(buf, "<Core (%lf, %lf, %lf) %s %s>",
	 &corepos[0], &corepos[1], &corepos[2],
	 qint->corefname, qint->corefmd5hash);
  if (strncmp(qint->corefmd5hash, md5hash(qint->corefname), 32)) {
    fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
    return -1;
  }
  readSlaterDetfromFile(qint->core, qint->corefname);

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

  qint->idx = (int*) malloc(q->A*sizeof(int));
  qint->ng = (int*) malloc(q->A*sizeof(int));
  qint->xi = (int*) malloc(q->A*sizeof(int));
  for (i=0; i<q->A; i++) {
    qint->idx[i] = 0;
    qint->ng[i] = 0;
    qint->xi[i] = 0;
  }
  
  q->n = q->ngauss*NGAUSS+NCORE;

  q->x = (double*) malloc(q->n*sizeof(double));

  q->x[0] = corepos[0]; q->x[1] = corepos[1]; q->x[2] = corepos[2];

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
    q->x[i*NGAUSS+ 0+NCORE] = chi0re; q->x[i*NGAUSS+ 1+NCORE] = chi0im;
    q->x[i*NGAUSS+ 2+NCORE] = chi1re; q->x[i*NGAUSS+ 3+NCORE] = chi1im;
    q->x[i*NGAUSS+ 4+NCORE] = are;    q->x[i*NGAUSS+ 5+NCORE] = aim;
    q->x[i*NGAUSS+ 6+NCORE] = b0re;   q->x[i*NGAUSS+ 7+NCORE] = b0im;
    q->x[i*NGAUSS+ 8+NCORE] = b1re;   q->x[i*NGAUSS+ 9+NCORE] = b1im;
    q->x[i*NGAUSS+10+NCORE] = b2re;   q->x[i*NGAUSS+11+NCORE] = b2im;
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


int CoreFMDwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  CoreFMDInternals* qint = q->internals;

  fprintf(fp, "<Para CoreFMD %s>\n", q->name);

  fprintf(fp, "<Core (%12.8f, %12.8f, %12.8f) %s %s\n", 
	  q->x[0], q->x[1], q->x[2],
	  qint->corefname, qint->corefmd5hash);

  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);

  for (k=0; k<q->A; k++)
    for (ki=0; ki<qint->ng[k]; ki++) {
      i = qint->idx[k]+ki;
      fprintf(fp,
	      "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	      k, qint->xi[k],
	      q->x[i*NGAUSS+ 0+NCORE], q->x[i*NGAUSS+ 1+NCORE],
	      q->x[i*NGAUSS+ 2+NCORE], q->x[i*NGAUSS+ 3+NCORE],
	      q->x[i*NGAUSS+ 4+NCORE], q->x[i*NGAUSS+ 5+NCORE],
	      q->x[i*NGAUSS+ 6+NCORE], q->x[i*NGAUSS+ 7+NCORE],
	      q->x[i*NGAUSS+ 8+NCORE], q->x[i*NGAUSS+ 9+NCORE],
	      q->x[i*NGAUSS+10+NCORE], q->x[i*NGAUSS+11+NCORE]);
  }
  fprintf(fp, "</Para>\n");

  return 0;
}


void CoreFMDinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki;
  CoreFMDInternals* qint = q->internals;

  Q->A = qint->core->A + q->A; 
  Q->Z = qint->core->Z + q->Z; 
  Q->N = qint->core->N + q->N;
  Q->ngauss = qint->core->ngauss + q->ngauss;

  Q->idx = (int*) malloc(Q->A*sizeof(int));
  for (k=0; k<qint->core->A; k++)
    Q->idx[k] = qint->core->idx[k];
  for (k=0; k<q->A; k++)
    Q->idx[qint->core->A + k] = qint->idx[k] + qint->core->ngauss;

  Q->ng = (int*) malloc(Q->A*sizeof(int));
  for (k=0; k<qint->core->A; k++)
    Q->ng[k] = qint->core->ng[k];
  for (k=0; k<q->A; k++)
    Q->ng[qint->core->A + k] = qint->ng[k];

  Q->G = (Gaussian*) malloc(Q->ngauss*sizeof(Gaussian));
  for (k=0; k<qint->core->A; k++) 
    for (ki=0; ki<qint->core->ng[k]; ki++)
      Q->G[qint->core->idx[k]+ki].xi = qint->core->G[k].xi; 
  for (k=0; k<q->A; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[qint->core->ngauss + qint->idx[k]+ki].xi = qint->xi[k]; 
}


void CoreFMDtoSlaterDet(const Para* q, SlaterDet* Q)
{
  CoreFMDInternals* qint = q->internals;
  int ngc = qint->core->ngauss;

  int i;

  // the core
  for (i=0; i<ngc; i++) {
    Q->G[i] = qint->core->G[i];
    moveGaussian(&Q->G[i], q->x);
  }

  // the nucleons
  for (i=0; i<q->ngauss; i++) {
    Q->G[ngc+i].chi[0] = q->x[i*NGAUSS+ 0+NCORE] + I* q->x[i*NGAUSS+ 1+NCORE];
    Q->G[ngc+i].chi[1] = q->x[i*NGAUSS+ 2+NCORE] + I* q->x[i*NGAUSS+ 3+NCORE];
    Q->G[ngc+i].a      = q->x[i*NGAUSS+ 4+NCORE] + I* q->x[i*NGAUSS+ 5+NCORE];
    Q->G[ngc+i].b[0]   = q->x[i*NGAUSS+ 6+NCORE] + I* q->x[i*NGAUSS+ 7+NCORE];
    Q->G[ngc+i].b[1]   = q->x[i*NGAUSS+ 8+NCORE] + I* q->x[i*NGAUSS+ 9+NCORE];
    Q->G[ngc+i].b[2]   = q->x[i*NGAUSS+10+NCORE] + I* q->x[i*NGAUSS+11+NCORE];
  }
}


void CoreFMDprojectgradSlaterDet(const Para* q,
				      const gradSlaterDet* dQ, double* dq)
{
  CoreFMDInternals* qint = q->internals;
  int ngc = qint->core->ngauss;
  int i,k;

  // gradients with respect to core
  for (i=0; i<3; i++)
    dq[i] = 0.0;

  for (k=0; k<qint->core->ngauss; k++)
    for (i=0; i<3; i++)
      dq[i] += 2.0*creal(dQ->gradval[k].b[i]);

  // gradients with respect to nucleons

  for (i=0; i<q->ngauss; i++) {
    dq[i*NGAUSS+ 0+NCORE] = 2.0*creal(dQ->gradval[ngc+i].chi[0]);
    dq[i*NGAUSS+ 1+NCORE] = 2.0*cimag(dQ->gradval[ngc+i].chi[0]);
    dq[i*NGAUSS+ 2+NCORE] = 2.0*creal(dQ->gradval[ngc+i].chi[1]);
    dq[i*NGAUSS+ 3+NCORE] = 2.0*cimag(dQ->gradval[ngc+i].chi[1]);
    dq[i*NGAUSS+ 4+NCORE] = 2.0*creal(dQ->gradval[ngc+i].a);
    dq[i*NGAUSS+ 5+NCORE] = 2.0*cimag(dQ->gradval[ngc+i].a);
    dq[i*NGAUSS+ 6+NCORE] = 2.0*creal(dQ->gradval[ngc+i].b[0]);
    dq[i*NGAUSS+ 7+NCORE] = 2.0*cimag(dQ->gradval[ngc+i].b[0]);
    dq[i*NGAUSS+ 8+NCORE] = 2.0*creal(dQ->gradval[ngc+i].b[1]);
    dq[i*NGAUSS+ 9+NCORE] = 2.0*cimag(dQ->gradval[ngc+i].b[1]);
    dq[i*NGAUSS+10+NCORE] = 2.0*creal(dQ->gradval[ngc+i].b[2]);
    dq[i*NGAUSS+11+NCORE] = 2.0*cimag(dQ->gradval[ngc+i].b[2]);
  }
}

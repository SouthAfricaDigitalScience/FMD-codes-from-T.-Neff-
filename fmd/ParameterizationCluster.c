/**

  \file ParameterizationCluster.c

  Parametrization of SlaterDet.

  FMD clusters


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
#include "ParameterizationCluster.h"

#include "misc/utils.h"
#include "misc/physics.h"

#include "numerics/rotationmatrices.h"


typedef struct {
  int ncluster;
  char (*clusterfname)[255];
  char (*clusterfmd5hash)[255];
  SlaterDet* cluster;
} ClusterInternals;


Parameterization ParameterizationCluster = {
  name : "Cluster",
  Pararead : Clusterread,
  Parawrite : Clusterwrite,
  Paraclone : Clusterclone,
  ParainitSlaterDet : ClusterinitSlaterDet,
  ParatoSlaterDet : ClustertoSlaterDet,
  ParaprojectgradSlaterDet : ClusterprojectgradSlaterDet
};


// number of parameters for Cluster (position + orientation)
#define NCLUSTER 6


void Clusterclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(ClusterInternals));
  ClusterInternals *qpint = qp->internals;
  ClusterInternals *qint = q->internals;

  qpint->ncluster = qint->ncluster;
  qpint->clusterfname = malloc(qint->ncluster*255*sizeof(char));
  qpint->clusterfmd5hash = malloc(qint->ncluster*255*sizeof(char));
  for (i=0; i<qint->ncluster; i++) {
    strcpy(qpint->clusterfname[i], qint->clusterfname[i]);
    strcpy(qpint->clusterfmd5hash[i], qint->clusterfmd5hash[i]);
  }

  qpint->cluster = malloc(qint->ncluster*sizeof(SlaterDet));
  for (i=0; i<qint->ncluster; i++)
    cloneSlaterDet(&qint->cluster[i], &qpint->cluster[i]);
}


#define BUFSIZE 255

int Clusterread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i,k;
  ClusterInternals* qint;

  qint = q->internals = malloc(sizeof(ClusterInternals));
					
  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "Cluster")) {
    fprintf(stderr, "not a Cluster parameter set\n");
    return -1;
  }
  strcpy(q->name, stripstr(qname, ">"));

  // read clusters
  int ncluster;

  fgets(buf, BUFSIZE, fp);
  if(sscanf(buf, "<Clusters %d>", &ncluster) != 1) {
    fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
    return -1;
  }
  qint->ncluster = ncluster;
  qint->clusterfname = malloc(ncluster*255*sizeof(char));
  qint->clusterfmd5hash = malloc(ncluster*255*sizeof(char));
  qint->cluster = malloc(ncluster*sizeof(SlaterDet));

  double clusterx[ncluster][6];
  for (i=0; i<ncluster; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "<Cluster (%lf, %lf, %lf) (%lf, %lf, %lf) %s %s>",
	       &clusterx[i][0], &clusterx[i][1], &clusterx[i][2],
	       &clusterx[i][3], &clusterx[i][4], &clusterx[i][5],
	       qint->clusterfname[i], qint->clusterfmd5hash[i]) != 8) {
      fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
      return -1;
    }
    if (strncmp(qint->clusterfmd5hash[i], md5hash(qint->clusterfname[i]), 32)) {
      fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
      return -1;
    }
    readSlaterDetfromFile(&qint->cluster[i], qint->clusterfname[i]);
  }
  
  q->n = ncluster*NCLUSTER;
  q->x = (double*) malloc(q->n*sizeof(double));

  for (k=0; k<ncluster; k++)
    for (i=0; i<NCLUSTER; i++)
      q->x[i+k*NCLUSTER] = clusterx[k][i];

  fgets(buf, BUFSIZE, fp);
  if (strncmp(buf, "</Para>", 7)) {
    fprintf(stderr, "didn't find </Para>\n");
    return -1;
  }	

  return 0;
}


int Clusterwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  ClusterInternals* qint = q->internals;

  fprintf(fp, "<Para Cluster %s>\n", q->name);

  fprintf(fp, "<Clusters %d>\n", qint->ncluster);

  for (k=0; k<qint->ncluster; k++) {
    fprintf(fp, "<Cluster (%12.8f, %12.8f, %12.8f) (%12.8f, %12.8f, %12.8f) %s %s\n", 
	    q->x[0+k*NCLUSTER], q->x[1+k*NCLUSTER], q->x[2+k*NCLUSTER],
	    q->x[3+k*NCLUSTER], q->x[4+k*NCLUSTER], q->x[5+k*NCLUSTER],
	    qint->clusterfname[k], qint->clusterfmd5hash[k]);
  }

  fprintf(fp, "</Para>\n");
  
  return 0;
}


void ClusterinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki, in;
  ClusterInternals* qint = q->internals;

  Q->A = 0; Q->Z = 0; Q->N = 0; Q->ngauss = 0;
  for (k=0; k<qint->ncluster; k++) {
    Q->A += qint->cluster[k].A; 
    Q->Z += qint->cluster[k].Z; 
    Q->N += qint->cluster[k].N;
    Q->ngauss += qint->cluster[k].ngauss;
  }

  // collect ngs
  Q->ng = (int*) malloc(Q->A*sizeof(int));
  in=0;
  for (k=0; k<qint->ncluster; k++) {
    for (ki=0; ki<qint->cluster[k].A; ki++) {
      Q->ng[in] = qint->cluster[k].ng[ki];
      in++;
    }
  }

  // build idx from ng
  Q->idx = (int*) malloc(Q->A*sizeof(int));
  Q->idx[0] = 0;
  for (k=1; k<Q->A; k++)
    Q->idx[k] = Q->idx[k-1]+Q->ng[k-1];

  // set xi of Gaussians
  Q->G = (Gaussian*) malloc(Q->ngauss*sizeof(Gaussian));
  in=0;
  for (k=0; k<qint->ncluster; k++) {
    for (ki=0; ki<qint->cluster[k].ngauss; ki++) {
      Q->G[in].xi = qint->cluster[k].G[ki].xi;
      in++;
    }
  }
}


void ClustertoSlaterDet(const Para* q, SlaterDet* Q)
{
  ClusterInternals* qint = q->internals;

  int i,k;
  double* euler;
  double* pos;
  complex double R2[2][2];
  double R3[3][3];

  // the clusters
  int gi = 0;
  for (k=0; k<qint->ncluster; k++) {
    pos = &q->x[0+k*NCLUSTER];
    euler = &q->x[3+k*NCLUSTER];
    rotatemat2(euler, R2);
    rotatemat3(euler, R3);
    for (i=0; i<qint->cluster[k].ngauss; i++) {
      Q->G[gi] = qint->cluster[k].G[i];
      rotateGaussian(&Q->G[gi], R3, R2);
      moveGaussian(&Q->G[gi], pos);
      gi++;
    }
  }
}


void ClusterprojectgradSlaterDet(const Para* q,
				 const gradSlaterDet* dQ, double* dq)
{
  ClusterInternals* qint = q->internals;

  int i,j,k,ki,gk;
  double* euler;
  complex double drotmat2[3][2][2];
  double drotmat3[3][3][3];
  

  complex double dchidang[2];
  complex double dbdang[3];
  // gradients with respect to clusters
  gk=0;
  for (k=0; k<qint->ncluster; k++) {
    for (i=0; i<NCLUSTER; i++)
      dq[i+k*NCLUSTER] = 0.0;

    // derivatives of the rotation matrices
    euler = &q->x[3+k*NCLUSTER];
    derivrotatemat2(euler, drotmat2);
    derivrotatemat3(euler, drotmat3);

    for (ki=0; ki<qint->cluster[k].ngauss; ki++) {

      // derivatives with respect to position
      for (i=0; i<3; i++)
	dq[i+k*NCLUSTER] += 2.0*creal(dQ->gradval[gk].b[i]);

      // derivatives with respect to angles
      for (i=0; i<3; i++) {
	for (j=0; j<2; j++)
	dchidang[j] = qint->cluster[k].G[ki].chi[j];
	cm2mult(drotmat2[i], dchidang);
	for (j=0; j<2; j++)
	  dq[3+i+k*NCLUSTER] += 2.0*creal(dQ->gradval[gk].chi[j]*conj(dchidang[j]));
	for (j=0; j<3; j++)
	  dbdang[j] = qint->cluster[k].G[ki].b[j];
	m3mult(drotmat3[i], dbdang);
	for (j=0; j<3; j++)
	  dq[3+i+k*NCLUSTER] += 2.0*creal(dQ->gradval[gk].b[j]*conj(dbdang[j]));
      }
      gk++;
    }
  }
}

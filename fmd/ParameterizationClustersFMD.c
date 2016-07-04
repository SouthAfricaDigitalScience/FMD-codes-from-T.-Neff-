/**

  \file ParameterizationClustersFMD.c

  Parametrization of SlaterDet.

  FMD clusters (position, velocity, orientation) plus FMD nucleons


  (c) 2006 Thomas Neff

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
#include "ParameterizationClustersFMD.h"

#include "misc/utils.h"
#include "misc/physics.h"

#include "numerics/rotationmatrices.h"


typedef struct {
  int ncluster;
  char (*clusterfname)[255];
  char (*clusterfmd5hash)[255];
  SlaterDet* cluster;
  int* ng;
  int* idx;
  int* xi;		///< isospin of nucleons
} ClustersFMDInternals;


Parameterization ParameterizationClustersFMD = {
  name : "ClustersFMD",
  Pararead : ClustersFMDread,
  Parawrite : ClustersFMDwrite,
  Paraclone : ClustersFMDclone,
  ParainitSlaterDet : ClustersFMDinitSlaterDet,
  ParatoSlaterDet : ClustersFMDtoSlaterDet,
  ParaprojectgradSlaterDet : ClustersFMDprojectgradSlaterDet
};


// number of parameters for Cluster (position, velocity, orientation)
#define NCLUSTER 9

// number of real parameter of a Gaussian
#define NGAUSS 12		


void ClustersFMDclone(const Para* q, Para* qp)
{
  int i;

  strcpy(qp->name, q->name);
  qp->A = q->A; qp->Z = q->Z; qp->N = q->N;
  qp->n = q->n;

  qp->x = malloc(q->n*sizeof(double));
  for (i=0; i<q->n; i++)
    qp->x[i] = q->x[i];

  qp->ngauss = q->ngauss;
  qp->internals = malloc(sizeof(ClustersFMDInternals));
  ClustersFMDInternals *qpint = qp->internals;
  ClustersFMDInternals *qint = q->internals;

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

int ClustersFMDread(FILE* fp, Para* q)
{
  char buf[BUFSIZE];
  int i,k;
  ClustersFMDInternals* qint;

  qint = q->internals = malloc(sizeof(ClustersFMDInternals));
					
  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Para ", 5) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Para ...>\n");
    return -1;
  }
  
  char Pname[80], qname[80];

  sscanf(buf, "<Para %s %s>", Pname, qname);
  if (strcmp(Pname, "ClustersFMD")) {
    fprintf(stderr, "not a ClustersFMD parameter set\n");
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

  double clusterx[ncluster][NCLUSTER];
  for (i=0; i<ncluster; i++) {
    fgets(buf, BUFSIZE, fp);
    if (sscanf(buf, "<Cluster (%lf, %lf, %lf) (%lf, %lf, %lf) (%lf, %lf, %lf) %s %s>",
	       &clusterx[i][0], &clusterx[i][1], &clusterx[i][2],
	       &clusterx[i][3], &clusterx[i][4], &clusterx[i][5],
	       &clusterx[i][6], &clusterx[i][7], &clusterx[i][8],
	       qint->clusterfname[i], qint->clusterfmd5hash[i]) != 11) {
      fprintf(stderr, "malformed line\n>>%s<<\n", stripstr(buf, "\n"));
      return -1;
    }
    if (strncmp(qint->clusterfmd5hash[i], md5hash(qint->clusterfname[i]), 32)) {
      fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
      return -1;
    }
    if (readSlaterDetfromFile(&qint->cluster[i], qint->clusterfname[i]))
      return -1;
  }

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
  
  q->n = q->ngauss*NGAUSS+ncluster*NCLUSTER;

  q->x = (double*) malloc(q->n*sizeof(double));

  for (k=0; k<ncluster; k++)
    for (i=0; i<NCLUSTER; i++)
      q->x[i+k*NCLUSTER] = clusterx[k][i];

  int nclust = ncluster*NCLUSTER;
  int xi;
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
    q->x[i*NGAUSS+ 0+nclust] = chi0re; q->x[i*NGAUSS+ 1+nclust] = chi0im;
    q->x[i*NGAUSS+ 2+nclust] = chi1re; q->x[i*NGAUSS+ 3+nclust] = chi1im;
    q->x[i*NGAUSS+ 4+nclust] = are;    q->x[i*NGAUSS+ 5+nclust] = aim;
    q->x[i*NGAUSS+ 6+nclust] = b0re;   q->x[i*NGAUSS+ 7+nclust] = b0im;
    q->x[i*NGAUSS+ 8+nclust] = b1re;   q->x[i*NGAUSS+ 9+nclust] = b1im;
    q->x[i*NGAUSS+10+nclust] = b2re;   q->x[i*NGAUSS+11+nclust] = b2im;
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


int ClustersFMDwrite(FILE* fp, const Para* q)
{
  int k, ki, i;
  ClustersFMDInternals* qint = q->internals;

  fprintf(fp, "<Para ClustersFMD %s>\n", q->name);

  fprintf(fp, "<Clusters %d>\n", qint->ncluster);

  for (k=0; k<qint->ncluster; k++) {
    fprintf(fp, "<Cluster (%12.8f, %12.8f, %12.8f) (%12.8f, %12.8f, %12.8f) (%12.8f, %12.8f, %12.8f) %s %s\n", 
	    q->x[0+k*NCLUSTER], q->x[1+k*NCLUSTER], q->x[2+k*NCLUSTER],
	    q->x[3+k*NCLUSTER], q->x[4+k*NCLUSTER], q->x[5+k*NCLUSTER],
	    q->x[6+k*NCLUSTER], q->x[7+k*NCLUSTER], q->x[8+k*NCLUSTER],
	    qint->clusterfname[k], qint->clusterfmd5hash[k]);
  }

  fprintf(fp, "%d %d %d\n", q->A, q->Z, q->N);
  fprintf(fp, "%d\n", q->ngauss);

  int nclust = qint->ncluster*NCLUSTER;

  for (k=0; k<q->A; k++)
    for (ki=0; ki<qint->ng[k]; ki++) {
      i = qint->idx[k]+ki;
      fprintf(fp,
	      "%2d  %+2d  (%12.8f,%12.8f)(%12.8f,%12.8f)  (%12.8f,%12.8f)  (%12.8f,%12.8f)(%12.8f,%12.8f)(%12.8f,%12.8f)\n",
	      k, qint->xi[k],
	      q->x[i*NGAUSS+ 0+nclust], q->x[i*NGAUSS+ 1+nclust],
	      q->x[i*NGAUSS+ 2+nclust], q->x[i*NGAUSS+ 3+nclust],
	      q->x[i*NGAUSS+ 4+nclust], q->x[i*NGAUSS+ 5+nclust],
	      q->x[i*NGAUSS+ 6+nclust], q->x[i*NGAUSS+ 7+nclust],
	      q->x[i*NGAUSS+ 8+nclust], q->x[i*NGAUSS+ 9+nclust],
	      q->x[i*NGAUSS+10+nclust], q->x[i*NGAUSS+11+nclust]);
  }
  fprintf(fp, "</Para>\n");

  return 0;
}


void ClustersFMDinitSlaterDet(const Para* q, SlaterDet* Q)
{
  int k, ki, in;
  ClustersFMDInternals* qint = q->internals;

  Q->A = q->A; Q->Z = q->Z; Q->N = q->N; Q->ngauss = q->ngauss;
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
  for (k=0; k<q->A; k++)
    Q->ng[in+k] = qint->ng[k];

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
  for (k=0; k<q->A; k++) 
    for (ki=0; ki<qint->ng[k]; ki++)
      Q->G[in + qint->idx[k]+ki].xi = qint->xi[k]; 

}


void ClustersFMDtoSlaterDet(const Para* q, SlaterDet* Q)
{
  ClustersFMDInternals* qint = q->internals;

  int i,k;
  double* euler;
  double* pos;
  double* vel;
  complex double R2[2][2];
  double R3[3][3];

  // the clusters
  int gi = 0;
  for (k=0; k<qint->ncluster; k++) {
    pos = &q->x[0+k*NCLUSTER];
    vel = &q->x[3+k*NCLUSTER];
    euler = &q->x[6+k*NCLUSTER];
    rotatemat2(euler, R2);
    rotatemat3(euler, R3);
    for (i=0; i<qint->cluster[k].ngauss; i++) {
      Q->G[gi] = qint->cluster[k].G[i];
      rotateGaussian(&Q->G[gi], R3, R2);
      boostGaussian(&Q->G[gi], vel);
      moveGaussian(&Q->G[gi], pos);
      gi++;
    }
  }

  // the nucleons
  int nclust = qint->ncluster*NCLUSTER;
  for (i=0; i<q->ngauss; i++) {
    Q->G[gi+i].chi[0] = q->x[i*NGAUSS+ 0+nclust] + I* q->x[i*NGAUSS+ 1+nclust];
    Q->G[gi+i].chi[1] = q->x[i*NGAUSS+ 2+nclust] + I* q->x[i*NGAUSS+ 3+nclust];
    Q->G[gi+i].a      = q->x[i*NGAUSS+ 4+nclust] + I* q->x[i*NGAUSS+ 5+nclust];
    Q->G[gi+i].b[0]   = q->x[i*NGAUSS+ 6+nclust] + I* q->x[i*NGAUSS+ 7+nclust];
    Q->G[gi+i].b[1]   = q->x[i*NGAUSS+ 8+nclust] + I* q->x[i*NGAUSS+ 9+nclust];
    Q->G[gi+i].b[2]   = q->x[i*NGAUSS+10+nclust] + I* q->x[i*NGAUSS+11+nclust];
  }
}


void ClustersFMDprojectgradSlaterDet(const Para* q,
				     const gradSlaterDet* dQ, double* dq)
{
  ClustersFMDInternals* qint = q->internals;
  Gaussian* G;

  int i,j,k,ki,gk;

  double* vel;
  double* euler;
  complex double dR2[3][2][2];
  double dR3[3][3][3];
  complex double R2[2][2];
  double R3[3][3];
  complex double Rb[3], Rchi[2];
  complex double phase;
  complex double dchidang[2];
  complex double dbdang[3];

  // gradients with respect to clusters
  gk=0;
  for (k=0; k<qint->ncluster; k++) {
    for (i=0; i<NCLUSTER; i++)
      dq[i+k*NCLUSTER] = 0.0;

    // velocity of cluster
    vel = &q->x[3+k*NCLUSTER];

    // rotation matrices and its derivatives
    euler = &q->x[6+k*NCLUSTER];
    rotatemat2(euler, R2);
    rotatemat3(euler, R3);
    derivrotatemat2(euler, dR2);
    derivrotatemat3(euler, dR3);

    for (ki=0; ki<qint->cluster[k].ngauss; ki++) {

      G = &qint->cluster[k].G[ki];

      for (j=0; j<2; j++)
	Rchi[j] = G->chi[j];
      cm2mult(R2, Rchi);

      for (j=0; j<3; j++)
	Rb[j] = G->b[j];
      m3mult(R3, Rb);

      // derivatives with respect to position
      for (i=0; i<3; i++)
	dq[i+k*NCLUSTER] += 2.0*creal(dQ->gradval[gk].b[i]);

      // derivatives with respect to velocity
      for (i=0; i<3; i++)
	dq[3+i+k*NCLUSTER] += 2.0*creal(conj(I*mass(G->xi)*G->a)*
					dQ->gradval[gk].b[i]);

      for (i=0; i<3; i++)
	for (j=0; j<2; j++)
	  dq[3+i+k*NCLUSTER] += 2.0*creal(conj(I*mass(G->xi)*Rb[i] - G->a*vel[i])*
					  dQ->gradval[gk].chi[j]);

      // derivatives with respect to angles

      phase = cexp(I*mass(G->xi)*(vel[0]*Rb[0]+vel[1]*Rb[1]+vel[2]*Rb[2])-
		   0.5*G->a*(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]));

      for (i=0; i<3; i++) {
	for (j=0; j<3; j++)
	  dbdang[j] = G->b[j];
	m3mult(dR3[i], dbdang);
	for (j=0; j<3; j++)
	  dq[6+i+k*NCLUSTER] += 2.0*creal(dQ->gradval[gk].b[j]*conj(dbdang[j]));

	for (j=0; j<2; j++)
	  dchidang[j] = phase* G->chi[j];
	cm2mult(dR2[i], dchidang);
	for (j=0; j<2; j++)
	  dchidang[j] += I*mass(G->xi)*
	    (vel[0]*dbdang[0]+vel[1]*dbdang[1]+vel[2]*dbdang[2])*phase*Rchi[j];
	
	for (j=0; j<2; j++)
	  dq[6+i+k*NCLUSTER] += 2.0*creal(phase*dQ->gradval[gk].chi[j]*conj(dchidang[j]));

      }

      gk++;
    }
  }

  // gradients with respect to nucleons
  int nclust = qint->ncluster*NCLUSTER;
  for (i=0; i<q->ngauss; i++) {
    dq[i*NGAUSS+ 0+nclust] = 2.0*creal(dQ->gradval[gk+i].chi[0]);
    dq[i*NGAUSS+ 1+nclust] = 2.0*cimag(dQ->gradval[gk+i].chi[0]);
    dq[i*NGAUSS+ 2+nclust] = 2.0*creal(dQ->gradval[gk+i].chi[1]);
    dq[i*NGAUSS+ 3+nclust] = 2.0*cimag(dQ->gradval[gk+i].chi[1]);
    dq[i*NGAUSS+ 4+nclust] = 2.0*creal(dQ->gradval[gk+i].a);
    dq[i*NGAUSS+ 5+nclust] = 2.0*cimag(dQ->gradval[gk+i].a);
    dq[i*NGAUSS+ 6+nclust] = 2.0*creal(dQ->gradval[gk+i].b[0]);
    dq[i*NGAUSS+ 7+nclust] = 2.0*cimag(dQ->gradval[gk+i].b[0]);
    dq[i*NGAUSS+ 8+nclust] = 2.0*creal(dQ->gradval[gk+i].b[1]);
    dq[i*NGAUSS+ 9+nclust] = 2.0*cimag(dQ->gradval[gk+i].b[1]);
    dq[i*NGAUSS+10+nclust] = 2.0*creal(dQ->gradval[gk+i].b[2]);
    dq[i*NGAUSS+11+nclust] = 2.0*cimag(dQ->gradval[gk+i].b[2]);
  }
}

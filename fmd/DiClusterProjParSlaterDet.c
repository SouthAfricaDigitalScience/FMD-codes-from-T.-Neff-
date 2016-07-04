/**

   \file DiClusterProjSlaterDet.c

   set of SlaterDet built out of two Clusters
   first Cluster has to be Jz eigenstate
   second Cluster is projected onto Jz eigenstate

   (c) 2005

*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "SlaterDet.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"
#include "DiClusterProjSlaterDet.h"

#include "misc/utils.h"
#include "misc/physics.h"
#include "numerics/zcw.h"
#include "numerics/wignerd.h"


typedef struct {
  SlaterDet* Qa;		///< SlaterDet of Cluster a
  SlaterDet* Qb;		///< intrinsic state of Cluster b
  double d;			///< distance of Clusters a and b
  int Ma;                       ///< Jz of Cluster a
  int *Jb;			///< angular momenta of Cluster b
  int *Mb;                      ///< Jz of Cluster b
  int *Kb;                      ///< K projection of Cluster b
  int *Pib;                     ///< Parity of Cluster b
  int izcw;                     ///< ZCW set usef for projection of cluster b
  int nangles;                  ///< number of angles
} DiClusterProjSlaterDetInternals;


MultiSlaterDet DiClusterProjSlaterDet = {
  A : 0,
  N : 0,
  n : 0,
  symmetry : DiClusterProjSlaterDetSymmetry,
  weight : DiClusterProjSlaterDetWeight,
  get : DiClusterProjSlaterDetGet,
  internals : NULL
};


// axial symmetry given by sum of K projections
Symmetry DiClusterProjSlaterDetSymmetry(const MultiSlaterDet* MSD, int iM)
{
  DiClusterProjSlaterDetInternals* internals = MSD->internals;

  Symmetry S=0;
  setSymmetry(&S, axial0 + internals->Ma + internals->Mb[iM]);
  
  return S;
}


complex double DiClusterProjSlaterDetWeight(const MultiSlaterDet* MSD, int iM, int i)
{
  DiClusterProjSlaterDetInternals* internals = MSD->internals;

  int izcw = internals->izcw;
  int nangles = internals->nangles;
  int J, M, K;
  J = internals->Jb[iM];
  M = internals->Mb[iM];
  K = internals->Kb[iM];

  double alpha,beta,gamma;
  getangles3(izcw, i, &alpha, &beta, &gamma);

  return 1.0/nangles* (J+1)* Djmkstar(J,M,K,alpha,beta,gamma);
}

// workspace
static SlaterDet *tmpQa, *tmpQb;

void joinintoSlaterDets(const SlaterDet* Qa, const SlaterDet* Qb, SlaterDet* Q)
{
  assert(Q->A == Qa->A+Qb->A);

  Q->Z = Qa->Z+Qb->Z;
  Q->N = Qa->N+Qb->N;
  Q->ngauss = Qa->ngauss+Qb->ngauss;

  int i;
  // copy Qa
  for (i=0; i<Qa->A; i++) {
    Q->idx[i] = Qa->idx[i];
    Q->ng[i] = Qa->ng[i];
  }
  for (i=0; i<Qa->ngauss; i++)
    Q->G[i] = Qa->G[i];

  // copy Qb
  for (i=0; i<Qb->A; i++) {
    Q->idx[i+Qa->A] = Qb->idx[i]+Qa->ngauss;
    Q->ng[i+Qa->A] = Qb->ng[i];
  }
  for (i=0; i<Qb->ngauss; i++) 
    Q->G[i+Qa->ngauss] = Qb->G[i];
}

// Parity not considered yet
void DiClusterProjSlaterDetGet(const MultiSlaterDet* MSD, int i, SlaterDet* Q)
{
  // initialize workspace
  if (!tmpQa) {
    tmpQa = malloc(sizeof(SlaterDet));
    tmpQb = malloc(sizeof(SlaterDet));
    allocateSlaterDet(tmpQa, MSD->A);
    allocateSlaterDet(tmpQb, MSD->A);
  }

  DiClusterProjSlaterDetInternals* internals = MSD->internals;

  int izcw = internals->izcw;
  double alpha,beta,gamma;
  getangles3(izcw, i, &alpha, &beta, &gamma);
  
  copySlaterDet(internals->Qa, tmpQa);
  copySlaterDet(internals->Qb, tmpQb);

  double massa = tmpQa->Z*mproton + tmpQa->N*mneutron;
  double massb = tmpQb->Z*mproton + tmpQb->N*mneutron;
  double mass = massa + massb;

  double xa[3] = {0,0,0};
  double xb[3] = {0,0,0};

  xa[2] = massb/mass* internals->d;
  xb[2] = -massa/mass* internals->d;

  moveSlaterDet(tmpQa, xa);

  rotateSlaterDet(tmpQb, alpha, beta, gamma);
  moveSlaterDet(tmpQb, xb);

  joinintoSlaterDets(tmpQa, tmpQb, Q);
}


#define BUFSIZE 255
int DiClusterProjSlaterDetRead(FILE* fp, MultiSlaterDet* MSD)
{
  *MSD = DiClusterProjSlaterDet;

  double d;
  int izcw;
  int N;
  SlaterDet Qa, Qb;
  char fnamea[255], fnameb[255];
  char md5hasha[255], md5hashb[255];
  int Ma;

  char buf[BUFSIZE];

  fgets(buf, BUFSIZE, fp);
  sscanf(buf, "<DiClusterProjSlaterDet d=%lf izcw=%d N=%d>", &d, &izcw, &N);
  
  int nangles=nangles3(izcw);

  int Jb[N];
  int Mb[N];
  int Kb[N];
  int Pib[N];

  int i;
  for (i=0; i<N; i++) {
    fgets(buf, BUFSIZE, fp);
    sscanf(buf, "<Projection %d %d %d %d>", &Jb[i], &Mb[i], &Kb[i], &Pib[i]);
  }

  fgets(buf, BUFSIZE, fp);
  sscanf(buf, "<Cluster %s %s %d>", fnamea, md5hasha, &Ma); 
  fgets(buf, BUFSIZE, fp);
  sscanf(buf, "<Cluster %s %s>", fnameb, md5hashb); 

  if (strncmp(md5hasha, md5hash(fnamea), 32)) {
    fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
    return -1;
  }
  if (strncmp(md5hashb, md5hash(fnameb), 32)) {
    fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
    return -1;
  }

  if (readSlaterDetfromFile(&Qa, fnamea)) {
    fprintf(stderr, "... couldn't open %s for reading\n", fnamea);
    return -1;
  }
  if (readSlaterDetfromFile(&Qb, fnameb)) {
    fprintf(stderr, "... couldn't open %s for reading\n", fnameb);
    return -1;
  }

  DiClusterProjSlaterDetInternals* internals = 
    malloc(sizeof(DiClusterProjSlaterDetInternals));

  internals->Qa = malloc(sizeof(SlaterDet));
  cloneSlaterDet(&Qa, internals->Qa);

  internals->Qb = malloc(sizeof(SlaterDet));
  cloneSlaterDet(&Qb, internals->Qb);

  internals->d = d;
  internals->Ma = Ma;

  internals->Jb = malloc(N*sizeof(int));
  internals->Mb = malloc(N*sizeof(int));
  internals->Kb = malloc(N*sizeof(int));
  internals->Pib = malloc(N*sizeof(int));

  for (i=0; i<N; i++) {
    internals->Jb[i] = Jb[i];
    internals->Mb[i] = Mb[i];
    internals->Kb[i] = Kb[i];
    internals->Pib[i] = Pib[i];
  }

  internals->izcw = izcw;
  internals->nangles = nangles;

  MSD->A = Qa.A + Qb.A;
  MSD->N = N;
  MSD->n = nangles;
  MSD->internals = internals;

  return 0;
}

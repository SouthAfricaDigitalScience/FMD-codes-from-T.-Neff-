/**

   \file DiClusterMultiProjSlaterDet.c

   set of MultiSlaterDet built out of two Clusters

   each cluster is a multiconfig state projected on a set of (J, Jz) values

   (c) 2006 Thomas Neff

*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "SlaterDet.h"
#include "Projection.h"
#include "Symmetry.h"

#include "MultiSlaterDet.h"
#include "DiClusterMultiProjSlaterDet.h"

#include "misc/utils.h"
#include "misc/physics.h"
#include "numerics/zcw.h"
#include "numerics/wignerd.h"


// angular momentum of clusters restricted by JMAX
// maximum number of SlaterDets per cluster CLMAX

#define CLMAX 5


typedef struct {
  SlaterDet* Qa;		///< SlaterDets of cluster a
  int na;			///< number of intrinsic states for cluster a
  int Na;			///< number of multiconfig states for cluster a
  int *Ja;
  int *Ma;
  int *Pia;
  int npia;			///< npia=1: no parity projection
  int izcwa;
  complex double (*wa)[CLMAX*(JMAX+1)];
  SlaterDet* Qb;		///< intrinsic states of Cluster b
  int nb;			///< number of intrinsic states for Cluster b
  int Nb;
  int *Jb;
  int *Mb;
  int *Pib;
  int npib;
  int izcwb;
  complex double (*wb)[CLMAX*(JMAX+1)];

  double d;			///< distance of Clusters a and b
} DiClusterMultiProjSlaterDetInternals;


MultiSlaterDet DiClusterMultiProjSlaterDet = {
  A : 0,
  N : 0,
  n : 0,
  symmetry : DiClusterMultiProjSlaterDetSymmetry,
  weight : DiClusterMultiProjSlaterDetWeight,
  get : DiClusterMultiProjSlaterDetGet,
  internals : NULL
};


// axial symmetry given by Ma+Mb 
Symmetry DiClusterMultiProjSlaterDetSymmetry(const MultiSlaterDet* MSD, int iM)
{
  DiClusterMultiProjSlaterDetInternals* internals = MSD->internals;
  int iMa = iM % internals->Na;
  int iMb = iM / internals->Na;

  Symmetry S=0;
  setSymmetry(&S, axial0 + internals->Ma[iMa] + internals->Mb[iMb]);

  return S;
}


complex double DiClusterMultiProjSlaterDetWeight(const MultiSlaterDet* MSD, int iM, int i)
{
  DiClusterMultiProjSlaterDetInternals* internals = MSD->internals;
  int iMa = iM % internals->Na;
  int iMb = iM / internals->Na;

  int npia = internals->npia;
  int izcwa = internals->izcwa;
  int nanglespa = nangles3(izcwa)*npia;

  int npib = internals->npib;
  int izcwb = internals->izcwb;
  int nanglespb = nangles3(izcwb)*npib;

  int Ja, Ma, Pia;
  Ja = internals->Ja[iMa];
  Ma = internals->Ma[iMa];
  Pia = internals->Pia[iMa]; 
  
  int Jb, Mb, Pib;
  Jb = internals->Jb[iMb];
  Mb = internals->Mb[iMb];
  Pib = internals->Pib[iMb]; 

  int ia = i % (nanglespa*internals->na);
  int ib = i / (nanglespa*internals->na);

  int iangpa = ia % nanglespa;
  int iangpb = ib % nanglespb;

  double alphaa,betaa,gammaa;
  getangles3(izcwa, iangpa/npia, &alphaa, &betaa, &gammaa);
  //  fprintf(stderr, "i: %d, iangpa: %d, alphaa: %5.2f, betaa: %5.2f, gammaa: %5.2f\n",
  //  i, iangpa, alphaa, betaa, gammaa);

  double alphab,betab,gammab;
  getangles3(izcwb, iangpb/npib, &alphab, &betab, &gammab);
  // fprintf(stderr, "i: %d, iangpb: %d, alphab: %5.2f, betab: %5.2f, gammab: %5.2f\n",
  //  i, iangpb, alphab, betab, gammab);

  complex double ca = 0.0;
  int Ka;
  for (Ka=-Ja; Ka<=Ja; Ka=Ka+2)
    ca += (npia==2 && Pia==-1 && ia%2 ? -1.0 : 1.0)/nanglespa* (Ja+1)* Djmkstar(Ja,Ma,Ka,alphaa,betaa,gammaa)*internals->wa[iMa][idxjm(Ja,Ka)+ia/nanglespa*(Ja+1)];

  complex double cb = 0.0;
  int Kb;
  for (Kb=-Jb; Kb<=Jb; Kb=Kb+2)
    cb += (npib==2 && Pib==-1 && ib%2 ? -1.0 : 1.0)/nanglespb* (Jb+1)* Djmkstar(Jb,Mb,Kb,alphab,betab,gammab)*internals->wb[iMb][idxjm(Jb,Kb)+ib/nanglespb*(Jb+1)];

  // fprintf(stderr, "Ja: %d, Ma: %d, ca: (%8.3f, %8.3f), Jb: %d, Mb: %d, cb: (%8.3f, %8.3f)\n",
  //  Ja, Ma, creal(ca), cimag(ca), Jb, Mb, creal(cb), cimag(cb));

  return ca*cb;

}

// workspace
static SlaterDet *tmpQa, *tmpQb;

static void joinintoSlaterDets(const SlaterDet* Qa, const SlaterDet* Qb, SlaterDet* Q)
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


void DiClusterMultiProjSlaterDetGet(const MultiSlaterDet* MSD, int i, SlaterDet* Q)
{
  // initialize workspace
  if (!tmpQa) {
    tmpQa = malloc(sizeof(SlaterDet));
    tmpQb = malloc(sizeof(SlaterDet));
    allocateSlaterDet(tmpQa, MSD->A);
    allocateSlaterDet(tmpQb, MSD->A);
  }

  DiClusterMultiProjSlaterDetInternals* internals = MSD->internals;

  int npia=internals->npia;
  int npib=internals->npib;

  int nanglespa = nangles3(internals->izcwa)*npia;
  int nanglespb = nangles3(internals->izcwb)*npib;

  int ia = i % (nanglespa*internals->na);
  int ib = i / (nanglespa*internals->na);

  // fprintf(stderr, "get: i=%d, na=%d, nb=%d, nanglespa=%d, nanglespb=%d, ia=%d, ib=%d\n",	i, internals->na, internals->nb, nanglespa, nanglespb, ia, ib);
  

  double alphaa,betaa,gammaa;
  getangles3(internals->izcwa, (ia % nanglespa)/npia, &alphaa, &betaa, &gammaa);

  double alphab,betab,gammab;
  getangles3(internals->izcwb, (ib % nanglespb)/npib, &alphab, &betab, &gammab);

  // fprintf(stderr, "tmpQa->A : %d\n", tmpQa->A);
  // fprintf(stderr, "Qa[%d].A : %d\n", ia/nanglespa, internals->Qa[ia/nanglespa].A);
  // fprintf(stderr, "tmpQb->A: %d\n", tmpQb->A);
  // fprintf(stderr, "Qb[%d].A : %d\n", ib/nanglespb, internals->Qb[ib/nanglespb].A);
  copySlaterDet(&internals->Qa[ia/nanglespa], tmpQa);
  copySlaterDet(&internals->Qb[ib/nanglespb], tmpQb);

  double massa = tmpQa->Z*mproton + tmpQa->N*mneutron;
  double massb = tmpQb->Z*mproton + tmpQb->N*mneutron;
  double mass = massa + massb;

  double xa[3] = {0,0,0};
  double xb[3] = {0,0,0};

  xa[2] = massb/mass* internals->d;
  xb[2] = -massa/mass* internals->d;

  rotateSlaterDet(tmpQa, alphaa, betaa, gammaa);
  if (npia==2 && ia%2)
    invertSlaterDet(tmpQa);
  moveSlaterDet(tmpQa, xa);

  rotateSlaterDet(tmpQb, alphab, betab, gammab);
  if (npib==2 && ib%2)
    invertSlaterDet(tmpQb);
  moveSlaterDet(tmpQb, xb);

  joinintoSlaterDets(tmpQa, tmpQb, Q);

}


static void sreadcvec(char* s, int n, complex double* a)
{
  char *c;
  double ref, imf;
  int l;

  c = strtok(s, " ,()");
  for (l=0; l<n; l++) {
    ref=atof(c); c=strtok(NULL, " ,()");
    imf=atof(c); c=strtok(NULL, " ,()");
    a[l] = ref + I*imf;
  }
}


#define BUFSIZE 4096
int DiClusterMultiProjSlaterDetRead(FILE* fp, MultiSlaterDet* MSD)
{
  *MSD = DiClusterMultiProjSlaterDet;

  DiClusterMultiProjSlaterDetInternals* internals = 
    malloc(sizeof(DiClusterMultiProjSlaterDetInternals));

  SlaterDet Qa, Qb;
  char fnamea[255], fnameb[255];
  char md5hasha[255], md5hashb[255];

  double d;
  int na, nb;
  int izcwa, izcwb;
  int npia, npib;
  int Na, Nb;

  char buf[BUFSIZE];
  fgets(buf, BUFSIZE, fp);
  sscanf(buf, "<DiClusterMultiProjSlaterDet d=%lf na=%d izcwa=%d npia=%d Na=%d nb=%d izcwb=%d npib=%d Nb=%d>", 
	 &d, &na, &izcwa, &npia, &Na, &nb, &izcwb, &npib, &Nb);

  internals->d = d;
  internals->na = na;
  internals->npia = npia;
  internals->izcwa = izcwa;
  internals->Na = Na;
  internals->nb = nb;
  internals->npib = npib;
  internals->izcwb = izcwb;
  internals->Nb = Nb;

  // Slater dets for Cluster a

  internals->Qa = malloc(na*sizeof(SlaterDet));

  int ia;
  for (ia=0; ia<na; ia++) {
    fgets(buf, BUFSIZE, fp);
    sscanf(buf, "<Cluster %s %s>", fnamea, md5hasha); 
    if (strncmp(md5hasha, md5hash(fnamea), 32)) {
      fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
      return -1;
    }
    if (readSlaterDetfromFile(&Qa, fnamea)) {
      fprintf(stderr, "... couldn't open %s for reading\n", fnamea);
      return -1;
    }
    cloneSlaterDet(&Qa, &internals->Qa[ia]);
  }

  internals->Ja = malloc(Na*sizeof(int));
  internals->Ma = malloc(Na*sizeof(int));
  internals->Pia = malloc(Na*sizeof(int));
  internals->wa = malloc(Na*CLMAX*(JMAX+1)*sizeof(complex double));

  int Ja;
  int Ma;
  int Pia;
  char ca[BUFSIZE];

  int Ia;
  for (Ia=0; Ia<Na; Ia++) {
    fgets(buf, BUFSIZE, fp);
    sscanf(buf, "<Projection %d %d %d %[() ,.0-9eE+-]>", &Ja, &Ma, &Pia, ca);

    // fprintf(stderr, "Ja: %d, Ma: %d, Pia: %d, coefficients: %s\n",
    //	    Ja, Ma, Pia, ca);

    internals->Ja[Ia] = Ja;
    internals->Ma[Ia] = Ma;
    internals->Pia[Ia] = Pia;
    sreadcvec(ca, na*(Ja+1), internals->wa[Ia]);
  }


  // Slater dets for Cluster b

  internals->Qb = malloc(nb*sizeof(SlaterDet));

  int ib;
  for (ib=0; ib<nb; ib++) {
    fgets(buf, BUFSIZE, fp);
    sscanf(buf, "<Cluster %s %s>", fnameb, md5hashb); 
    if (strncmp(md5hashb, md5hash(fnameb), 32)) {
      fprintf(stderr, "...  SlaterDetFile does not match with existing hash\n");
      return -1;
    }
    if (readSlaterDetfromFile(&Qb, fnameb)) {
      fprintf(stderr, "... couldn't open %s for reading\n", fnameb);
      return -1;
    }
    cloneSlaterDet(&Qb, &internals->Qb[ib]);
  }

  internals->Jb = malloc(Nb*sizeof(int));
  internals->Mb = malloc(Nb*sizeof(int));
  internals->Pib = malloc(Nb*sizeof(int));
  internals->wb = malloc(Nb*CLMAX*(JMAX+1)*sizeof(complex double));

  int Jb;
  int Mb;
  int Pib;
  char cb[BUFSIZE];

  int Ib;
  for (Ib=0; Ib<Nb; Ib++) {
    fgets(buf, BUFSIZE, fp);
    sscanf(buf, "<Projection %d %d %d %[() ,.0-9eE+-]>", &Jb, &Mb, &Pib, cb);

    internals->Jb[Ib] = Jb;
    internals->Mb[Ib] = Mb;
    internals->Pib[Ib] = Pib;
    sreadcvec(cb, nb*(Jb+1), internals->wb[Ib]);
  }

  MSD->A = internals->Qa[0].A + internals->Qb[0].A;
  MSD->N = Na*Nb;
  MSD->n = na*npia*nangles3(izcwa)*nb*npib*nangles3(izcwb);
  MSD->internals = internals;

  fprintf(stderr, "N: %d, n: %d\n", MSD->N, MSD->n);

  return 0;
}

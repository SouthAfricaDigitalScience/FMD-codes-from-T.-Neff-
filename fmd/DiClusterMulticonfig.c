/**

   \file DiClusterMulticonfig.c

   set of MultiSlaterDet built out of two Clusters

   each cluster is a multiconfig state projected on a set of (J, Jz, pi) values

   (c) 2009 Thomas Neff

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
#include "DiClusterMulticonfig.h"

#include "misc/utils.h"
#include "misc/physics.h"
#include "numerics/zcw.h"
#include "numerics/gaussquad.h"
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
  int nanga;
  double *alphaa, *betaa, *gammaa, *weighta;
  complex double (*wa)[CLMAX*(JMAX+1)];
  SlaterDet* Qb;		///< intrinsic states of Cluster b
  int nb;			///< number of intrinsic states for Cluster b
  int Nb;
  int *Jb;
  int *Mb;
  int *Pib;
  int npib;
  int nangb;
  double *alphab, *betab, *gammab, *weightb;
  complex double (*wb)[CLMAX*(JMAX+1)];

  double d;			///< distance of Clusters a and b
} DiClusterMulticonfigInternals;


MultiSlaterDet DiClusterMulticonfig = {
  A : 0,
  N : 0,
  n : 0,
  symmetry : DiClusterMulticonfigSymmetry,
  weight : DiClusterMulticonfigWeight,
  get : DiClusterMulticonfigGet,
  internals : NULL
};


// axial symmetry given by Ma+Mb 
Symmetry DiClusterMulticonfigSymmetry(const MultiSlaterDet* MSD, int iM)
{
  DiClusterMulticonfigInternals* internals = MSD->internals;
  int iMa = iM % internals->Na;
  int iMb = iM / internals->Na;

  Symmetry S=0;
  setSymmetry(&S, axial0 + internals->Ma[iMa] + internals->Mb[iMb]);

  return S;
}


complex double DiClusterMulticonfigWeight(const MultiSlaterDet* MSD, int iM, int i)
{
  DiClusterMulticonfigInternals* internals = MSD->internals;
  int iMa = iM % internals->Na;
  int iMb = iM / internals->Na;

  int npia = internals->npia;
  int nanglespa = internals->nanga*npia;

  int npib = internals->npib;
  int nanglespb = internals->nangb*npib;

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

  double alphaa, betaa, gammaa, weighta;
  alphaa = internals->alphaa[iangpa/npia];
  betaa = internals->betaa[iangpa/npia];
  gammaa = internals->gammaa[iangpa/npia];
  weighta = internals->weighta[iangpa/npia];
  // fprintf(stderr, "i: %d, iangpa: %d, alphaa: %5.2f, betaa: %5.2f, gammaa: %5.2f\n",
  //     i, iangpa, alphaa, betaa, gammaa);

  double alphab, betab, gammab, weightb;
  alphab = internals->alphab[iangpb/npib];
  betab = internals->betab[iangpb/npib];
  gammab = internals->gammab[iangpb/npib];
  weightb = internals->weightb[iangpb/npib];
  //  fprintf(stderr, "i: %d, iangpb: %d, alphab: %5.2f, betab: %5.2f, gammab: %5.2f\n",
  //      i, iangpb, alphab, betab, gammab);

  complex double ca = 0.0;
  int Ka;
  for (Ka=-Ja; Ka<=Ja; Ka=Ka+2)
    ca += 1.0/npia* (npia==2 && Pia==-1 && ia%2 ? -1.0 : 1.0)* weighta* (Ja+1)/(8*M_PI*M_PI)* Djmkstar(Ja,Ma,Ka,alphaa,betaa,gammaa)*internals->wa[iMa][idxjm(Ja,Ka)+ia/nanglespa*(Ja+1)];

  complex double cb = 0.0;
  int Kb;
  for (Kb=-Jb; Kb<=Jb; Kb=Kb+2)
    cb += 1.0/npib* (npib==2 && Pib==-1 && ib%2 ? -1.0 : 1.0)* weightb* (Jb+1)/(8*M_PI*M_PI)* Djmkstar(Jb,Mb,Kb,alphab,betab,gammab)*internals->wb[iMb][idxjm(Jb,Kb)+ib/nanglespb*(Jb+1)];

  //fprintf(stderr, "Ja: %d, Ma: %d, ca: (%8.3f, %8.3f), Jb: %d, Mb: %d, cb: (%8.3f, %8.3f)\n",
  //        Ja, Ma, creal(ca), cimag(ca), Jb, Mb, creal(cb), cimag(cb));

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


void DiClusterMulticonfigGet(const MultiSlaterDet* MSD, int i, SlaterDet* Q)
{
  // initialize workspace
  if (!tmpQa) {
    tmpQa = malloc(sizeof(SlaterDet));
    tmpQb = malloc(sizeof(SlaterDet));
    allocateSlaterDet(tmpQa, MSD->A);
    allocateSlaterDet(tmpQb, MSD->A);
  }

  DiClusterMulticonfigInternals* internals = MSD->internals;

  int npia = internals->npia;
  int nanglespa = internals->nanga*npia;

  int npib = internals->npib;
  int nanglespb = internals->nangb*npib;

  int ia = i % (nanglespa*internals->na);
  int ib = i / (nanglespa*internals->na);

  // fprintf(stderr, "get: i=%d, na=%d, nb=%d, nanglespa=%d, nanglespb=%d, ia=%d, ib=%d\n",	i, internals->na, internals->nb, nanglespa, nanglespb, ia, ib);
  
  double alphaa, betaa, gammaa;
  alphaa = internals->alphaa[(ia % nanglespa)/npia];
  betaa = internals->betaa[(ia % nanglespa)/npia];
  gammaa = internals->gammaa[(ia % nanglespa)/npia];

  double alphab, betab, gammab;
  alphab = internals->alphab[(ib % nanglespb)/npib];
  betab = internals->betab[(ib % nanglespb)/npib];
  gammab = internals->gammab[(ib % nanglespb)/npib];

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


static int initangles(const char* projpar, int* nang,
                      double** alpha, double** beta, double** gamma, double** weight)
{       
  char *c = projpar;

  c = strtok(c, "-");

  if (!strncmp(c, "none", 4)) {
    *nang = 1;
    *alpha  = malloc(sizeof(double));
    *beta   = malloc(sizeof(double));
    *gamma  = malloc(sizeof(double));
    *weight = malloc(sizeof(double));
    *alpha[0] = 0.0; *beta[0] = 0.0; *gamma[0] = 0.0; *weight[0] = 8*M_PI*M_PI;
    // weight is off by factor of 2J+1 !
  } else if (!strncmp(c, "zcw", 3)) {
    c=strtok(NULL, "-");
    int zcwidx = atoi(c);
    *nang = nangles3(zcwidx);
    *alpha  = malloc(*nang* sizeof(double));
    *beta   = malloc(*nang* sizeof(double));
    *gamma  = malloc(*nang* sizeof(double));
    *weight = malloc(*nang* sizeof(double));
    int i;
    for (i=0; i<*nang; i++) {
      getangles3(zcwidx, i, &(*alpha)[i], &(*beta)[i], &(*gamma)[i]);
      (*weight)[i] = (8*M_PI*M_PI)/(*nang);
    }
  } else {
    int reflectalpha=0;
    int reflectbeta=0;
    int reflectgamma=0;

    if (*c == 'r') {
      reflectalpha = 1;
      c++;
    }
    int nalpha = atoi(c);

    c=strtok(NULL, "-");
    if (*c == 'r') {
      reflectbeta = 1;
      c++;
    }
    int nbeta = atoi(c);

    c=strtok(NULL, "-");
    if (*c == 'r') {
      reflectgamma = 1;
      c++;
    }
    int ngamma = atoi(c);

    double aalpha[nalpha], walpha[nalpha];
    double acosb[2*nbeta], wcosb[2*nbeta];
    double agamma[ngamma], wgamma[ngamma];

    if (reflectalpha) 
      ShiftedPeriodicTrapezoidalPoints(nalpha, 0, M_PI, 0, aalpha, walpha);
    else
      ShiftedPeriodicTrapezoidalPoints(nalpha, 0, 2*M_PI, 0, aalpha, walpha);
 
    if (reflectbeta) 
      GaussLegendrePoints(2*nbeta, -1.0, 1.0, acosb, wcosb);
    else
      GaussLegendrePoints(nbeta, -1.0, 1.0, acosb, wcosb);

    if (reflectgamma) 
      ShiftedPeriodicTrapezoidalPoints(ngamma, 0, M_PI, 0, agamma, wgamma);
    else
      ShiftedPeriodicTrapezoidalPoints(ngamma, 0, 2*M_PI, 0, agamma, wgamma);

    *nang = nalpha*nbeta*ngamma;
    *alpha  = malloc(*nang* sizeof(double));
    *beta   = malloc(*nang* sizeof(double));
    *gamma  = malloc(*nang* sizeof(double));
    *weight = malloc(*nang* sizeof(double));

    int ia,ib,ig;
    for (ia=0; ia<nalpha; ia++)
      for (ib=0; ib<nbeta; ib++)
        for (ig=0; ig<ngamma; ig++) {
          (*alpha)[ia+ib*nalpha+ig*(nalpha*nbeta)] = aalpha[ia];
          (*gamma)[ia+ib*nalpha+ig*(nalpha*nbeta)] = agamma[ig];
          if (reflectbeta) {
            (*beta)[ia+ib*nalpha+ig*(nalpha*nbeta)] = acos(acosb[2*nbeta-ib-1]);
            (*weight)[ia+ib*nalpha+ig*(nalpha*nbeta)] = walpha[ia]*2*wcosb[2*nbeta-ib-1]*wgamma[ig];
          } else {
            (*beta)[ia+ib*nalpha+ig*(nalpha*nbeta)] = acos(acosb[nbeta-ib-1]);
            (*weight)[ia+ib*nalpha+ig*(nalpha*nbeta)] = walpha[ia]*wcosb[nbeta-ib-1]*wgamma[ig];
          }
        }
  }

  return 0;
}



#define BUFSIZE 4096
int DiClusterMulticonfigRead(FILE* fp, MultiSlaterDet* MSD)
{
  *MSD = DiClusterMulticonfig;

  DiClusterMulticonfigInternals* internals = 
    malloc(sizeof(DiClusterMulticonfigInternals));

  SlaterDet Qa, Qb;
  char fnamea[255], fnameb[255];
  char md5hasha[255], md5hashb[255];

  double d;
  int na, nb;
  char anga[15], angb[15];
  int npia, npib;
  int Na, Nb;

  char buf[BUFSIZE];
  fgets(buf, BUFSIZE, fp);
  sscanf(buf, "<DiClusterMulticonfig d=%lf na=%d anga=%s npia=%d Na=%d nb=%d angb=%s npib=%d Nb=%d>", 
	 &d, &na, anga, &npia, &Na, &nb, angb, &npib, &Nb);

  internals->d = d;
  internals->na = na;
  internals->npia = npia;
  internals->Na = Na;
  internals->nb = nb;
  internals->npib = npib;
  internals->Nb = Nb;

  // decode angular projection parameters

  initangles(anga, &internals->nanga,
             &internals->alphaa, &internals->betaa, &internals->gammaa, &internals->weighta);
  initangles(angb, &internals->nangb,
             &internals->alphab, &internals->betab, &internals->gammab, &internals->weightb);

  {
    int i;
    fprintf(stderr, "Angles Cluster A:\n");
    for (i=0; i<internals->nanga; i++)
      fprintf(stderr, "(alpha, beta, gamma): weight = (%8.3f, %8.3f, %8.3f): %8.5f\n",
              internals->alphaa[i], internals->betaa[i], internals->gammaa[i],
              internals->weighta[i]);
    fprintf(stderr, "Angles Cluster B:\n");
    for (i=0; i<internals->nangb; i++)
      fprintf(stderr, "(alpha, beta, gamma): weight = (%8.3f, %8.3f, %8.3f): %8.5f\n",
              internals->alphab[i], internals->betab[i], internals->gammab[i],
              internals->weightb[i]);
  }

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
  MSD->n = na*npia*internals->nanga*nb*npib*internals->nangb;
  MSD->internals = internals;

  fprintf(stderr, "N: %d, n: %d\n", MSD->N, MSD->n);

  return 0;
}

/**

   \file Formfactors.c

   calculate formfactors
   center of mass motion is treated exactly


   (c) 2004-2007 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "CenterofMass.h"

#include "Formfactors.h"

#include "numerics/cmath.h"
#include "numerics/gaussquad.h"
#include "numerics/sphericalharmonics.h"
#include "numerics/sphericalbessel.h"

#include "misc/utils.h"
#include "misc/physics.h"

#define SQR(x)	((x)*(x))


// The following Many-Body Operators have to initialized with
// int dim and FormfactorPara par before use

char* fflabel[NMULTIPOLES] = { "EMonopole", 
			       "EDipole", 
			       "EQuadrupole", 
			       "EOctupole" };

// Array of Multipole Formfactor Operators
ManyBodyOperator OpMultipoleFormfactor[NMULTIPOLES];


ManyBodyOperators OpMultipoleFormfactors = {
 n : NMULTIPOLES,
 Op : OpMultipoleFormfactor,
 me : calcMultipoleFormfactorsod
};


void initOpFormfactors(FormfactorPara* par)
{
  int i;

  OpMultipoleFormfactors.dim = 2*par->npoints;
  OpMultipoleFormfactors.size = 2*par->npoints;
  OpMultipoleFormfactors.par = par;

  for (i=0; i<NMULTIPOLES; i++) {
    char* name = malloc(60);
    sprintf(name, "%sFormfactors-%05.2f-%d-%d-%d%s",
	    fflabel[i],
	    par->qmax, par->npoints, par->nalpha, par->ncosb, 
	    par->recoil ? "-recoil" : "");
    OpMultipoleFormfactor[i].name = name;
    OpMultipoleFormfactor[i].rank = 2*i;
    OpMultipoleFormfactor[i].pi = i % 2;
    OpMultipoleFormfactor[i].dim = 2*par->npoints;
    OpMultipoleFormfactor[i].size = 2*par->npoints;
    OpMultipoleFormfactor[i].par = par;
    OpMultipoleFormfactor[i].me = NULL;
  }
}


// static SlaterDet and SlaterDetAux workspace
static SlaterDet* Qboost;
static SlaterDetAux* Xboost;

static void ob_formfactorq(double q[3],
			   const Gaussian* G1, const Gaussian* G2,
			   const GaussianAux* X, complex double ffactor[2])
{
  complex double f;

  f = X->R*X->S* cexp(I*(X->rho[0]*q[0]+X->rho[1]*q[1]+X->rho[2]*q[2]) - 
		      0.5*X->alpha*(SQR(q[0])+SQR(q[1])+SQR(q[2])));

  if (G1->xi == 1)
    ffactor[0] += f* X->T;
  if (G1->xi == -1)
    ffactor[1] += f* X->T;
}
  

// calculate off-diagonal proton and neutron point formfactor F(q)
void calcFormfactorod(double q[3], int recoil,
		      const SlaterDet* Q, const SlaterDet* Qp,
		      const SlaterDetAux* X,
		      complex double ffactor[2])
{
  // inititialize workspace if necessary
  if (recoil && !Qboost) {
    Qboost = malloc(sizeof(SlaterDet));
    initSlaterDet(Qp, Qboost);
    Xboost = malloc(sizeof(SlaterDetAux));
    initSlaterDetAux(Qp, Xboost);
  }

  OneBodyOperator op_ob_formfactorq = {dim: 2, opt: 1, par: q,
				       me : ob_formfactorq};

  if (recoil) { 
    double mass = Q->Z*mproton+Q->N*mneutron;
    double Vcm[3] = {-q[0]/mass, -q[1]/mass, -q[2]/mass};

    copySlaterDet(Qp, Qboost);
    boostSlaterDet(Qboost, Vcm);
    calcSlaterDetAuxod(Q, Qboost, Xboost);
    calcSlaterDetOBMEod(Q, Qboost, Xboost, &op_ob_formfactorq, ffactor);
  } else 
    calcSlaterDetOBMEod(Q, Qp, X, &op_ob_formfactorq, ffactor);
}


int il[NMULTIPOLES] = {0, 1, 4, 9};
int diml[NMULTIPOLES] = {1, 3, 5, 7};

void calcMultipoleFormfactorsod(FormfactorPara* par,
				const SlaterDet* Q, const SlaterDet* Qp,
				const SlaterDetAux* X,
				complex double* ffactor)
{
  double qmax = par->qmax;
  int ialpha, nalpha = par->nalpha;
  int icosb, ncosb = par->ncosb;
  int npoints = par->npoints;
  int recoil = par->recoil;
  int m, i, t, l;
  double alpha, beta, weight;
  double q, k[3];
  complex double ff[2];

  for (l=0; l<NMULTIPOLES; l++)
    for (t=0; t<2; t++)
      for (i=0; i<npoints; i++)
	for (m=0; m<diml[l]; m++)
	  ffactor[m + i*diml[l] + t*diml[l]*npoints + il[l]*2*npoints] = 0.0;

  double alphal[nalpha], walphal[nalpha];
  double cosbl[ncosb], wcosbl[ncosb];
  
  ShiftedPeriodicTrapezoidalPoints(nalpha, 0, 2*M_PI, 0.0, alphal, walphal);
  GaussLegendrePoints(ncosb, -1, 1, cosbl, wcosbl); 

  for (i=0; i<npoints; i++) {
    q = qmax* i/(npoints-1);
    for (ialpha=0; ialpha<nalpha; ialpha++)
      for (icosb=0; icosb<ncosb; icosb++) {
	alpha = alphal[ialpha];
	beta = acos(cosbl[icosb]);
	weight = walphal[ialpha]*wcosbl[icosb];

	k[0] = q*cos(alpha)*sin(beta);
	k[1] = q*sin(alpha)*sin(beta);
	k[2] = q*cos(beta);

	calcFormfactorod(k, recoil, Q, Qp, X, ff);
	for (l=0; l<NMULTIPOLES; l++)
	  for (t=0; t<2; t++)
	    for (m=-l; m<=l; m++)
	      ffactor[m+l + i*diml[l] + t*diml[l]*npoints + il[l]*2*npoints] += 
		weight* ff[t]* Y(2*l,2*m,beta,alpha);
	}
    }
}


// outpout routines

// Formfactors are saved as F(q) and not as |F(q)|^2 !

void writePointFormfactors(FILE* fp,
			   const Projection* P,
			   const FormfactorPara* p,
			   int l,
			   int j, int pi, int a,
			   void* ffactorme,
			   const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = ffctme[ipj][idx];
  int i;
  double q, ml, mlp, mln;
  
  for (i=0; i<npoints; i++) {
    q = qmax*i/(npoints-1);
    mlp = ff[0][i]/norm;
    mln = ff[1][i]/norm;
    ml = (ff[0][i]+ff[1][i])/norm;

    fprintf(fp, "%6.3f\t%13.8g\t%13.8g\t%13.8g\n",
	    q,
	    sqrt(j+1.0)*ml,
	    sqrt(j+1.0)*mlp,
	    sqrt(j+1.0)*mln);
  }
}


void writeChargeFormfactors(FILE* fp,
			    const Projection* P,
			    const FormfactorPara* p,
			    int l,
			    int j, int pi, int a,
			    void* ffactorme,
			    const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  double norm = E->norm[ipj][idx];
  complex double (**ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = ffctme[ipj][idx];
  int i;
  double q;
  double ml;
  
  for (i=0; i<npoints; i++) {
    q = qmax*i/(npoints-1);
    
    ml = creal((ff[0][i]*fproton(q) + ff[1][i]*fneutron(q))/norm);    

    fprintf(fp, "%6.3f\t%13.5g\n", q, sqrt(j+1.0)*ml);
  }
}


void writeTransitionPointFormfactors(FILE* fp,
                                     const Projection* P,
                                     const FormfactorPara* p,
                                     int l,
                                     int jfin, int pfin, int afin,
                                     int jini, int pini, int aini,
                                     void* ffactorme,
                                     const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipjfin = idxpij(jmax,pfin,jfin);
  int ipjini = idxpij(jmax,pini,jini);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idxfin = E->index[ipjfin][afin];
  int idxini = E->index[ipjini][aini];
  double normfin = E->norm[ipjfin][idxfin];
  double normini = E->norm[ipjini][idxini];
  complex double (****ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = 
    ffctme[ipjfin][ipjini][idxfin][idxini];
  int i;
  double q;
  double mlp, mln, ml;
  complex double ff0p, ff0n, phase;
  
  // get rid of possible phase factors
  // formfactor should be real and positive for q=0
  
  int ifix=npoints/3;

  // average over proton and neutron phasese (is this correct ?)
  ff0p = ff[0][ifix]; ff0n = ff[1][ifix];
  phase = 0.5*(ff0p/cabs(ff0p)+ff0n/cabs(ff0n));

  for (i=0; i<npoints; i++) {
    q = qmax*i/(npoints-1);

    mlp = creal(ff[0][i]/phase)/sqrt(normfin*normini);
    mln = creal(ff[1][i]/phase)/sqrt(normfin*normini);
    ml = creal((ff[0][i]+ff[1][i])/phase)/sqrt(normfin*normini);

    fprintf(fp, "%6.3f\t%13.8g\t%13.8g\t%13.8g\n",
	    q,
	    sqrt(jfin+1.0)*ml,
	    sqrt(jfin+1.0)*mlp,
	    sqrt(jfin+1.0)*mln);
  }
}


void writeTransitionChargeFormfactors(FILE* fp,
				      const Projection* P,
				      const FormfactorPara* p,
				      int l,
				      int jfin, int pfin, int afin,
				      int jini, int pini, int aini,
				      void* ffactorme,
				      const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipjfin = idxpij(jmax,pfin,jfin);
  int ipjini = idxpij(jmax,pini,jini);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idxfin = E->index[ipjfin][afin];
  int idxini = E->index[ipjini][aini];
  double normfin = E->norm[ipjfin][idxfin];
  double normini = E->norm[ipjini][idxini];
  complex double (****ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = 
    ffctme[ipjfin][ipjini][idxfin][idxini];
  int i;
  double q;
  double ml;
  complex double ff0, phase;
  
  // get rid of possible phase factors
  // formfactor should be real and positive for q=0
  
  int ifix=npoints/3;
  double qfix=qmax*ifix/(npoints-1);

  ff0 = ff[0][ifix]*fproton(qfix)+ff[1][ifix]*fneutron(qfix);
  phase = ff0/cabs(ff0);

  for (i=0; i<npoints; i++) {
    q = qmax*i/(npoints-1);
    
    ml = creal((ff[0][i]*fproton(q) + ff[1][i]*fneutron(q))/phase)/
      sqrt(normfin*normini);

    fprintf(fp, "%6.3f\t%13.5g\n", q, sqrt(jfin+1)*ml);
  }
}


void writePointDensities(FILE* fp,
			 const Projection* P,
			 const FormfactorPara* p,
			 double rmax, int npointsr,
			 int A, int l,
			 int j, int pi, int a,
			 void* ffactorme,
			 const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = ffctme[ipj][idx];
  complex double anu[npoints][2];

  int i, nu;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writePointDensities not yet implemented for l != 0\n");
    exit(-1);
  }

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu][0] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*ff[0][nu]/norm;
    anu[nu][1] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*ff[1][nu]/norm;
  }

  // calculate densities 
      
  double deltar=rmax/(npointsr-1);
  double rho, rhop, rhon;

  for (i=0; i<npointsr; i++) {

    rho = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rho += (anu[nu][0]+anu[nu][1])*besselj0(qnu[nu]*i*deltar);

    rhop = 0.0;
    for (nu=1; nu<npoints; nu++)
      rhop += anu[nu][0]* besselj0(qnu[nu]*i*deltar);

    rhon = 0.0;
    for (nu=1; nu<npoints; nu++)
      rhon += anu[nu][1]* besselj0(qnu[nu]*i*deltar);

    fprintf(fp, "%6.3f\t%15.8g\t%15.8g\t%15.8g\n", 
	    i*deltar, rho, rhop, rhon);
  }
}  
  

// use proton charge form factor to smear out protons and neutrons

void writeMatterDensities(FILE* fp,
                          const Projection* P,
                          const FormfactorPara* p,
                          double rmax, int npointsr,
                          int A, int l,
                          int j, int pi, int a,
                          void* ffactorme,
                          const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = ffctme[ipj][idx];
  complex double anu[npoints][2];

  int i, nu;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writeMatterDensities not yet implemented for l != 0\n");
    exit(-1);
  }

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu][0] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*ff[0][nu]*fproton(qnu[nu])/norm;
    anu[nu][1] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*ff[1][nu]*fproton(qnu[nu])/norm;
  }

  // calculate densities 
      
  double deltar=rmax/(npointsr-1);
  double rho, rhop, rhon;

  for (i=0; i<npointsr; i++) {

    rho = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rho += (anu[nu][0]+anu[nu][1])*besselj0(qnu[nu]*i*deltar);

    rhop = 0.0;
    for (nu=1; nu<npoints; nu++)
      rhop += anu[nu][0]* besselj0(qnu[nu]*i*deltar);

    rhon = 0.0;
    for (nu=1; nu<npoints; nu++)
      rhon += anu[nu][1]* besselj0(qnu[nu]*i*deltar);

    fprintf(fp, "%6.3f\t%15.8g\t%15.8g\t%15.8g\n", 
	    i*deltar, rho, rhop, rhon);
  }
}  


void writeChargeDensities(FILE* fp,
			  const Projection* P,
			  const FormfactorPara* p,
			  double rmax, int npointsr,
			  int Z, int l,
			  int j, int pi, int a,
			  void* ffactorme,
			  const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = ffctme[ipj][idx];
  double anu[npoints];

  int i, nu;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writeChargeDensities not yet implemented for l != 0\n");
    exit(-1);
  }

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*
      (ff[0][nu]*fproton(qnu[nu])+ff[1][nu]*fneutron(qnu[nu]))/norm;
  }

  // calculate densities 
      
  double deltar=rmax/(npointsr-1);
  double rho;

  // charge densities
  for (i=0; i<npointsr; i++) {
    rho = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rho += anu[nu]*besselj0(qnu[nu]*i*deltar);
    fprintf(fp, "%6.3f\t%15.8g\n", i*deltar, rho);
  }
}  


void writeChargeFourierBessel(FILE* fp,
			      const Projection* P,
			      const FormfactorPara* p,
			      double rmax, int npointsr,
			      int Z, int l,
			      int j, int pi, int a,
			      void* ffactorme,
			      const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = ffctme[ipj][idx];
  double anu[npoints];

  int nu;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writeChargeFourierBessel not yet implemented for l != 0\n");
    exit(-1);
  }

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*
      (ff[0][nu]*fproton(qnu[nu])+ff[1][nu]*fneutron(qnu[nu]))/norm;
  }

  // write fourier bessel coefficients

  fprintf(fp, "# Rc      = %8.3f\n", Rc);
  fprintf(fp, "# qmax    = %8.3f\n", qmax);
  fprintf(fp, "# npoints = %d\n", npoints-1);

  fprintf(fp, "\n#nu\tq[nu]\ta[nu]\n");

  for (nu=1; nu<npoints; nu++)
    fprintf(fp, "%d\t%8.3f\t%15.8g\n", nu, qnu[nu], anu[nu]);

}  


void writeTransitionPointDensities(FILE* fp,
				   const Projection* P,
				   const FormfactorPara* p,
				   double rmax, int npointsr,
				   int A, int l,
				   int jfin, int pfin, int afin,
				   int jini, int pini, int aini,
				   void* ffactorme,
				   const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipjini = idxpij(jmax,pini,jini);
  int ipjfin = idxpij(jmax,pfin,jfin);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idxfin = E->index[ipjfin][afin];
  int idxini = E->index[ipjini][aini];
  double normfin = E->norm[ipjfin][idxfin];
  double normini = E->norm[ipjini][idxini];
  complex double (****ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = 
    ffctme[ipjfin][ipjini][idxfin][idxini];
  complex double anu[npoints][2];

  int i, nu;
  complex double ff0, phase;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writeTransitionPointDensities not yet implemented for l != 0\n");
    exit(-1);
  }

  // get rid of possible phase factor

  int ifix=npoints/3;
  double qfix=qmax*ifix/(npoints-1);

  ff0 = ff[0][ifix]*fproton(qfix)+ff[1][ifix]*fneutron(qfix);
  phase = ff0/cabs(ff0);

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu][0] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*
      creal(ff[0][nu]/phase)/sqrt(normfin*normini);
    anu[nu][1] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*
      creal(ff[1][nu]/phase)/sqrt(normfin*normini);
  }

  // calculate densities 
      
  double deltar=rmax/(npointsr-1);
  double rho, rhop, rhon;

  for (i=0; i<npointsr; i++) {

    rho = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rho += (anu[nu][0]+anu[nu][1])*besselj0(qnu[nu]*i*deltar);

    rhop = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rhop += anu[nu][0]*besselj0(qnu[nu]*i*deltar);

    rhon = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rhon += anu[nu][1]*besselj0(qnu[nu]*i*deltar);

    fprintf(fp, "%6.3f\t%15.8g\t%15.8g\t%15.8g\n", 
	    i*deltar, rho, rhop, rhon);
  }
}  


void writeTransitionChargeDensities(FILE* fp,
				    const Projection* P,
				    const FormfactorPara* p,
				    double rmax, int npointsr,
				    int Z, int l,
				    int jfin, int pfin, int afin,
				    int jini, int pini, int aini,
				    void* ffactorme,
				    const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipjini = idxpij(jmax,pini,jini);
  int ipjfin = idxpij(jmax,pfin,jfin);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idxfin = E->index[ipjfin][afin];
  int idxini = E->index[ipjini][aini];
  double normfin = E->norm[ipjfin][idxfin];
  double normini = E->norm[ipjini][idxini];
  complex double (****ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = 
    ffctme[ipjfin][ipjini][idxfin][idxini];
  complex double anu[npoints];

  int i, nu;
  complex double ff0, phase;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writeChargeDensities not yet implemented for l != 0\n");
    exit(-1);
  }

  // get rid of possible phase factor

  int ifix=npoints/3;
  double qfix=qmax*ifix/(npoints-1);

  ff0 = ff[0][ifix]*fproton(qfix)+ff[1][ifix]*fneutron(qfix);
  phase = ff0/cabs(ff0);

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*
      creal((ff[0][nu]*fproton(qnu[nu])+ff[1][nu]*fneutron(qnu[nu]))/phase)/
      sqrt(normfin*normini);
  }

  // calculate densities 
      
  double deltar=rmax/(npointsr-1);
  double rho;

  // charge densities
  for (i=0; i<npointsr; i++) {
    rho = 0.0; 
    for (nu=1; nu<npoints; nu++)
      rho += anu[nu]*besselj0(qnu[nu]*i*deltar);

    fprintf(fp, "%6.3f\t%15.8g\n", i*deltar, rho);
  }
}  


void writeTransitionChargeFourierBessel(FILE* fp,
					const Projection* P,
					const FormfactorPara* p,
					double rmax, int npointsr,
					int Z, int l,
					int jfin, int pfin, int afin,
					int jini, int pini, int aini,
					void* ffactorme,
					const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipjini = idxpij(jmax,pini,jini);
  int ipjfin = idxpij(jmax,pfin,jfin);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idxfin = E->index[ipjfin][afin];
  int idxini = E->index[ipjini][aini];
  double normfin = E->norm[ipjfin][idxfin];
  double normini = E->norm[ipjini][idxini];
  complex double (****ffctme)[2*npoints] = ffactorme;
  complex double (*ff)[npoints] = 
    ffctme[ipjfin][ipjini][idxfin][idxini];
  double anu[npoints];

  int nu;
  complex double ff0, phase;
  double qnu[npoints];
  double Rc = M_PI*(npoints-1)/qmax;

  if (l != 0) {
    fprintf(stderr, "writeChargeDensities not yet implemented for l != 0\n");
    exit(-1);
  }

  // get rid of possible phase factor

  int ifix=npoints/3;
  double qfix=qmax*ifix/(npoints-1);

  ff0 = ff[0][ifix]*fproton(qfix)+ff[1][ifix]*fneutron(qfix);
  phase = ff0/cabs(ff0);

  // calculate coefficients anu
  for (nu=1; nu<npoints; nu++) {
    qnu[nu] = M_PI*nu/Rc;
    anu[nu] = 1/sqrt(4*M_PI)*SQR(qnu[nu])/(2*M_PI*Rc)*
      creal((ff[0][nu]*fproton(qnu[nu])+ff[1][nu]*fneutron(qnu[nu]))/phase)/
      sqrt(normfin*normini);
  }

  // write fourier bessel coefficients

  fprintf(fp, "# Rc      = %8.3f\n", Rc);
  fprintf(fp, "# qmax    = %8.3f\n", qmax);
  fprintf(fp, "# npoints = %d\n", npoints-1);

  fprintf(fp, "\n#nu\tq[nu]\ta[nu]\n");

  for (nu=1; nu<npoints; nu++)
    fprintf(fp, "%d\t%8.3f\t%15.8g\n", nu, qnu[nu], anu[nu]);

}  


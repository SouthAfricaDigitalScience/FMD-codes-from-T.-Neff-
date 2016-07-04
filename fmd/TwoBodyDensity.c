/**

  \file TwoBodyDensity.c

  calculate matrix elements of two-body density operators


  (c) 2010 Thomas Neff

*/

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "TwoBodyDensity.h"

#include "numerics/cmath.h"
#include "numerics/sphericalbessel.h"
#include "numerics/legendrep.h"
#include "numerics/clebsch.h"

#include "misc/physics.h"


#define MAX(x,y) ((x) > (y) ? (x) : (y))

inline double sqr(double x) { return x*x; }

// (S,T)=(0,0),(0,1),(1,0),(1,1)
// summed over MS and MT

ManyBodyOperator OpPairs = {
 name : "Pairs",
 rank : 0,
 pi : 0,
 dim : 4,
 size : 4,
 par : NULL,
 me : calcPairsod
};


// (S,T)=(0,0),(0,1),(1,0),(1,1)
// summed over MS and MT

// we might also be interested in pp,nn channels

ManyBodyOperator OpTwoBodyDensityR = {
  name : NULL,
  rank : 0,
  pi : 0,
  dim : 0,
  size : 0,
  par : NULL,
  me : calcTBDensRod
};

ManyBodyOperator OpTwoBodyDensityQ = {
  name : NULL,
  rank : 0,
  pi : 0,
  dim : 0,
  size : 0,
  par : NULL,
  me : calcTBDensQod
};

// two-body densities in different L channels
// total density for each S,T channel and individual contributions up to Lmax

ManyBodyOperator OpTwoBodyDensityRL = {
  name : NULL,
  rank : 0,
  pi : 0,
  dim : 0,
  size : 0,
  par : NULL,
  me : calcTBDensRLod
};

ManyBodyOperator OpTwoBodyDensityQL = {
  name : NULL,
  rank : 0,
  pi : 0,
  dim : 0,
  size : 0,
  par : NULL,
  me : calcTBDensQLod
};


void initOpTwoBodyDensityR(TBDensRPara* par)
{
  char* name = malloc(60);
  sprintf(name, "TwoBodyDensityR-%05.2f-%d", 
	  par->rmax, par->npoints);
  OpTwoBodyDensityR.name = name;
  OpTwoBodyDensityR.dim = 4*par->npoints;
  OpTwoBodyDensityR.size = 4*par->npoints;
  OpTwoBodyDensityR.par = par;
}  


void initOpTwoBodyDensityQ(TBDensQPara* par)
{
  char* name = malloc(60);
  sprintf(name, "TwoBodyDensityQ-%05.2f-%d", par->qmax, par->npoints);
  OpTwoBodyDensityQ.name = name;
  OpTwoBodyDensityQ.dim = 4*par->npoints;
  OpTwoBodyDensityQ.size = 4*par->npoints;
  OpTwoBodyDensityQ.par = par;
}  
  

void initOpTwoBodyDensityRL(TBDensRLPara* par)
{
  char* name = malloc(60);
  sprintf(name, "TwoBodyDensityRL-%d-%d-%05.2f-%d", 
	  par->lmax, par->lambdamax, par->rmax, par->npoints);
  OpTwoBodyDensityRL.name = name;
  OpTwoBodyDensityRL.dim = (2*par->lmax+6)*par->npoints;
  OpTwoBodyDensityRL.size = (2*par->lmax+6)*par->npoints;
  OpTwoBodyDensityRL.par = par;
}  


void initOpTwoBodyDensityQL(TBDensQLPara* par)
{
  char* name = malloc(60);
  sprintf(name, "TwoBodyDensityQL-%d-%d-%05.2f-%d", 
	  par->lmax, par->lambdamax, par->qmax, par->npoints);
  OpTwoBodyDensityQL.name = name;
  OpTwoBodyDensityQL.dim = (2*par->lmax+6)*par->npoints;
  OpTwoBodyDensityQL.size = (2*par->lmax+6)*par->npoints;
  OpTwoBodyDensityQL.par = par;
}  


static void tb_pairs(void* par,
		     const Gaussian* G1, const Gaussian* G2, 
		     const Gaussian* G3, const Gaussian* G4, 
		     const GaussianAux* X13, const GaussianAux* X24, 
		     complex double pairs[4])
{	
  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
    (G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  complex double SS, sigsig;
  complex double RR;
  
  SS = X13->S * X24->S; 
  sigsig = cvec3mult(X13->sig, X24->sig);

  RR = X13->R * X24->R;

  // S,T=0,0
  pairs[0] += 0.25*(SS-sigsig)*0.25*(TT-tautau)*RR;
  // S,T=0,1
  pairs[1] += 0.25*(SS-sigsig)*0.25*(3*TT+tautau)*RR;
  // S,T=1,0
  pairs[2] += 0.25*(3*SS+sigsig)*0.25*(TT-tautau)*RR;
  // S,T=1,1
  pairs[3] += 0.25*(3*SS+sigsig)*0.25*(3*TT+tautau)*RR;
}

/*

// static workspace for angular grid

static int nalpha = 0;
static int ncosb = 0;
static double* alphal;
static double* walphal;
static double* cosbl;
static double* wcosbl;

static void initializeangulargrid(int na, int nb)
{
  // already initialized ?
  if (nalpha == na && ncosb == nb)
    return;

  nalpha = na;
  ncosb = nb;

  alphal = malloc(nalpha*sizeof(double));
  walphal = malloc(nalpha*sizeof(double));
  cosbl = malloc(ncosb*sizeof(double));
  wcosbl = malloc(ncosb*sizeof(double));

  ShiftedPeriodicTrapezoidalPoints(nalpha, 0, 2*M_PI, 0.0, alphal, walphal);
  GaussLegendrePoints(ncosb, -1, 1, cosbl, wcosbl);
}


static void tb_densr(TBDensRPara* par,
                     const Gaussian* G1, const Gaussian* G2, 
		     const Gaussian* G3, const Gaussian* G4, 
		     const GaussianAux* X13, const GaussianAux* X24, 
		     complex double* densr)
{
  int npoints = par->npoints;
  double rmax= par->rmax;

  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
    (G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  complex double SS, sigsig;
  complex double RR;
  
  SS = X13->S * X24->S; 
  sigsig = cvec3mult(X13->sig, X24->sig);

  RR = X13->R * X24->R;

  complex double alphar, rho[3];
  int i;

  alphar = X13->alpha + X24->alpha;
  for (i=0; i<3; i++)
    rho[i] = X13->rho[i]-X24->rho[i];

  initializeangulargrid(par->nalpha, par->ncosb);
  
  int ir, ialpha, icosb;
  double r, alpha, beta, x[3], weight;
  complex double xmrho[3], xmrho2, dens;

  for (ir=0; ir<npoints; ir++) {
    r = rmax* ir/(npoints-1);

    for (ialpha=0; ialpha<nalpha; ialpha++)
      for (icosb=0; icosb<ncosb; icosb++) {
	alpha = alphal[ialpha];
	beta = acos(cosbl[icosb]);
	weight = walphal[ialpha]*wcosbl[icosb];

	x[0] = r*cos(alpha)*sin(beta);
	x[1] = r*sin(alpha)*sin(beta);
	x[2] = r*cos(beta);

	for (i=0; i<3; i++)
	  xmrho[i] = x[i] - rho[i];
	xmrho2 = cvec3sqr(xmrho);
	
	dens = RR/cpow(2*M_PI*alphar, 1.5)*cexp(-0.5*xmrho2/alphar)* weight;

	// S,T=0,0
	densr[ir          ] += 0.25*(SS-sigsig)*0.25*(TT-tautau)*dens;
	// S,T=0,1
	densr[ir+1*npoints] += 0.25*(SS-sigsig)*0.25*(3*TT+tautau)*dens;
	// S,T=1,0
	densr[ir+2*npoints] += 0.25*(3*SS+sigsig)*0.25*(TT-tautau)*dens;
	// S,T=1,1
	densr[ir+3*npoints] += 0.25*(3*SS+sigsig)*0.25*(3*TT+tautau)*dens;
      }
  }
}

*/

// analytic integration over angles

static void tb_densr(TBDensRPara* par,
                     const Gaussian* G1, const Gaussian* G2, 
		     const Gaussian* G3, const Gaussian* G4, 
		     const GaussianAux* X13, const GaussianAux* X24, 
		     complex double* densr)
{
  double rmax = par->rmax;
  int npoints = par->npoints;

  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
    (G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  complex double SS, sigsig;
  complex double RR;
  
  SS = X13->S * X24->S; 
  sigsig = cvec3mult(X13->sig, X24->sig);

  RR = X13->R * X24->R;

  complex double alphar, rho[3], rho2, beta[3], beta2;
  int i;

  alphar = X13->alpha + X24->alpha;
  for (i=0; i<3; i++)
    rho[i] = X13->rho[i]-X24->rho[i];
  rho2 = cvec3sqr(rho);
  for (i=0; i<3; i++)
    beta[i] = rho[i]/alphar;
  beta2 = cvec3sqr(beta);

  complex double dens;

  int ir;
  double r;
  for (ir=0; ir<npoints; ir++) {
    r = rmax* ir/(npoints-1);

    dens = RR/cpow(2*M_PI*alphar, 1.5)*cexp(-0.5/alphar*(r*r+rho2))* 
      4*M_PI*cbesseli0(r*csqrt(beta2));

    // S,T=0,0
    densr[ir          ] += 0.25*(SS-sigsig)*0.25*(TT-tautau)*dens;
    // S,T=0,1
    densr[ir+1*npoints] += 0.25*(SS-sigsig)*0.25*(3*TT+tautau)*dens;
    // S,T=1,0
    densr[ir+2*npoints] += 0.25*(3*SS+sigsig)*0.25*(TT-tautau)*dens;
    // S,T=1,1
    densr[ir+3*npoints] += 0.25*(3*SS+sigsig)*0.25*(3*TT+tautau)*dens;
  }
}


static void tb_densq(TBDensQPara* par,
                     const Gaussian* G1, const Gaussian* G2, 
		     const Gaussian* G3, const Gaussian* G4, 
		     const GaussianAux* X13, const GaussianAux* X24, 
		     complex double* densq)
{
  int npoints = par->npoints;
  double qmax= par->qmax;

  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
    (G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  complex double SS, sigsig;
  complex double RR;
  
  SS = X13->S * X24->S; 
  sigsig = cvec3mult(X13->sig, X24->sig);

  RR = X13->R * X24->R;

  complex double alphaq, pi[3], pi2, beta[3], beta2;
  int i;

  alphaq = 4.0/(X13->lambda + X24->lambda);
  for (i=0; i<3; i++)
    pi[i] = 0.5*(X13->pi[i] - X24->pi[i]);
  pi2 = cvec3sqr(pi);
  for (i=0; i<3; i++)
    beta[i] = alphaq*pi[i];
  beta2 = cvec3sqr(beta);
    
  complex double dens;
  
  int iq;
  double q;
  for (iq=0; iq<npoints; iq++) {
    q = qmax* iq/(npoints-1);

    dens = RR*cpow(alphaq/(2*M_PI), 1.5)*cexp(-0.5*alphaq*(q*q+pi2))*
      4*M_PI*cbesseli0(q*csqrt(beta2));

    // S,T=0,0
    densq[iq          ] += 0.25*(SS-sigsig)*0.25*(TT-tautau)*dens;
    // S,T=0,1
    densq[iq+1*npoints] += 0.25*(SS-sigsig)*0.25*(3*TT+tautau)*dens;
    // S,T=1,0
    densq[iq+2*npoints] += 0.25*(3*SS+sigsig)*0.25*(TT-tautau)*dens;
    // S,T=1,1
    densq[iq+3*npoints] += 0.25*(3*SS+sigsig)*0.25*(3*TT+tautau)*dens;
  }
}


// static workspace for clebsch gordans
// cg2[lmax+1][lambdamax+1][lmax+lambdamax+1]

double* cg2tab = NULL;

void initcg2(int lmax, int lambdamax)
{
  if (cg2tab) return;

  cg2tab = malloc((lmax+1)*(lambdamax+1)*(lmax+lambdamax+1)*sizeof(double));
  double (*cg2)[lambdamax+1][lmax+lambdamax+1] = cg2tab;
  
  int l, lambda, lp;
  for (l=0; l<=lmax; l++)
    for (lambda=0; lambda<=lambdamax; lambda++)
      for (lp=abs(l-lambda); lp<=l+lambda; lp++)
	cg2[l][lambda][lp] = sqr(clebsch(2*l,2*lp,2*lambda,0,0,0));
}

	    
// bug for lambdamax > 0: number of pairs is too large !?

static void tb_densrl(TBDensRLPara* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const Gaussian* G3, const Gaussian* G4, 
		      const GaussianAux* X13, const GaussianAux* X24, 
		      complex double* densrl)
{
  double rmax= par->rmax;
  int npoints = par->npoints;
  int lmax = par->lmax;
  int lambdamax = par->lambdamax;

  initcg2(lmax, lambdamax);
  double (*cg2)[lambdamax+1][lmax+lambdamax+1] = cg2tab;

  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
    (G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  double PiT0, PiT1;
  PiT0 = 0.25*(TT-tautau);
  PiT1 = 0.25*(3*TT+tautau);

  complex double SS, sigsig;
  complex double PiS0, PiS1;
  SS = X13->S * X24->S; 
  sigsig = cvec3mult(X13->sig, X24->sig);
  PiS0 = 0.25*(SS-sigsig); 
  PiS1 = 0.25*(3*SS+sigsig);

  complex double RR;
  RR = X13->R * X24->R;

  complex double alphar, rho[3], rho2, pi[3];
  int i;

  alphar = 1.0/(X13->alpha + X24->alpha);
  for (i=0; i<3; i++)
    rho[i] = X13->rho[i]-X24->rho[i];
  rho2 = cvec3sqr(rho);
  for (i=0; i<3; i++)
    pi[i] = 0.5*(X13->pi[i]-X24->pi[i]);

  complex double betav[3], beta;
  for (i=0; i<3; i++)
    betav[i] = alphar*rho[i];
  beta = csqrt(cvec3sqr(betav));

  complex double lambda14, lambda23; 
  lambda14 = 1.0/(conj(G1->a)+G4->a);
  lambda23 = 1.0/(conj(G2->a)+G3->a);

  complex double Arp, App;
  Arp = (conj(G1->a)-G4->a)*lambda14+(conj(G2->a)-G3->a)*lambda23;
  App = 2.0*(conj(G1->a)*G4->a*lambda14+conj(G2->a)*G3->a*lambda23);
  
  complex double Axy;
  complex double betaxv[3], betayv[3], betax, betay, betaxy;
  Axy = -0.25*alphar+0.5/App;
  for (i=0; i<3; i++) {
    betaxv[i] = (0.5*alphar-0.5*Arp/App)*rho[i] - I*pi[i] ;
    betayv[i] = (0.5*alphar+0.5*Arp/App)*rho[i] + I*pi[i] ;
  }
  betax = csqrt(cvec3sqr(betaxv));
  betay = csqrt(cvec3sqr(betayv));
  // if length of vectors is zero set cosine to 1
  if (cabs(betax) < 1.E-8 || cabs(betay) < 1.E-8) {
    betaxy = 1.0;
  } else {
    betaxy = cvec3mult(betaxv,betayv)/(betax*betay);
  }

  complex double dens, densl[lmax+1];
  complex double ilx[lmax+lambdamax+1], ily[lmax+lambdamax+1], 
    ilxy[lambdamax+1], Lpxy[lmax+lambdamax+1];

  clegendreps(lmax+lambdamax, betaxy, Lpxy);

  int ir;
  double r;
  complex double Cr;
  for (ir=0; ir<npoints; ir++) {
    r = rmax* ir/(npoints-1);
    
    // common factor
    Cr = RR*cpow(alphar/(2*M_PI), 1.5)*cexp(-0.5*alphar*(r*r+rho2));

    // total density 
    dens = Cr* 4*M_PI*cbesseli0(beta*r);

    // densities in individual channels
    cbesselis(lmax+lambdamax, betax*r, ilx);
    cbesselis(lmax+lambdamax, betay*r, ily);
    cbesselis(lambdamax, Axy*sqr(r), ilxy);

    int l, lambda, lp;
    for (l=0; l<=lmax; l++) {
      densl[l] = 0.0;
      for (lambda=0; lambda<=lambdamax; lambda++)
	for (lp=abs(l-lambda); lp<=l+lambda; lp++)
	  densl[l] += Cr* 4*M_PI* (2*l+1)*(2*lp+1)*
	    ilx[lp]*ily[lp]*Lpxy[lp]*
	    ilxy[lambda]*cg2[l][lambda][lp];
    }

    densrl[ir          ] += PiS0*PiT0*dens;
    densrl[ir+1*npoints] += PiS0*PiT1*dens;
    densrl[ir+2*npoints] += PiS1*PiT0*dens;
    densrl[ir+3*npoints] += PiS1*PiT1*dens;

    int idx=4;
    for (l=0; l<=lmax; l++)
      if (l%2) { // odd channels
	densrl[ir+idx*npoints] += PiS0*PiT0*densl[l]; idx++;
	densrl[ir+idx*npoints] += PiS1*PiT1*densl[l]; idx++;
      } else { // even channels
	densrl[ir+idx*npoints] += PiS0*PiT1*densl[l]; idx++;
	densrl[ir+idx*npoints] += PiS1*PiT0*densl[l]; idx++;
      }
  }
}


static void tb_densql(TBDensQLPara* par,
		      const Gaussian* G1, const Gaussian* G2, 
		      const Gaussian* G3, const Gaussian* G4, 
		      const GaussianAux* X13, const GaussianAux* X24, 
		      complex double* densql)
{
  double qmax= par->qmax;
  int npoints = par->npoints;
  int lmax = par->lmax;
  int lambdamax = par->lambdamax;

  initcg2(lmax, lambdamax);
  double (*cg2)[lambdamax+1][lmax+lambdamax+1] = cg2tab;

  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+
    (G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  double PiT0, PiT1;
  PiT0 = 0.25*(TT-tautau);
  PiT1 = 0.25*(3*TT+tautau);

  complex double SS, sigsig;
  complex double PiS0, PiS1;
  SS = X13->S * X24->S; 
  sigsig = cvec3mult(X13->sig, X24->sig);
  PiS0 = 0.25*(SS-sigsig); 
  PiS1 = 0.25*(3*SS+sigsig);

  complex double RR;
  RR = X13->R * X24->R;

  complex double alphaq, rho[3], pi[3], pi2;
  int i;

  alphaq = 4.0/(X13->lambda+X24->lambda);
  for (i=0; i<3; i++)
    rho[i] = X13->rho[i]-X24->rho[i];
  for (i=0; i<3; i++)
    pi[i] = 0.5*(X13->pi[i]-X24->pi[i]);
  pi2 = cvec3sqr(pi);

  complex double betav[3], beta;
  for (i=0; i<3; i++)
    betav[i] = alphaq*pi[i];
  beta = csqrt(cvec3sqr(betav));

  complex double lambda14, lambda23; 
  lambda14 = 1.0/(conj(G1->a)+G4->a);
  lambda23 = 1.0/(conj(G2->a)+G3->a);

  complex double Arr, Arp;
  Arr = 0.5*(lambda14+lambda23);
  Arp = (conj(G1->a)-G4->a)*lambda14+(conj(G2->a)-G3->a)*lambda23;
  
  complex double Axy;
  complex double betaxv[3], betayv[3], betax, betay, betaxy;
  Axy = -0.25*alphaq+0.5/Arr;
  for (i=0; i<3; i++) {
    betaxv[i] = (0.5*alphaq+0.5*Arp/Arr)*pi[i] + I*rho[i] ;
    betayv[i] = (0.5*alphaq-0.5*Arp/Arr)*pi[i] - I*rho[i] ;
  }
  betax = csqrt(cvec3sqr(betaxv));
  betay = csqrt(cvec3sqr(betayv));
  // if length of vectors is zero set cosine to 1
  if (cabs(betax) < 1.E-8 || cabs(betay) < 1.E-8) {
    betaxy = 1.0;
  } else {
    betaxy = cvec3mult(betaxv,betayv)/(betax*betay);
  }

  complex double dens, densl[lmax+1];
  complex double ilx[lmax+lambdamax+1], ily[lmax+lambdamax+1], 
    ilxy[lambdamax+1], Lpxy[lmax+lambdamax+1];

  clegendreps(lmax+lambdamax, betaxy, Lpxy);

  int iq;
  double q;
  complex double Cq;
  for (iq=0; iq<npoints; iq++) {
    q = qmax* iq/(npoints-1);
    
    // common factor
    Cq = RR*cpow(alphaq/(2*M_PI), 1.5)*cexp(-0.5*alphaq*(q*q+pi2));

    // total density 
    dens = Cq* 4*M_PI*cbesseli0(beta*q);

    // densities in individual channels
    cbesselis(lmax+lambdamax, betax*q, ilx);
    cbesselis(lmax+lambdamax, betay*q, ily);
    cbesselis(lambdamax, Axy*sqr(q), ilxy);
   
    int l, lambda, lp;
    for (l=0; l<=lmax; l++) {
      densl[l] = 0.0;
      for (lambda=0; lambda<=lambdamax; lambda++)
	for (lp=abs(l-lambda); lp<=l+lambda; lp++)
	  densl[l] += Cq* 4*M_PI* (2*l+1)*(2*lp+1)*
	    ilx[lp]*ily[lp]*Lpxy[lp]*
	    ilxy[lambda]*cg2[l][lambda][lp];
    }

    densql[iq          ] += PiS0*PiT0*dens;
    densql[iq+1*npoints] += PiS0*PiT1*dens;
    densql[iq+2*npoints] += PiS1*PiT0*dens;
    densql[iq+3*npoints] += PiS1*PiT1*dens;

    int idx=4;
    for (l=0; l<=lmax; l++)
      if (l%2) { // odd channels
	densql[iq+idx*npoints] += PiS0*PiT0*densl[l]; idx++;
	densql[iq+idx*npoints] += PiS1*PiT1*densl[l]; idx++;
      } else { // even channels
	densql[iq+idx*npoints] += PiS0*PiT1*densl[l]; idx++;
	densql[iq+idx*npoints] += PiS1*PiT0*densl[l]; idx++;
      }
  }
}


void calcPairs(const SlaterDet* Q, const SlaterDetAux* X, double pairs[4])
{
  TwoBodyOperator op_tb_pairs = {dim: 4, opt: 0, par: NULL, me: tb_pairs};

  calcSlaterDetTBME(Q, X, &op_tb_pairs, pairs);
}


void calcPairsod(void* Par,
		 const SlaterDet* Q, const SlaterDet* Qp,
		 const SlaterDetAux* X, complex double pairs[4])
{
  TwoBodyOperator op_tb_pairs = {dim: 4, opt: 0, par: NULL, me: tb_pairs};

  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_pairs, pairs);
}


void writeprojectedPairs(FILE* fp, const Projection* P, 
			 const complex double (**pairsexp)[4],
			 const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i;

  char prefix[8];

  int* idx;
  int ngood;
  complex double *norm, *H;
  complex double (*pairs)[4];
  

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

        if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
        else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

	idx = E->index[idxpij(jmax,p,j)];
        norm = E->norm[idxpij(jmax,p,j)];
        H = E->v[idxpij(jmax,p,j)];
 
        pairs = pairsexp[idxpij(jmax,p,j)];

        fprintf(fp, "\n%s N           = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
        fprintf(fp, "\n%s H           = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));
        fprintf(fp, "\n%s (S,T)=(0,0) = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.3f", creal(pairs[idx[i]][0]/norm[idx[i]]));
        fprintf(fp, "\n%s (S,T)=(0,1) = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.3f", creal(pairs[idx[i]][1]/norm[idx[i]]));
        fprintf(fp, "\n%s (S,T)=(1,0) = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.3f", creal(pairs[idx[i]][2]/norm[idx[i]]));
        fprintf(fp, "\n%s (S,T)=(1,1) = ", prefix); 
        for (i=0; i<ngood; i++) 
          fprintf(fp, "   %8.3f", creal(pairs[idx[i]][3]/norm[idx[i]]));

	fprintf(fp, "\n");        
      }
    }
}


void calcTBDensRod(TBDensRPara* par,
                   const SlaterDet* Q, const SlaterDet* Qp,
                   const SlaterDetAux* X, complex double* densr)
{
  int npoints=par->npoints;

  TwoBodyOperator op_tb_densr = {dim: 4*npoints, opt: 0, par: par, me: tb_densr};

  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_densr, densr);
}


void writeTBDensR(FILE* fp,
		  const Projection* P,
		  const TBDensRPara* p,
		  int j, int pi, int a,
		  void* tbdensme,
		  const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double rmax = p->rmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**densme)[4*npoints] = tbdensme;
  complex double (*dens)[npoints] = densme[ipj][idx];
  int i;
  double r, d00, d01, d10, d11;

  for (i=0; i<npoints; i++) {
    r = rmax*i/(npoints-1);
    d00 = dens[0][i]/norm;
    d01 = dens[1][i]/norm;
    d10 = dens[2][i]/norm;
    d11 = dens[3][i]/norm;

    fprintf(fp, "%6.3f\t%13.8g\t%13.8g\t%13.8g\t%13.8g\n",
	    r, d00, d01, d10, d11);
  }
}


void calcTBDensQod(TBDensQPara* par,
                   const SlaterDet* Q, const SlaterDet* Qp,
                   const SlaterDetAux* X, complex double* densq)
{
  int npoints=par->npoints;

  TwoBodyOperator op_tb_densq = {dim: 4*npoints, opt: 0, par: par, me: tb_densq};

  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_densq, densq);
}


void writeTBDensQ(FILE* fp,
		  const Projection* P,
		  const TBDensQPara* p,
		  int j, int pi, int a,
		  void* tbdensme,
		  const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**densme)[4*npoints] = tbdensme;
  complex double (*dens)[npoints] = densme[ipj][idx];
  int i;
  double q, d00, d01, d10, d11;

  for (i=0; i<npoints; i++) {
    q = qmax*i/(npoints-1);
    d00 = dens[0][i]/norm;
    d01 = dens[1][i]/norm;
    d10 = dens[2][i]/norm;
    d11 = dens[3][i]/norm;

    fprintf(fp, "%6.3f\t%13.8g\t%13.8g\t%13.8g\t%13.8g\n",
	    q, d00, d01, d10, d11);
  }
}


void calcTBDensRLod(TBDensRLPara* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X, complex double* densr)
{
  int npoints=par->npoints;
  int dim=2*par->lmax+6;

  TwoBodyOperator op_tb_densrl = {dim: dim*npoints, opt: 0, par: par, me: tb_densrl};

  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_densrl, densr);
}


void writeTBDensRL(FILE* fp,
		   const Projection* P,
		   const TBDensRLPara* p,
		   int j, int pi, int a,
		   void* tbdensme,
		   const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  int lmax = p->lmax;
  int dim = 2*lmax+6;
  double rmax = p->rmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**densme)[dim*npoints] = tbdensme;
  complex double (*dens)[npoints] = densme[ipj][idx];
  int i, l, d;
  double r, dd;

  for (i=0; i<npoints; i++) {
    r = rmax*i/(npoints-1);
    fprintf(fp, "  %6.3f", r);
    for (d=0; d<dim; d++) {
      dd = dens[d][i]/norm;
      fprintf(fp, "\t%13.8g", dd);
    }
    fprintf(fp, "\n");
  }
}


void calcTBDensQLod(TBDensQLPara* par,
		    const SlaterDet* Q, const SlaterDet* Qp,
		    const SlaterDetAux* X, complex double* densq)
{
  int npoints=par->npoints;
  int dim=2*par->lmax+6;

  TwoBodyOperator op_tb_densql = {dim: dim*npoints, opt: 0, par: par, me: tb_densql};

  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_densql, densq);
}


void writeTBDensQL(FILE* fp,
		   const Projection* P,
		   const TBDensQLPara* p,
		   int j, int pi, int a,
		   void* tbdensme,
		   const Eigenstates* E)
{
  int jmax = P->jmax;
  int ipj = idxpij(jmax,pi,j);
  int npoints = p->npoints;
  int lmax = p->lmax;
  int dim = 2*lmax+6;
  double qmax = p->qmax;
  int idx = E->index[ipj][a];
  complex double norm = E->norm[ipj][idx];
  complex double (**densme)[dim*npoints] = tbdensme;
  complex double (*dens)[npoints] = densme[ipj][idx];
  int i, l, d;
  double q, dd;

  for (i=0; i<npoints; i++) {
    q = qmax*i/(npoints-1);
    fprintf(fp, "  %6.3f", q);
    for (d=0; d<dim; d++) {
      dd = dens[d][i]/norm;
      fprintf(fp, "\t%13.8g", dd);
    }
    fprintf(fp, "\n");
  }
}

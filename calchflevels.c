/**

  \file calcHFlevels.c

  calculate the HF single-particle states of a FMD Slater determinant


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/CenterofMass.h"
#include "fmd/Potential.h"
#include "fmd/Hamiltonian.h"
#include "fmd/AngularMomenta.h"
#include "fmd/Isospin.h"
#include "fmd/Parity.h"
#include "fmd/NOscillator.h"

#include "numerics/cmat.h"
#include "misc/utils.h"
#include "misc/physics.h"


void expect(int n, complex double* A, complex double* V, double* a)
{
  int i;
  int k,m;

  for (i=0; i<n; i++) {
    a[i] = 0.0;
    for (m=0; m<n; m++)
      for (k=0; k<n; k++)
	a[i] += conj(V[k+i*n])*A[k+m*n]*V[m+i*n];
  }
}


double angroot(double j2)
{
  return (-0.5+sqrt(j2+0.25));
}


typedef struct {
  double rank;
  double e;
  double j2;
  double l2;
  double hosci;
  double t3;
  double pi;
} SingleParticleState;


int cmpSingleParticleState(SingleParticleState* spa, SingleParticleState* spb)
{
  if (spa->rank < spb->rank)
    return +1;
  else
    return -1;
}


int main(int argc, char *argv[])
{
  createinfo(argc, argv);
  
  int c;
  int nljorder=0;
  int potentialonly=0;
  int cm=0;
  double omega=0.0;

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction slaterdetfile\n"
	    "\n   -s           sort sp states according to n, l, j"
	    "\n   -v           use only potential"
	    "\n   -o OMEGA     use oscillator constant [MeV]"
	    "\n   -c           T-Tcm\n",
	    argv[0]);
    exit(-1);
  }
  
  while((c = getopt(argc, argv, "svco:")) != -1)
    switch (c) {
    case 's':
      nljorder = 1;
      break;
    case 'v':
      potentialonly=1;
      break;
    case 'c':
      cm = 1;
      break;
    case 'o':
      omega = atof(optarg)/hbc;
      break;
    }
  
  char* interactionfile = argv[optind];
  char* slaterdetfile = argv[optind+1];

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    exit(-1);
  Int.cm = cm;

  SlaterDet Q;
  if (readSlaterDetfromFile(&Q, slaterdetfile))
    exit(-1);

  // number of nucleons
  int A = Q.A;

  // single-particle overlap matrix
  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  // Tcm
  double tcm;
  calcTCM(&Q, &X, &tcm);

  // derive oscillator constant from center of mass motion
  double omegacm = 4.0/3.0*tcm;

  if (omega == 0.0)
    omega = omegacm;

  complex double n[A*A];
  copycmat(A, X.n, n);

  // hartree-fock matrix in Gaussian single-particle basis

  complex double hhf[A*A];
  if (potentialonly)	calcPotentialHF(&Int, &Q, &X, hhf);
  else			calcHamiltonianHF(&Int, &Q, &X, hhf);
  
  complex double l2hf[A*A];
  calcl2HF(&Q, &X, l2hf);
  
  complex double j2hf[A*A];
  calcj2HF(&Q, &X, j2hf);

  complex double hoscihf[A*A];
  calcHOsciHF(&Q, &X, omega, hoscihf);

  complex double t3hf[A*A];
  calct3HF(&Q, &X, t3hf);

  complex double parhf[A*A];
  calcparHF(&Q, &X, parhf);

  // solve eigenvalue problem
  complex double lambda[A];
  complex double V[A*A];
  int dim;

  if (nljorder) {
    complex double nljhf[A*A];
    int i;
    for (i=0; i<A*A; i++)
      nljhf[i] = 10000*hoscihf[i]/omega - 100*j2hf[i] + l2hf[i];

      generalizedeigensystem(nljhf, n, A, 0.0, lambda, V, &dim);
  } else
    generalizedeigensystem(hhf, n, A, 0.0, lambda, V, &dim);

  // calculate expectation values
  double norm[A];
  expect(A, n, V, norm);

  double h[A];
  expect(A, hhf, V, h);

  double l2[A];
  expect(A, l2hf, V, l2);
  
  double j2[A];
  expect(A, j2hf, V, j2);

  double hosci[A];
  expect(A, hoscihf, V, hosci);

  double t3[A];
  expect(A, t3hf, V, t3);

  double par[A];
  expect(A, parhf, V, par);

  fprintinfo(stdout);

  // sort single-particle states by energy
  SingleParticleState sp[A];

  int i;
  for (i=0; i<A; i++) {
    sp[i].rank = lambda[i];
    sp[i].e = h[i]/norm[i];
    sp[i].j2 = j2[i]/norm[i];
    sp[i].l2 = l2[i]/norm[i];
    sp[i].hosci = hosci[i]/norm[i];
    sp[i].t3 = t3[i]/norm[i];
    sp[i].pi = par[i]/norm[i];
  }	
  
  qsort(sp, A, sizeof(SingleParticleState), cmpSingleParticleState);
  
  // print proton and neutron levels

  fprintf(stdout, "\nusing oscillator constant: hbar Omega = %8.3f MeV\n",
	  omega*hbc);

  fprintf(stdout, "\nproton levels:\n");
  for (i=0; i<A; i++)
    if (sp[i].t3 > 0.0) {
      fprintf(stdout, "e: %8.3f MeV, j: %5.3f, l: %5.3f, pi: %+5.2f, nosci: %5.3f\n",
	      hbc*sp[i].e, angroot(sp[i].j2), angroot(sp[i].l2), sp[i].pi,
	      sp[i].hosci/omega-1.5);
    }

  fprintf(stdout, "\nneutron levels:\n");
  for (i=0; i<A; i++)
    if (sp[i].t3 < 0.0) {
      fprintf(stdout, "e: %8.3f MeV, j: %5.3f, l: %5.3f, pi: %+5.2f, nosci: %5.3f\n",
	      hbc*sp[i].e, angroot(sp[i].j2), angroot(sp[i].l2), sp[i].pi,
	      sp[i].hosci/omega-1.5);
    }

  return 0;
}

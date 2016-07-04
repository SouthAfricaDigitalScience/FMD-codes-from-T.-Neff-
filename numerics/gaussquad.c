/**

  \file gaussquad.c

  Knots and Weights for Gaussian Quadrature
  calling routines from IQPACK

  (c) 2005 Thomas Neff

*/


#include <stdio.h>
#include <math.h>
#include <float.h>
#include "fortranc.h"

#include "gaussquad.h"

static inline double dabs(double x)
{
  return (x < 0 ? -x : x);
}

static inline double dsqr(double x)
{
  return (x*x);
}


void FORTRAN(cgqf)(int* nt, double* t, double* wts, int* kind,
		   double* alpha, double* beta, double* a, double* b,
		   int* lo, int* nwf, double* wf, int* niwf, int* iwf,
		   int* ier);


void FORTRAN(machep)(double* x)
{
  *x = DBL_EPSILON;
}


/// weight function 1

void GaussLegendrePoints(int n, double a, double b, 
			 double x[], double w[])
{
  int KIND=1;
  double T[n];
  double WTS[n];
  int LO=0;
  int NWF=2*n;
  double WF[NWF];
  int NIWF=2*n;
  int IWF[NIWF];
  int IER;

  FORTRAN(cgqf)(&n, T, WTS, &KIND, NULL, NULL, &a, &b, &LO,
		&NWF, WF, &NIWF, IWF, &IER);

  int i;
  for (i=0; i<n; i++) {
    x[i] = T[i];
    w[i] = WTS[i];
  }

  if (IER) {
    fprintf(stderr, "CGQF error code: %d\n", IER);
  }
  
}


/// Divide weights by Chebyshev weight function

void GaussChebyshevPoints(int n, double a, double b, 
                          double x[], double w[])
{
  int KIND=2;
  double T[n];
  double WTS[n];
  int LO=0;
  int NWF=2*n;
  double WF[NWF];
  int NIWF=2*n;
  int IWF[NIWF];
  int IER;

  FORTRAN(cgqf)(&n, T, WTS, &KIND, NULL, NULL, &a, &b, &LO,
		&NWF, WF, &NIWF, IWF, &IER);

  int i;
  for (i=0; i<n; i++) {
    x[i] = T[i];
    w[i] = WTS[i]/(pow((b-x[i])*(x[i]-a), -0.5));
  }

  if (IER) {
    fprintf(stderr, "CGQF error code: %d\n", IER);
  }
  
}


/// Divide weights by Jacobi weight function

void GaussJacobiPoints(int n, double a, double b, 
			     double alpha, double beta,
			     double x[], double w[])
{
  int KIND=4;
  double T[n];
  double WTS[n];
  int LO=0;
  int NWF=2*n;
  double WF[NWF];
  int NIWF=2*n;
  int IWF[NIWF];
  int IER;

  FORTRAN(cgqf)(&n, T, WTS, &KIND, &alpha, &beta, &a, &b, &LO,
		&NWF, WF, &NIWF, IWF, &IER);

  int i;
  for (i=0; i<n; i++) {
    x[i] = T[i];
    w[i] = WTS[i]/(pow(b-x[i], alpha)*pow(x[i]-a, beta));
  }

  if (IER) {
    fprintf(stderr, "CGQF error code: %d\n", IER);
  }
  
}


/// Divide weights by Exponential weight function

void GaussExponentialPoints(int n, double a, double b, 
			    double alpha,
			    double x[], double w[])
{
  int KIND=7;
  double T[n];
  double WTS[n];
  int LO=0;
  int NWF=2*n;
  double WF[NWF];
  int NIWF=2*n;
  int IWF[NIWF];
  int IER;

  FORTRAN(cgqf)(&n, T, WTS, &KIND, &alpha, NULL, &a, &b, &LO,
		&NWF, WF, &NIWF, IWF, &IER);

  int i;
  for (i=0; i<n; i++) {
    x[i] = T[i];
    // don't use inf for x = (a+b)/2
    w[i] = (dabs(x[i]-0.5*(a+b))<1e-9 ? 0.0 : WTS[i]/(pow(dabs(x[i]-0.5*(a+b)), alpha)));
  }

  if (IER) {
    fprintf(stderr, "CGQF error code: %d\n", IER);
  }
  
}

#define GPOINTS 6

const double rroots[GPOINTS][GPOINTS] = {
  { 1.1283791670955126, 0.0, 0.0, 0.0, 0.0, 0.0 },
  { 0.7539869079898889, 1.734055298879165, 0.0, 0.0, 0.0, 0.0 },
  { 0.5507554932505945, 1.2873090298702823, 2.220484543409835, 0.0, 0.0, 0.0 },
  { 0.4238628193901623, 1.0143321045669047, 1.7424373751621773, 2.639813333635729, 0.0, 0.0 },
  { 0.3384095960706952, 0.8266625437756976, 1.4328543726759924, 2.145408222424107, 3.0141724674934722, 0.0 },
  { 0.27778943815369034, 0.6899111472703859, 1.208987335931981, 1.813832458324578, 2.510418043632598, 3.355569460516944 }};

const double rweights[GPOINTS][GPOINTS] = {
  { 1.2432708030495752, 0.0, 0.0, 0.0, 0.0, 0.0 },
  { 0.8504841015767846, 1.138571648468604, 0.0, 0.0, 0.0, 0.0 },
  { 0.644019215153679, 0.8243863791634771, 1.0777615246427668, 0.0, 0.0, 0.0 },
  { 0.5095472549555211, 0.6622918325848185, 0.7976347900889224, 1.0344858511856345, 0.0,  0.0 },
  { 0.4150111323268987, 0.5527176371663822, 0.6576037440453348, 0.774073858814447, 1.0009856479417851, 0.0 },
  { 0.3456032994608327, 0.47120445227330776, 0.5633988416389077, 0.6470640287975326, 0.7537960884832555, 0.9737666201169521 }};


void GaussSqrHermitePoints(int n, double b, double x[], double w[])
{
  int ir;

  if (n>GPOINTS) {
    fprintf(stderr, "GaussSqrHermitePoints not yet implemented for n>%d\n", GPOINTS);
    return;
  }

  for (ir=0; ir<n; ir++) {
    x[ir] = rroots[n-1][ir]/sqrt(b);
    w[ir] = dsqr(rroots[n-1][ir])/b* rweights[n-1][ir]/sqrt(b);
  }
}


/// Divide weights by Gen. Hermite weight function
/// integrate from 0 to infinity, therefore take only half of the points

void GaussGenHermitePoints(int n, double b, 
			   double x[], double w[])
{
  int KIND=6;
  int N=2*n;
  double A=0.0;
  double ALPHA=2.0;
  double T[N];
  double WTS[N];
  int LO=0;
  int NWF=2*N;
  double WF[NWF];
  int NIWF=2*N;
  int IWF[NIWF];
  int IER;

  FORTRAN(cgqf)(&N, T, WTS, &KIND, &ALPHA, NULL, &A, &b, &LO,
		&NWF, WF, &NIWF, IWF, &IER);

  int i;
  for (i=0; i<n; i++) {
    x[i] = T[n+i];
    w[i] = WTS[n+i]/(pow(x[i],ALPHA)*exp(-b*dsqr(x[i])));
  }

  if (IER) {
    fprintf(stderr, "CGQF error code: %d\n", IER);
  }
  
}
  

void ExtendedTrapezoidalPoints(int n, double a, double b,
			       double x[], double w[])
{
  int i;
  double h=(b-a)/(n-1);
  
  x[0] = a; x[n-1] = b;
  w[0] = 0.5*h; w[n-1] = 0.5*h;

  for (i=1; i<n-1; i++) {
    x[i] = a+i*h;
    w[i] = h;
  }
}


void ShiftedPeriodicTrapezoidalPoints(int n, double a, double b, double shift,
				      double x[], double w[])
{
  int i;
  double h = (b-a)/n;

  for (i=0; i<n; i++) {
    x[i] = a+ (i+shift)*h;
    w[i] = h;
  }
}

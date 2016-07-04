/**

  \file gaussquad.h

  Knots and Weights for Gaussian Quadrature
  calling routines from IQPACK

  (c) 2005 Thomas Neff

*/


#ifndef _GAUSSQUAD_H
#define _GAUSSQUAD_H


void GaussLegendrePoints(int n, double a, double b, 
			 double x[], double w[]);

void GaussChebyshevPoints(int n, double a, double b, 
                          double x[], double w[]);

void GaussJacobiPoints(int n, double a, double b, 
		       double alpha, double beta,
		       double x[], double w[]);

void GaussExponentialPoints(int n, double a, double b, 
			    double alpha,
			    double x[], double w[]);

void GaussSqrHermitePoints(int n, double b,
			   double x[], double w[]);

void GaussGenHermitePoints(int n, double b,
			   double x[], double w[]);

void ExtendedTrapezoidalPoints(int n, double a, double b,
			       double x[], double w[]);

void ShiftedPeriodicTrapezoidalPoints(int n, double a, double b, double shift,
				      double x[], double w[]);

#endif
  
  
  

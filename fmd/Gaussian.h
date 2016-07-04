/**

  \file Gaussian.h

  Parameters of a Gaussian Wavepacket


  (c) 2003 Thomas Neff

*/


#ifndef _GAUSSIAN_H
#define _GAUSSIAN_H

#include <complex.h>


typedef struct {
  int xi;			///< isospin +/- 1 (proton/neutron)	
  complex double chi[2];	///< complex spinor, additional freedom to give norm and phase
  complex double a;		///< complex width	
  complex double b[3];		///< complex position in phase space
} Gaussian;


/// Auxiliary quantities for matrix elements with Gaussians.
/// Definitions as given in the diploma thesis
typedef struct {
  complex double sig[3];
  complex double lambda;
  complex double alpha;
  complex double rho[3];
  complex double rho2;
  complex double pi[3];
  complex double pi2;
  complex double rhopi;
  complex double rhoxpi[3];
  int	         T;
  complex double S;
  complex double R;
  complex double Q;
} GaussianAux;


/// One-body operator.
/// generic definition of an operator that calculates matrix elements
/// with Gaussian one-body states. matrix elements will be added to 
/// val[dim]. If the matrix element is proportional to the isospin
/// overlap T opt should be set to true.
typedef struct {
  int dim;
  int opt;
  void* par;
  void (*me)(void* parameters,
	     const Gaussian* G1, const Gaussian* G2, const GaussianAux* X12, 
	     complex double val[]);
} OneBodyOperator;


/// Two-body operator.
/// generic definition of an operator that calculates matrix elements
/// with Gaussian one-body states. matrix elements will be added to 
/// val[dim]. If the matrix element is proportional to the isospin
/// overlap T opt should be set to true.
typedef struct {
  int dim;
  int opt;
  void* par;
  void (*me)(void* parameters,
	     const Gaussian* G1, const Gaussian* G2, 
	     const Gaussian* G3, const Gaussian* G4,
	     const GaussianAux* X13, const GaussianAux* X24,
	     complex double val[]);
} TwoBodyOperator;


/// calculate auxiliary quantities for two Gaussians
void calcGaussianAux(const Gaussian* ga, const Gaussian* gb, GaussianAux* aux);


/// move Gaussian by d
void moveGaussian(Gaussian* g, double d[3]);

/// boost Gaussian by v
void boostGaussian(Gaussian* g, double v[3]);

/// rotate Gaussian using rotation matrices R3 in coordinate and
/// R2 in spin space 
void rotateGaussian(Gaussian* G, 
		    const double R3[3][3], const complex double R2[2][2]);

/// parity operation with respect to origin
void invertGaussian(Gaussian* g);

/// time reversal operator
void timerevertGaussian(Gaussian* g);

/// scale lengths by kappa
void scaleGaussian(Gaussian* g, double kappa);

/// spin flip
void spinflipGaussian(Gaussian* g);


#endif

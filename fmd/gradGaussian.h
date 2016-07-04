/**

  \file gradGaussian.h

  Derivatives of a Gaussian Wavepacket


  (c) 2003 Thomas Neff

*/


#ifndef _GRADGAUSSIAN_H
#define _GRADGAUSSIAN_H

#include <complex.h>

#include "Gaussian.h"
#include "../numerics/rotationmatrices.h"


/// gradient of a spinor
typedef struct {
  complex double chi[2];
} gradSpinor;

/// vector of spinor gradients
typedef struct {
  complex double chi[2][3];
} gradVectorSpinor;

/// gradient of a scalar quantity
typedef struct {
  complex double a;
  complex double b[3];
} gradScalar;

/// gradient of a vector quantity
typedef struct {
  complex double a[3];
  complex double b;
} gradVector;


/// gradient with respect to spin and spatial parameters
typedef struct {
  complex double chi[2];	///< derivatives with respect to spinor
  complex double a;		///< derivative with respect to a
  complex double b[3];		///< derivatives with respect to b
} gradGaussian;


/// Derivatives of spin and spatial overlap of Gaussians
typedef struct {
  gradVectorSpinor dsig;
  complex double dlambda;
  complex double dalpha;
  gradVector drho;
  gradScalar drho2;
  gradVector dpi;
  gradScalar dpi2;

  gradSpinor dS;	     	///< derivatives of spin overlap S
  gradScalar dR;		///< derivatives of spatial overlap R
  gradGaussian dQ;		///< derivatives of overlap Q
} gradGaussianAux;


/// One-body operator
typedef struct {
  int opt;
  void* par;
  void (*me)(void* parameters,
	     const Gaussian* G1, const Gaussian* G2, 
	     const GaussianAux* X12, const gradGaussianAux* dX12, 
	     complex double* val, gradGaussian* dval);
} gradOneBodyOperator;


/// Two-body operator
typedef struct {
  int opt;
  void* par;
  void (*me)(void* parameters,
	     const Gaussian* G1, const Gaussian* G2, 
	     const Gaussian* G3, const Gaussian* G4,
	     const GaussianAux* X13, const GaussianAux* X24,
	     const gradGaussianAux* dX13,
	     complex double* val, gradGaussian* dval);
} gradTwoBodyOperator;


///
inline static void zerogradGaussian(gradGaussian* dG)
{
  int i;

  for (i=0; i<2; i++)
    dG->chi[i] = 0.0;
  dG->a = 0.0;
  for (i=0; i<3; i++)
    dG->b[i] = 0.0;
}

inline static void addtogradGaussian(gradGaussian* dG, const gradGaussian* ddG)
{
  int i;

  for (i=0; i<2; i++)
    dG->chi[i] += ddG->chi[i];
  dG->a += ddG->a;
  for (i=0; i<3; i++)
    dG->b[i] += ddG->b[i];
}

inline static void multgradGaussian(gradGaussian* dG, complex double x)
{
  int i;

  for (i=0; i<2; i++)
    dG->chi[i] *= x;
  dG->a *= x;
  for (i=0; i<3; i++)
    dG->b[i] *= x;
}

inline static void addmulttogradGaussian(gradGaussian* dG, const gradGaussian* ddG, const complex double x)
{
  int i;

  for (i=0; i<2; i++)
    dG->chi[i] += x*ddG->chi[i];
  dG->a += x*ddG->a;
  for (i=0; i<3; i++)
    dG->b[i] += x*ddG->b[i];
}

inline static void rotategradGaussian(gradGaussian* dG, 
				      const double R3[3][3], const complex double R2[2][2])
{
  cm2mult(R2, dG->chi);
  m3mult(R3, dG->b);
}


/// calculate auxiliary quantities for two Gaussians
void calcgradGaussianAux(const Gaussian* G1, const Gaussian* G2, 
			 const GaussianAux* X, gradGaussianAux* dX);


#endif

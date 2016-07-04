/*

  cmath.h

  collection of complex algebra routines


  (c) 2003 Thomas Neff

*/


#ifndef _CALGEBRA_H
#define _CALGEBRA_H


#include <stdio.h>
#include <complex.h>


/// complex x^2
inline static complex double csqr(complex double x)
{
  return(x*x);
}

/// complex x^3
inline static complex double cpow3(complex double x)
{
  return(x*x*x);
}

/// complex x^3/2
inline static complex double cpow32(complex double x)
{
  return(x*csqrt(x));
}

/// vector squared
inline static double vec3sqr(const double v[3])
{
  return(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

/// complex vector squared
inline static complex double cvec3sqr(const complex double v[3])
{
  return(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

/// scalar product
inline static double vec3mult(const double a[3],
                         const double b[3])
{
  return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

/// complex scalar product
inline static complex double cvec3mult(const complex double a[3], 
			 const complex double b[3])
{
  return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

/// cross product c = a x b
inline static void vec3cross(const double a[3],
			     const double b[3], double c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/// complex cross product c = a x b
inline static void cvec3cross(const complex double a[3], 
			      const complex double b[3], complex double c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/// spat product (axb).c
inline static double vec3spat(const double a[3],
			      const double b[3],
			      const double c[3])
{
  return ((a[1]*b[2]-a[2]*b[1])*c[0] +
          (a[2]*b[0]-a[0]*b[2])*c[1] +
          (a[0]*b[1]-a[1]*b[0])*c[2]);
}

/// complex spat product (axb).c
inline static complex double cvec3spat(const complex double a[3], 
				       const complex double b[3], 
				       const complex double c[3])
{
  return ((a[1]*b[2]-a[2]*b[1])*c[0] +
	  (a[2]*b[0]-a[0]*b[2])*c[1] +
	  (a[0]*b[1]-a[1]*b[0])*c[2]);
}


#endif

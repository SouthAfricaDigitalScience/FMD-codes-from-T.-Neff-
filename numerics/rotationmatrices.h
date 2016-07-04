/**

  \file rotationmatrices.h

  Rotation matrices in Spinor and Coordinate space
  and derivatives of rotation matrices with respect to Euler angles.


  (c) 2004 Thomas Neff

*/


#ifndef _ROTATIONMATRICES_H
#define _ROTATIONMATRICES_H


#include <math.h>
#include <complex.h>


/// multiply complex 3-vector with real 3x3-matrix
void m3mult(const double A[3][3], complex double x[3]);


/// multiply complex 2-vector with complex 2x2-matrix
void cm2mult(const complex double A[2][2], complex double x[2]);


/// spinor rotation matrix
void rotatemat2(const double euler[3], complex double rmat2[2][2]);


/// derivatives of spinor rotation matrix
void derivrotatemat2(const double euler[3], complex double drmat2[3][2][2]);

/// vector rotation matrix
void rotatemat3(const double euler[3], double rmat3[3][3]);


/// derivatives of vector rotation matrix
void derivrotatemat3(const double euler[3], double drmat3[3][3][3]);


#endif

/**

  \file rotationmatrices.c

  Rotation matrices in Spinor and Coordinate space
  and derivatives of rotation matrices with respect to Euler angles.


  (c) 2004 Thomas Neff

*/


#include <math.h>
#include <complex.h>

#include "rotationmatrices.h"


void m3mult(const double A[3][3], complex double x[3])
{
  complex double y[3];
  int i,j;

  for (i=0; i<3; i++) {
    y[i] = 0.0;
    for (j=0; j<3; j++) 
      y[i] += A[i][j]*x[j];
  }
  for (i=0; i<3; i++)
    x[i] = y[i];
}


void cm2mult(const complex double A[2][2], complex double x[2])
{
  complex double y[2];
  int i,j;

  for (i=0; i<2; i++) {
    y[i] = 0.0;
    for (j=0; j<2; j++) 
      y[i] += A[i][j]*x[j];
  }
  for (i=0; i<2; i++)
    x[i] = y[i];
}


void rotatemat2(const double euler[3], complex double rmat2[2][2])
{
  double a = euler[0]; double b = euler[1]; double g = euler[2];
  rmat2[0][0] = cos(0.5*b)*cexp(-I*0.5*(a+g));
  rmat2[0][1] = -sin(0.5*b)*cexp(-I*0.5*(a-g));
  rmat2[1][0] = sin(0.5*b)*cexp(I*0.5*(a-g));
  rmat2[1][1] = cos(0.5*b)*cexp(I*0.5*(a+g));
}


void derivrotatemat2(const double euler[3], complex double drmat2[3][2][2])
{
  double a = euler[0]; double b = euler[1]; double g = euler[2];

  // derivative with respect to a
  drmat2[0][0][0] = -0.5*I*cos(0.5*b)*cexp(-I*0.5*(a+g));
  drmat2[0][0][1] = +0.5*I*sin(0.5*b)*cexp(-I*0.5*(a-g));
  drmat2[0][1][0] = +0.5*I*sin(0.5*b)*cexp(I*0.5*(a-g));
  drmat2[0][1][1] = +0.5*I*cos(0.5*b)*cexp(I*0.5*(a+g));

  // derivative with respect to b
  drmat2[1][0][0] = -0.5*sin(0.5*b)*cexp(-I*0.5*(a+g));
  drmat2[1][0][1] = -0.5*cos(0.5*b)*cexp(-I*0.5*(a-g));
  drmat2[1][1][0] = +0.5*cos(0.5*b)*cexp(I*0.5*(a-g));
  drmat2[1][1][1] = -0.5*sin(0.5*b)*cexp(I*0.5*(a+g));

  // derivative with respect to g
  drmat2[2][0][0] = -0.5*I*cos(0.5*b)*cexp(-I*0.5*(a+g));
  drmat2[2][0][1] = -0.5*I*sin(0.5*b)*cexp(-I*0.5*(a-g));
  drmat2[2][1][0] = -0.5*I*sin(0.5*b)*cexp(I*0.5*(a-g));
  drmat2[2][1][1] = +0.5*I*cos(0.5*b)*cexp(I*0.5*(a+g));
}


void rotatemat3(const double euler[3], double rmat3[3][3])
{
  double a = euler[0]; double b = euler[1]; double g = euler[2];

  rmat3[0][0] = cos(g)*cos(b)*cos(a)-sin(g)*sin(a);
  rmat3[0][1] = -sin(g)*cos(b)*cos(a)-cos(g)*sin(a);
  rmat3[0][2] = sin(b)*cos(a);
  rmat3[1][0] = cos(g)*cos(b)*sin(a)+sin(g)*cos(a);
  rmat3[1][1] = -sin(g)*cos(b)*sin(a)+cos(g)*cos(a);
  rmat3[1][2] = sin(b)*sin(a);
  rmat3[2][0] = -cos(g)*sin(b);
  rmat3[2][1] =  sin(g)*sin(b);
  rmat3[2][2] = cos(b);
}


void derivrotatemat3(const double euler[3], double drmat3[3][3][3])
{
  double a = euler[0]; double b = euler[1]; double g = euler[2];

  // derivative with respect to a
  drmat3[0][0][0] = -cos(g)*cos(b)*sin(a)-sin(g)*cos(a);
  drmat3[0][0][1] = sin(g)*cos(b)*sin(a)-cos(g)*cos(a);
  drmat3[0][0][2] = -sin(b)*sin(a);
  drmat3[0][1][0] = cos(g)*cos(b)*cos(a)-sin(g)*sin(a);
  drmat3[0][1][1] = -sin(g)*cos(b)*cos(a)-cos(g)*sin(a);
  drmat3[0][1][2] = sin(b)*cos(a);
  drmat3[0][2][0] = 0.0;
  drmat3[0][2][1] = 0.0;
  drmat3[0][2][2] = 0.0;

  // derivative with respect to b
  drmat3[1][0][0] = -cos(g)*sin(b)*cos(a);
  drmat3[1][0][1] = sin(g)*sin(b)*cos(a);
  drmat3[1][0][2] = cos(b)*cos(a);
  drmat3[1][1][0] = -cos(g)*sin(b)*sin(a);
  drmat3[1][1][1] = sin(g)*sin(b)*sin(a);
  drmat3[1][1][2] = cos(b)*sin(a);
  drmat3[1][2][0] = -cos(g)*cos(b);
  drmat3[1][2][1] =  sin(g)*cos(b);
  drmat3[1][2][2] = -sin(b);

  // derivative with respect to g
  drmat3[2][0][0] = -sin(g)*cos(b)*cos(a)-cos(g)*sin(a);
  drmat3[2][0][1] = -cos(g)*cos(b)*cos(a)+sin(g)*sin(a);
  drmat3[2][0][2] = 0.0;
  drmat3[2][1][0] = -sin(g)*cos(b)*sin(a)+cos(g)*cos(a);
  drmat3[2][1][1] = -cos(g)*cos(b)*sin(a)-sin(g)*cos(a);
  drmat3[2][1][2] = 0.0;
  drmat3[2][2][0] = sin(g)*sin(b);
  drmat3[2][2][1] = cos(g)*sin(b);
  drmat3[2][2][2] = 0.0;
}

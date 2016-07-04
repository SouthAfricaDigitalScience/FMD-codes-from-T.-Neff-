/**

  \file Quadrupole.c

  calculate quadrupole parameters beta and gamma


  (c) 2003 Thomas Neff

*/

#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Quadrupole.h"
#include "Radii.h"

#include "numerics/cmath.h"


static void ob_quadrupole(const SlaterDet* Q,
			  const Gaussian*G1, const Gaussian* G2, 
			  const GaussianAux* X, 
			  complex double q[9][3])
{
  int A=Q->A, Z=Q->Z, N=Q->N;
  int i,j;
  complex double qme;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      qme = 3*X->rho[i]*X->rho[j]*X->Q;
      q[i+j*3][0] += (1.0-1.0/A)* qme;
      q[i+j*3][1] += ((1.0-2.0/A)*(1+G1->xi)/2 + 1.0*Z/(A*A))* qme;
      q[i+j*3][2] += ((1.0-2.0/A)*(1-G1->xi)/2 + 1.0*N/(A*A))* qme;
    }
    
    qme = X->rho2* X->Q;
    q[j+j*3][0] += -(1.0-1.0/A)* qme;
    q[j+j*3][1] += -(1.0*Z/(A*A) + (1.0-2.0/A)*(1+G1->xi)/2)* qme;
    q[j+j*3][2] += -(1.0*N/(A*A) + (1.0-2.0/A)*(1-G1->xi)/2)* qme;
  }
}


static void tb_quadrupole(const SlaterDet* Q,
			  const Gaussian* G1, const Gaussian* G2, 
			  const Gaussian* G3, const Gaussian* G4, 
			  const GaussianAux* X13, const GaussianAux* X24, 
			  complex double q[9][3])
{
  int A=Q->A, Z=Q->Z, N=Q->N;
  int i,j;
  complex double qme;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++) {
      qme = 3.0*(X13->rho[i]*X24->rho[j]+X13->rho[j]*X24->rho[i])* X13->Q* X24->Q;
      q[i+j*3][0] += -1.0/A* qme;
      q[i+j*3][1] += (1.0*Z/(A*A) -1.0/A*((1+G1->xi)/2+(1+G2->xi)/2))* qme;
      q[i+j*3][2] += (1.0*N/(A*A) -1.0/A*((1-G1->xi)/2+(1-G2->xi)/2))* qme;
    }

    qme = 2.0* cvec3mult(X13->rho, X24->rho)* X13->Q* X24->Q;
    q[j+j*3][0] += 1.0/A* qme;
    q[j+j*3][1] += -(1.0*Z/(A*A) - 1.0/A*((1+G1->xi)/2+(1+G2->xi)/2))* qme;
    q[j+j*3][2] += -(1.0*N/(A*A) - 1.0/A*((1-G1->xi)/2+(1-G2->xi)/2))* qme;
  }
}

#define SQR(x) (x)*(x)

static double det(double R[3][3])
{
  return (R[0][0]*R[1][1]*R[2][2] + R[0][1]*R[1][2]*R[2][0] + 
          R[0][2]*R[1][0]*R[2][1] - R[0][0]*R[1][2]*R[2][1] - 
          R[0][1]*R[1][0]*R[2][2] - R[0][2]*R[1][1]*R[2][0]);
}


void calcQuadrupole(const SlaterDet* Q, const SlaterDetAux* X, 
		    double beta[3], double gamma[3])
{
  double qone[9][3], qtwo[9][3], q[9], qsq, detq;
  int i,j;

  OneBodyOperator op_ob_quadrupole = {dim: 27, opt: 1, par: Q, me: ob_quadrupole};
  TwoBodyOperator op_tb_quadrupole = {dim: 27, opt: 1, par: Q, me: tb_quadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_quadrupole, qone);
  calcSlaterDetTBME(Q, X, &op_tb_quadrupole, qtwo);

  for (j=0; j<3; j++) {
    qsq = 0.0;
    for (i=0; i<9; i++) {
      q[i] = qone[i][j] + qtwo[i][j];
      // q[i] = qone[i][j];
      qsq += SQR(q[i]);
    }
    detq = det(q);

    beta[j] = sqrt(5/(24*M_PI)*qsq);
    gamma[j] = 60/M_PI*acos(5*sqrt(5)/(16*pow(M_PI,1.5))*detq/pow(beta[j],3));
  }
}


void calcQuadrupolePrime(const SlaterDet* Q, const SlaterDetAux* X, 
			 double beta[3])
{
  double qone[9][3], qtwo[9][3], q[9], qsq, detq;
  double r2[3];
  int i,j;

  calcRadii2(Q, X, &r2[0], &r2[1], &r2[2]);
  
  OneBodyOperator op_ob_quadrupole = {dim: 27, opt: 1, par: Q, me: ob_quadrupole};
  TwoBodyOperator op_tb_quadrupole = {dim: 27, opt: 1, par: Q, me: tb_quadrupole};

  calcSlaterDetOBME(Q, X, &op_ob_quadrupole, qone);
  calcSlaterDetTBME(Q, X, &op_tb_quadrupole, qtwo);

  for (j=0; j<3; j++) {
    qsq = 0.0;
    for (i=0; i<9; i++) {
      q[i] = qone[i][j] + qtwo[i][j];
      // q[i] = qone[i][j];
      qsq += SQR(q[i]);
    }
    detq = det(q);

    beta[j] = sqrt(2*M_PI/15*qsq)/r2[j];
  }
  beta[0] /= Q->A;
  beta[1] /= Q->Z;
  beta[2] /= Q->N;

}

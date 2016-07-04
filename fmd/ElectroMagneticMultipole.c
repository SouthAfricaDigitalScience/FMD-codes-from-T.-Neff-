/**

  \file ElectroMagneticMultipole.c 

  Electro-Magnetic Multipole Operators

  calculates electric monopole, magnetic dipole and
  electric quadrupole moment in spherical basis


  (c) 2003 Thomas Neff

*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Projection.h"
#include "ElectroMagneticMultipole.h"

#include "numerics/cmath.h"
#include "misc/utils.h"
#include "misc/physics.h"


// EMonopole
ManyBodyOperator OpEMonopole = {
  name : "EMonopole",
  rank : 0,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcEMonopoleod
};

// EDipole
ManyBodyOperator OpEDipole = {
  name : "EDipole",
  rank : 2,
  pi : 1,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcEDipoleod
};

// MDipole
ManyBodyOperator OpMDipole = {
  name : "MDipole",
  rank : 2,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcMDipoleod
};

// EQuadrupole
ManyBodyOperator OpEQuadrupole = {
  name : "EQuadrupole",
  rank : 4,
  pi : 0,
  dim : 1,
  size : 1,
  par : NULL,
  me : calcEQuadrupoleod
};


// identical to charge radius squared (modulo factor Z)

static void ob_emonopole(const SlaterDet* Q,
			 const Gaussian* G1, const Gaussian* G2, 
			 const GaussianAux* X, complex double *emonopole)
{
  complex double mr2;
  int A=Q->A, Z=Q->Z;

  mr2 = (3.0*X->alpha + X->rho2)* X->Q;
 
  *emonopole += (1.0*Z/(A*A) + (1+G1->xi)/2 *(1.0-2.0/A))* mr2;
}


static void tb_emonopole(const SlaterDet* Q,
			 const Gaussian* G1, const Gaussian* G2, 
			 const Gaussian* G3, const Gaussian* G4, 
			 const GaussianAux* X13, const GaussianAux* X24, 
			 complex double *emonopole)
{	
  complex double mr2;
  int A=Q->A, Z=Q->Z;

  mr2 = cvec3mult(X13->rho, X24->rho)*X13->Q*X24->Q;

  *emonopole +=  2.0*(1.0*Z/(A*A) - 
		      ((1+G1->xi)/2+(1+G2->xi)/2) *1.0/A)* mr2;
}


static void ob_edipole(const SlaterDet* Q,
		       const Gaussian* G1, const Gaussian* G2,
		       const GaussianAux* X, complex double ed[3])
{
  int A=Q->A, Z=Q->Z;
  double M = Z*mproton+(A-Z)*mneutron;
  int i;

  for (i=0; i<3; i++)	
    ed[i] += ((1+G1->xi)/2 - Z*mass(G1->xi)/M)* X->rho[i]* X->Q;
}


static void ob_mdipole(const SlaterDet* Q,
		       const Gaussian* G1, const Gaussian* G2, 
		       const GaussianAux* X, complex double m[3])
{
  int A=Q->A, Z=Q->Z;
  complex double muspin[3], muorbital[3];
  int i;

  for (i=0; i<3; i++)
    muspin[i] = X->T*((1+G1->xi)/2*muproton + (1-G1->xi)/2*muneutron)*
      X->sig[i]*X->R;
  
  for (i=0; i<3; i++)
    muorbital[i] = ((1+G1->xi)/2*(1.0-2.0/A)+1.0*Z/(A*A))* X->rhoxpi[i]* X->Q;
	
  for (i=0; i<3; i++)
    m[i] += muspin[i] + muorbital[i]; 
}


static void tb_mdipole(const SlaterDet* Q,
		       const Gaussian* G1, const Gaussian* G2, 
		       const Gaussian* G3, const Gaussian* G4, 
		       const GaussianAux* X13, const GaussianAux* X24, 
		       complex double m[3])
{
  int A=Q->A, Z=Q->Z;
  int i;

  for (i=0; i<3; i++)
    m[i] += 2*(-((1+G1->xi)/2+(1+G2->xi)/2)*1.0/A + 1.0*Z/(A*A))*
      (X13->rho[(i+1)%3]*X24->pi[(i+2)%3]-X13->rho[(i+2)%3]*X24->pi[(i+1)%3])*
      X13->Q* X24->Q;
}


static void ob_equadrupole(const SlaterDet* Q,
			   const Gaussian* G1, const Gaussian* G2, 
			   const GaussianAux* X, 
			   complex double q[9])
{
  int A=Q->A, Z=Q->Z;
  int i,j;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++)
      q[i+j*3] += ((1.0-2.0/A)*(1+G1->xi)/2 + 1.0*Z/(A*A))*
	X->rho[i]*X->rho[j]*X->Q;

    // q[j+j*3] += - 1.0/3.0* ((1.0-2.0/A)*(1+G1->xi)/2 + 1.0*Z/(A*A))*
    //  X->rho2* X->Q; 
  }
}


static void tb_equadrupole(const SlaterDet* Q,
			   const Gaussian* G1, const Gaussian* G2, 
			   const Gaussian* G3, const Gaussian* G4, 
			   const GaussianAux* X13, const GaussianAux* X24, 
			   complex double q[9])
{
  int A=Q->A, Z=Q->Z;
  int i,j;

  for (j=0; j<3; j++) {
    for (i=0; i<3; i++)
      q[i+j*3] += (-2.0/A*(1+G1->xi)/2 + 1.0*Z/(A*A))*
	(X13->rho[i]*X24->rho[j]+X13->rho[j]*X24->rho[i])* X13->Q* X24->Q;

    //  q[j+j*3] += -2.0/3.0* (-2.0/A*(1+G1->xi)/2 + 1.0*Z/(A*A))* 
    //  cvec3mult(X13->rho, X24->rho)* X13->Q* X24->Q;
  }
}


void calcEMonopoleod(void* Par,
		     const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X, 
		     complex double* emonopole)
{
  complex double emonoone, emonotwo;

  OneBodyOperator op_ob_emono = {dim: 1, opt: 1, par: Q, me: ob_emonopole};
  TwoBodyOperator op_tb_emono = {dim: 1, opt: 1, par: Q, me: tb_emonopole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_emono, &emonoone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_emono, &emonotwo);

  *emonopole = emonoone + emonotwo;
}


void calcEDipoleod(void* Par,
		   const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double edipole[3])
{
  complex double ed[3];

  OneBodyOperator op_ob_edipole = {dim: 3, opt: 1, par: Q, me: ob_edipole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_edipole, ed);

  // spherical components -1, 0, +1
  edipole[0] = 	 0.5*sqrt(1.5/M_PI)*(ed[0] - I*ed[1]);
  edipole[1] =   0.5*sqrt(3.0/M_PI)*(ed[2]);
  edipole[2] = - 0.5*sqrt(1.5/M_PI)*(ed[0] + I*ed[1]);
}



void calcMDipoleod(void* Par,
		   const SlaterDet* Q, const SlaterDet* Qp,
		   const SlaterDetAux* X, 
		   complex double mdipole[3])
{
  complex double mone[3], mtwo[3], m[3];
  int i;

  OneBodyOperator op_ob_mdi = {dim: 3, opt: 1, par: Q, me: ob_mdipole};
  TwoBodyOperator op_tb_mdi = {dim: 3, opt: 1, par: Q, me: tb_mdipole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_mdi, mone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_mdi, mtwo);

  for (i=0; i<3; i++)
    m[i] = mone[i] + mtwo[i];

  // spherical components -1, 0, +1
  mdipole[0] = 	 0.5*sqrt(1.5/M_PI)*(m[0] - I*m[1]);
  mdipole[1] =   0.5*sqrt(3.0/M_PI)*(m[2]);
  mdipole[2] = - 0.5*sqrt(1.5/M_PI)*(m[0] + I*m[1]);
}


void calcEQuadrupoleod(void* Par,
		       const SlaterDet* Q, const SlaterDet* Qp,
		       const SlaterDetAux* X, 
		       complex double equadrupole[5])
{
  complex double qone[9], qtwo[9], q[9];
  int i;

  OneBodyOperator op_ob_equad = {dim: 9, opt: 1, par: Q, me: ob_equadrupole};
  TwoBodyOperator op_tb_equad = {dim: 9, opt: 1, par: Q, me: tb_equadrupole};

  calcSlaterDetOBMEod(Q, Qp, X, &op_ob_equad, qone);
  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_equad, qtwo);

  for (i=0; i<9; i++)
    q[i] = qone[i] + qtwo[i];

  // spherical components -2, -1, 0, 1, 2
  equadrupole[0] =  0.25*sqrt(7.5/M_PI)*(q[0]-I*(q[1]+q[3])-q[4]);
  equadrupole[1] =  0.25*sqrt(7.5/M_PI)*(q[2]+q[6]-I*(q[5]+q[7]));
  equadrupole[2] = -0.25*sqrt(5.0/M_PI)*(q[0]+q[4]-2*q[8]);
  equadrupole[3] = -0.25*sqrt(7.5/M_PI)*(q[2]+q[6]+I*(q[5]+q[7]));
  equadrupole[4] =  0.25*sqrt(7.5/M_PI)*(q[0]+I*(q[1]+q[3])-q[4]);
}
					 

void writeprojectedEMMultipoles(FILE* fp,
				const Projection* P,
				const complex double** emonopole,
				const complex double** edipole,
				const complex double** mdipole,
				const complex double** equadrupole,
				const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int p,j,i;

  char prefix[8];

  int* idx;
  int ngood;
  complex double* norm;
  complex double *H, *E0, *M1, *E2;
  

  for (p=0; p<=1; p++)
    for (j=odd; j<jmax; j=j+2) {
      
      ngood = E->ngood[idxpij(jmax,p,j)];

      if (ngood) {

	if(odd) sprintf(prefix, "[%d/2%c]", j, p ? '-' : '+'); 
	else    sprintf(prefix, "[%d%c]", j/2, p ? '-' : '+'); 

	idx = E->index[idxpij(jmax,p,j)];
	norm = E->norm[idxpij(jmax,p,j)];
	H = E->v[idxpij(jmax,p,j)];
 
	E0 = emonopole[idxpij(jmax,p,j)];
	M1 = mdipole[idxpij(jmax,p,j)];
	E2 = equadrupole[idxpij(jmax,p,j)];

	fprintf(fp, "\n%s N     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.5f", creal(norm[idx[i]]));
	fprintf(fp, "\n%s H     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", hbc*creal(H[idx[i]]));
	fprintf(fp, "\n%s mu    = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", 
		  sqrt(4.0*M_PI/3.0)*sqrt(1.0*j/(j+2))*
		  creal(M1[idx[i]]/norm[idx[i]]));
	fprintf(fp, "\n%s Q     = ", prefix); 
	for (i=0; i<ngood; i++) 
	  fprintf(fp, "   %8.3f", 
		  sqrt(16.0*M_PI/5.0)*sqrt(1.0*j*(j-1)/((j+3)*(j+2)))*
		  creal(E2[idx[i]]/norm[idx[i]]));

	fprintf(fp, "\n");	  
      }
    }
}


#define SQR(x) ((x)*(x))


void writeprojectedtransitions(FILE* fp,
			       const Projection* P,
			       int rank, int pi,
			       const char* label, const char* unit,
			       const complex double**** transme,
			       const Eigenstates* E)
{
  int odd=P->odd;
  int jmax=P->jmax;

  int pini, jini, pfin, jfin;
  int idxini, idxfin;
  int *idxi, *idxf;
  int i, f;
  double B;

  for (pini=0; pini<=1; pini++)
    for (jini=odd; jini<jmax; jini=jini+2) {
      pfin = (pini+pi)%2;
      for (jfin=abs(jini-rank); jfin<=min(jini+rank,jmax-1); jfin=jfin+2) {

	idxini = idxpij(jmax,pini,jini);
	idxfin = idxpij(jmax,pfin,jfin);

	if (E->ngood[idxfin] == 0 || 
	    E->ngood[idxini] == 0)
	  continue;

	idxi = E->index[idxini];
	idxf = E->index[idxfin];

	if (P->odd) fprintf(fp, "\n\n %s [%d/2%c] <- [%d/2%c] (%s)\n\n",
			    label, jfin, pfin ? '-' : '+', 
			    jini, pini ? '-' : '+', unit);
	else	    fprintf(fp, "\n\n %s [%d%c] <- [%d%c] (%s)\n\n",
			    label, jfin/2, pfin ? '-' : '+', 
			    jini/2, pini ? '-' : '+', unit);

	fprintf(fp, "%14c", ' ');
	for (i=0; i<E->ngood[idxini]; i++) 
	  fprintf(fp, "%10.2f MeV", hbc*creal(E->v[idxini][idxi[i]]));

	for (f=0; f<E->ngood[idxfin]; f++) {
	  fprintf(fp, "\n%10.2f MeV", hbc*creal(E->v[idxfin][idxf[f]]));
	  for (i=0; i<E->ngood[idxini]; i++) {
	    if (rank == 0) {
	      B = sqrt(jfin+1.0)*
		cabs(transme[idxfin][idxini][idxf[f]][idxi[i]])/
		sqrt(cabs(E->norm[idxfin][idxf[f]]*E->norm[idxini][idxi[i]]));
	    } else {	
	      B = (jfin+1.0)/(jini+1.0)*
		SQR(cabs(transme[idxfin][idxini][idxf[f]][idxi[i]]))/
		cabs(E->norm[idxfin][idxf[f]]*E->norm[idxini][idxi[i]]);
	    }
	    if (creal(E->v[idxfin][idxf[f]]) < creal(E->v[idxini][idxi[i]]))
	      fprintf(fp, "%10.3f    ", B);
	    else
	      fprintf(fp, "%10.3f *  ", B);
	  }
	}
	fprintf(fp, "\n");	
      }
    }
}


void writeprojectedtransitionEMMultipoles(FILE* fp,
					  const Projection* P,
					  const complex double**** emonopole,
					  const complex double**** edipole,
					  const complex double**** mdipole,
					  const complex double**** equadrupole,
					  const Eigenstates* E)
{
  fprintf(fp, "\n ########### Electric Monopole Transitions\n");
  writeprojectedtransitions(fp, P, 0, 0, "M(E0)", "fm^2", emonopole, E);

  fprintf(fp, "\n ########### Electric Dipole Transitions\n");
  writeprojectedtransitions(fp, P, 2, 1, "B(E1)", "e^2 fm^2", edipole, E);

  fprintf(fp, "\n ########### Magnetic Dipole Transitions\n");
  writeprojectedtransitions(fp, P, 2, 0, "B(M1)", "mu0^2", mdipole, E);

  fprintf(fp, "\n ########### Electric Quadrupole Transitions\n");
  writeprojectedtransitions(fp, P, 4, 0, "B(E2)", "e^2 fm^4", equadrupole, E);

}

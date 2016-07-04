/**

  \file Potential.c

  calculate matrix element of Interaction


  (c) 2003 Thomas Neff

  Changes (HH):
  06/08/03
  + Added factor 0.5 to the gradients of prVTrp implicit in the prescription
    used for the super operator basis (hermitization).
  
  04/12/21
  + Added matrix element of {L2 S12(p,p)}_H
  + Introduced new auxiliaries: psi, thelakapi 
  
  04/09/03
  + Added matrix elements for S12(p,p) and (p_r v(r) + v(r)p_r)S12(r,p)

  04/07/20
  + Added matrix elements of L2, L2LS, and S12(L,L) interaction terms.
  + Slightly reorganized calculation of auxiliaries.

*/

#include <complex.h>

#include "Gaussian.h"
#include "SlaterDet.h"
#include "Interaction.h"
#include "Potential.h"

#include "numerics/cmath.h"
#include "numerics/coulomb.h"


static void tb_pot(Interaction* P,
		   const Gaussian* G1, const Gaussian* G2, 
		   const Gaussian* G3, const Gaussian* G4, 
		   const GaussianAux* X13, const GaussianAux* X24, 
		   complex double v[])
{	
  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+(G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!TT && !tautau) return;

  complex double SS, RR;
  complex double lambda, alpha, psi;
  complex double rho[3], rho2, pi[3], pi2, sigsig;
  complex double rhoxpi[3], S[3];
  complex double L2, LS, S12, beta, theta, rhopi;
  complex double S12pipi, S12ll, S12rhopi;
  
  int i;
  
  SS = X13->S * X24->S; RR = X13->R * X24->R;
  lambda = X13->lambda + X24->lambda;
  alpha = X13->alpha + X24->alpha;

  for (i=0; i<3; i++) rho[i] = X13->rho[i] - X24->rho[i];
  rho2 = cvec3sqr(rho);
  for (i=0; i<3; i++) pi[i] = 0.5*(X13->pi[i] - X24->pi[i]);

  sigsig = cvec3mult(X13->sig, X24->sig);

  if (P->l2 || P->spinorbit || P->l2ls || P->tll || P->tpp || P->l2tpp){
    cvec3cross(rho, pi, rhoxpi);
    L2 = cvec3sqr(rhoxpi);
  }
  
  if (P->spinorbit || P->l2ls) {
    for (i=0; i<3; i++)
      S[i] = 0.5*(X13->sig[i]*X24->S+X13->S*X24->sig[i]);
    LS = cvec3mult(rhoxpi, S);
  }

  if (P->momentump2 || P->momentumpr2 || P->l2 || P->l2ls || P->tll || P->tpp || P->prtrp || P->l2tpp) {
    beta = I*((conj(G1->a) - G3->a)*X13->lambda + (conj(G2->a) - G4->a)*X24->lambda);
    theta = (conj(G1->a)*X13->lambda + conj(G2->a)*X24->lambda)*
            (G3->a*X13->lambda + G4->a*X24->lambda);
    psi = theta - alpha*lambda;	    
    pi2 = cvec3sqr(pi);
    rhopi = cvec3mult(rho, pi);
  }

  if (P->tensor || P->tll || P->tpp || P->prtrp || P->l2tpp)
    S12 = 3*cvec3mult(X13->sig, rho)*cvec3mult(X24->sig, rho) - sigsig*rho2;

  if (P->tll || P->tpp || P->prtrp || P->l2tpp) {
    S12rhopi = 1.5*(cvec3mult(X13->sig, rho)*cvec3mult(X24->sig, pi) + cvec3mult(X13->sig, pi)*cvec3mult(X24->sig, rho))
		  - sigsig * rhopi;
    S12pipi = 3*cvec3mult(X13->sig, pi)*cvec3mult(X24->sig, pi) - sigsig*pi2;
  }

  if (P->tll || P->tpp || P->l2tpp) {
    S12ll = 3*cvec3mult(X13->sig, rhoxpi)*cvec3mult(X24->sig, rhoxpi) - sigsig*L2;
  }

  InteractionType type;
  double gamma, kappa;
  int ilabel;
  complex double iakapi, kapakapi, thelakapi;
  complex double Gi, PiPi, PiPiG, PirPir;
  complex double LL, L2LS, TLL;
  complex double TPP, TRP, L2TPP;
  
  complex double l2tpppipi, l2tpprhopi, l2tpprr, l2tppll;
  
  int pipidone=0, pipigdone=0, pirpirdone=0;
  int l2done=0, l2lsdone=0;
  int tlldone=0, tppdone=0;
  int trpdone=0, l2tppdone=0;
  int idx;

  for (idx=0; idx<P->ncomponents; idx++) {

    type = P->c[idx].type;
    gamma = P->c[idx].gamma;
    kappa = P->c[idx].kappa;
    ilabel = P->c[idx].ilabel;
    
    if ((idx==0 || kappa != P->c[idx-1].kappa) && type != VC) {
      iakapi = 1/(alpha+kappa);
      kapakapi = kappa*iakapi;
      if (P->l2 || P->l2ls || P->tll || P->tpp || P->prtrp || P->l2tpp)
        thelakapi = theta*iakapi + lambda*kapakapi;
	
      Gi = cpow32(kapakapi)*cexp(-0.5*rho2*iakapi);
      pipidone = 0;
      pipigdone = 0;
      pirpirdone = 0;
      l2done = 0;
      l2lsdone = 0;
      tlldone = 0;
      tppdone = 0;
      l2tppdone = 0;
      trpdone = 0;
    }


    // Calculation of non-overlap and non-Gaussian parts of the matrix elements
 
    // p2V
    if (!pipidone && p2V <= type && type <= tsp2V) {
      PiPi = pi2-0.5*beta*iakapi*rhopi+
	0.25*theta*csqr(iakapi)*rho2+
	0.75*(lambda-theta*iakapi);

      pipidone = 1;
    }

    // Vp2
    if (!pipigdone && Vp2 <= type && type <= tsVp2) {
      PiPiG = pi2-0.5*beta*iakapi*rhopi+
	(0.25*theta-0.5)*csqr(iakapi)*rho2+
	0.75*lambda-(0.75*theta-1.5)*iakapi;

      pipigdone = 1;
    }
        
    // pr2V
    if (!pirpirdone && pr2V <= type && type <= tspr2V) {
      PirPir = csqr(kapakapi)*(csqr(rhopi)-0.5*beta*iakapi*rhopi*rho2+
			       0.25*theta*csqr(iakapi)*csqr(rho2))+
	       alpha*kapakapi*(pi2-0.5*beta*iakapi*rhopi+
                               0.25*theta*csqr(iakapi)*rho2) -
	       2*csqr(kapakapi)*(-beta*rhopi+theta*iakapi*rho2) +
	       0.25*(csqr(kapakapi)*rho2+3*alpha*kapakapi)*(lambda-theta*iakapi) +
	       3*csqr(kapakapi)*theta;

      pirpirdone = 1;
    }

    // Vl2
    if (!l2done && Vl2 <= type && type <= tsVl2) {
          LL = kapakapi*(kapakapi*L2 + 2*alpha*pi2 - beta*rhopi + 0.5*thelakapi*rho2
    	     - 1.5*psi);
      
      l2done = 1;
    }


    // Vl2ls
    if (!l2lsdone && (type == Vl2ls || type == tVl2ls)) {
      L2LS = LS*csqr(kapakapi)*(kapakapi*L2 + 4*alpha*pi2 - 2*beta*rhopi + thelakapi*rho2
      		    - 5*psi) + 2*LS*kapakapi;
      
      l2lsdone = 1;
    } 

    // VTll - S12(l,l)
    if (!tlldone && (type == VTll || type == tVTll)) {
      TLL = csqr(kapakapi)*S12ll - alpha*kapakapi*S12pipi
      	    - 0.25*kapakapi*thelakapi*S12
	    + 0.5*kapakapi*beta*S12rhopi;

      tlldone = 1;
    }
    
    // VTpp - S12(p,p)
    if (!tppdone && (type == VTpp || type == tVTpp)) {
      TPP = csqr(kapakapi)*(kapakapi*S12ll*(5*alpha + kapakapi*rho2)
      			    + S12pipi*(9*csqr(alpha) + 13*alpha*kapakapi*rho2 + 2*csqr(kapakapi*rho2))
			    - S12rhopi*(4.5*alpha*beta + 16*alpha*kapakapi*rhopi + 2.5*beta*kapakapi*rho2
			    		+ 4*csqr(kapakapi)*rhopi*rho2)
			    + S12*(5.25*kapakapi*psi + 2.25*theta - 4.5
			    	   + 2*csqr(kapakapi*rhopi) + 4*kapakapi*beta*rhopi
				   - 0.75*kapakapi*thelakapi*rho2));
      tppdone = 1;
    }
    
    // NOTE: 06/08/03 Added factor 0.5 which was implied by hermitizing the super operator basis ...
    // prVTrp (p_r v(r) +v(r) p_r) S12(r,p)
    if (!trpdone && (type == prVTrp || type == tprVTrp)) {
      TRP = 0.5*csqr(kapakapi)*(S12pipi*2*alpha*(3*alpha + kapakapi*rho2)
      			    + S12*(1.5*theta - 2.625*kapakapi*csqr(beta) - 3
			    	   + 0.5*csqr(kapakapi)*beta*iakapi*rho2*rhopi
				   - 2*kapakapi*(alpha*pi2 + kapakapi*csqr(rhopi))
				   - 0.5*kapakapi*beta*(1 + 9*kapakapi)*rhopi
				   + 0.375*kapakapi*iakapi*csqr(beta)*rho2)
		            + S12rhopi*(-1.5*alpha*beta*(2 - 7*kapakapi)
			    	   - 0.5*csqr(kapakapi)*beta*iakapi*csqr(rho2)
				   + 2*kapakapi*(3*alpha + kapakapi*rho2)*rhopi
				   - 0.5*kapakapi*beta*(5 - 12*kapakapi)*rho2));
      trpdone = 1;
    }

    // Vl2Tpp - {L2 S12(p,p)}_H
    if (!l2tppdone && (type == Vl2Tpp || type == tVl2Tpp)) {
      l2tpppipi = 60*cpow(alpha, 3)*kapakapi*pi2 + 117*csqr(alpha*kapakapi)*pi2*rho2
      		   - 33*csqr(alpha*kapakapi*rhopi) 
      		   - 21*alpha*cpow(kapakapi,3)*rho2*csqr(rhopi) 
		   + 33*alpha*cpow(kapakapi,3)*pi2*csqr(rho2)
		   + 2*cpow(kapakapi,4)*L2*csqr(rho2)
		   + 3*cpow(kapakapi,3)*cpow(rho2,3)*thelakapi
		   - 6*beta*cpow(kapakapi,3)*csqr(rho2)*rhopi
		   + 3*csqr(alpha)*(9 - 35*kapakapi*psi)
		   - 42*alpha*beta*csqr(kapakapi)*rho2*rhopi
		   - 30*csqr(alpha)*beta*kapakapi*rhopi
		   + alpha*kapakapi*(74 - 196.5*kapakapi*psi + 15*theta)*rho2
		   + csqr(kapakapi)*(16 + 21*theta - 52.5*kapakapi*psi)*csqr(rho2);
		   
      l2tpprhopi = 24*alpha*cpow(kapakapi*rhopi,3)
		   - 96*csqr(alpha*kapakapi)*pi2*rhopi
		   - 48*alpha*cpow(kapakapi,3)*rhopi*rho2*pi2
		   - 4*cpow(kapakapi,4)*rho2*rhopi*L2
		   - 6*cpow(kapakapi,3)*thelakapi*csqr(rho2)*rhopi
		   - 4.5*beta*cpow(kapakapi,3)*rho2*L2
		   + 12*beta*cpow(kapakapi,3)*rho2*csqr(rhopi)
		   - 34.5*alpha*beta*csqr(kapakapi)*L2
		   + 30*alpha*beta*csqr(kapakapi*rhopi)
		   - 30*csqr(alpha)*beta*kapakapi*pi2
		   - 1.5*alpha*beta*(9 - 35*kapakapi*psi)
		   + alpha*kapakapi*(-138 + 60*theta + 201*kapakapi*psi)*rhopi
		   + csqr(kapakapi)*(-68 + 87*kapakapi*psi + 12*theta)*rho2*rhopi
		   - 4.5*beta*csqr(kapakapi)*thelakapi*csqr(rho2)
		   + beta*kapakapi*(-17.5 - 7.5*theta + 48*kapakapi*psi)*rho2;
		   
      l2tppll = 42*csqr(alpha*kapakapi)*pi2 - 9*alpha*cpow(kapakapi,3)*csqr(rhopi)
      		   + 15*alpha*cpow(kapakapi,3)*pi2*rho2 
		   + cpow(kapakapi,4)*rho2*L2
		   + 1.5*cpow(kapakapi,3)*thelakapi*csqr(rho2)
		   - 3*beta*cpow(kapakapi,3)*rhopi*rho2
		   - 21*alpha*beta*csqr(kapakapi)*rhopi
		   + alpha*kapakapi*(30 - 73.5*kapakapi*psi)
		   + csqr(kapakapi)*(6 - 24*kapakapi*psi + 10.5*theta)*rho2;
		   
      l2tpprr = 12*alpha*cpow(kapakapi,3)*pi2*csqr(rhopi)
      		   + 2*cpow(kapakapi,4)*csqr(rhopi)*L2
		   + 3.75*cpow(kapakapi,3)*thelakapi*rho2*csqr(rhopi)
		   - 0.75*cpow(kapakapi,3)*thelakapi*pi2*csqr(rho2)
		   - 0.75*csqr(kapakapi*thelakapi*rho2)
		   + 6*beta*cpow(kapakapi,3)*rhopi*(L2 - csqr(rhopi))
		   + 24*alpha*beta*csqr(kapakapi)*pi2*rhopi
		   + alpha*kapakapi*(-35 + 34.5*kapakapi*psi + 15*theta)*pi2
		   + csqr(kapakapi)*(78.5 - 41.25*kapakapi*psi - 56.25*theta)*csqr(rhopi)
		   + csqr(kapakapi)*(-14.5 + 12.75*kapakapi*psi + 5.25*theta)*pi2*rho2
		   + kapakapi*thelakapi*(-9.75 + 3.75*theta + 13.5*kapakapi*psi)*rho2
		   + 7.5*beta*csqr(kapakapi)*thelakapi*rho2*rhopi
		   + beta*kapakapi*(37 - 67.5*kapakapi*psi - 7.5*theta)*rhopi
		   - 13.5 + 6.75*theta - 47.25*csqr(kapakapi*psi)
		   + 5.25*kapakapi*(13 - 5*theta)*psi;
		   
      L2TPP = csqr(kapakapi)*(S12pipi*l2tpppipi + S12rhopi*l2tpprhopi + S12*l2tpprr + S12ll*l2tppll);
      
      l2tppdone = 1;
		    
    }
        
    // Central Potentials
    if (type == V && TT)
      v[ilabel] += gamma*TT*SS*RR*Gi;
    else if (type == sV && TT)
      v[ilabel] += gamma*TT*sigsig*RR*Gi;
    else if (type == tV && tautau)
      v[ilabel] += gamma*tautau*SS*RR*Gi;
    else if (type == tsV && tautau)
      v[ilabel] += gamma*tautau*sigsig*RR*Gi;

    // p2V Potentials
    else if (type == p2V && TT)
      v[ilabel] += gamma*TT*SS*RR*PiPi*Gi;
    else if (type == sp2V && TT)
      v[ilabel] += gamma*TT*sigsig*RR*PiPi*Gi;
    else if (type == tp2V && tautau)
      v[ilabel] += gamma*tautau*SS*RR*PiPi*Gi;
    else if (type == tsp2V && tautau)
      v[ilabel] += gamma*tautau*sigsig*RR*PiPi*Gi;

    // Vp2 Potentials
    else if (type == Vp2 && TT)
      v[ilabel] += gamma*TT*SS*RR*PiPiG*Gi;
    else if (type == sVp2 && TT)
      v[ilabel] += gamma*TT*sigsig*RR*PiPiG*Gi;
    else if (type == tVp2 && tautau)
      v[ilabel] += gamma*tautau*SS*RR*PiPiG*Gi;
    else if (type == tsVp2 && tautau)
      v[ilabel] += gamma*tautau*sigsig*RR*PiPiG*Gi;

    // pr2 Potentials
    else if (type == pr2V && TT)
      v[ilabel] += gamma*TT*SS*RR*PirPir*Gi;
    else if (type == spr2V && TT)
      v[ilabel] += gamma*TT*sigsig*RR*PirPir*Gi;
    else if (type == tpr2V && tautau)
      v[ilabel] += gamma*tautau*SS*RR*PirPir*Gi;
    else if (type == tspr2V && tautau)
      v[ilabel] += gamma*tautau*sigsig*RR*PirPir*Gi;

    // l2 Potentials
    else if (type == Vl2 && TT)
      v[ilabel] += gamma*TT*SS*RR*Gi*LL;
    else if (type == sVl2 && TT)
      v[ilabel] += gamma*TT*sigsig*RR*Gi*LL;
    else if (type == tVl2 && tautau)
      v[ilabel] += gamma*tautau*SS*RR*Gi*LL;
    else if (type == tsVl2 && tautau)
      v[ilabel] += gamma*tautau*sigsig*RR*Gi*LL;

    // Spin-Orbit Potentials
    else if (type == Vls && TT)
      v[ilabel] += gamma*TT*LS*kapakapi*RR*Gi;
    else if (type == tVls && tautau)
      v[ilabel] += gamma*tautau*LS*kapakapi*RR*Gi;

    // l2ls Potentials
    else if (type == Vl2ls && TT)
      v[ilabel] += gamma*TT*RR*Gi*L2LS;
    else if (type == tVl2ls && tautau)
      v[ilabel] += gamma*tautau*RR*Gi*L2LS;

    // Tensor Potentials
    else if (type == VT && TT)
      v[ilabel] += gamma*TT*S12*csqr(kapakapi)*RR*Gi;
    else if (type == tVT && tautau)
      v[ilabel] += gamma*tautau*S12*csqr(kapakapi)*RR*Gi;

    // S12(l,l) Potentials
    else if (type == VTll && TT)
      v[ilabel] += gamma*TT*RR*Gi*TLL;
    else if (type == tVTll && tautau)
      v[ilabel] += gamma*tautau*RR*Gi*TLL;

    // S12(p,p) Potentials
    else if (type == VTpp && TT)
      v[ilabel] += gamma*TT*RR*Gi*TPP;
    else if (type == tVTpp && tautau)
      v[ilabel] += gamma*tautau*RR*Gi*TPP;

    // {L2 S12(p,p)}_H Potentials
    else if (type == Vl2Tpp && TT)
      v[ilabel] += gamma*TT*RR*Gi*L2TPP;
    else if (type == tVl2Tpp && tautau)
      v[ilabel] += gamma*tautau*RR*Gi*L2TPP;
    
    // (p_r v(r) + v(r) p_r) S12(r,p) 
    else if (type == prVTrp && TT)
      v[ilabel] += gamma*TT*RR*Gi*TRP;
    else if (type == tprVTrp && tautau)
      v[ilabel] += gamma*tautau*RR*Gi*TRP;
    
    // Coulomb Potential
    else if (type == VC && TT && G1->xi == 1 && G2->xi == 1)
      v[ilabel] += gamma*TT*SS*RR/csqrt(2*alpha)*zcoulomb(0.5*rho2/alpha);
  }
}


void calcPotential(const Interaction *P,
		   const SlaterDet* Q, const SlaterDetAux* X, double v[])
{
  int i;
  TwoBodyOperator op_tb_pot = {dim: P->n, opt: 0, par: P, me: tb_pot};

  calcSlaterDetTBME(Q, X, &op_tb_pot, v);
  
  v[0] = 0.0;
  for (i=1; i<P->n; i++)
    v[0] += v[i];
}


void calcPotentialrowcol(const Interaction* P,
			 const SlaterDet* Q, const SlaterDetAux* X, double v[], 
			 int k, int l)
{
  int i;
  TwoBodyOperator op_tb_pot = {dim: P->n, opt: 0, par: P, me: tb_pot};

  calcSlaterDetTBMErowcol(Q, X, &op_tb_pot, v, k, l);
  
  v[0] = 0.0;
  for (i=1; i<P->n; i++)
    v[0] += v[i];
} 


void calcPotentialod(const Interaction *P,
		     const SlaterDet* Q, const SlaterDet* Qp,
		     const SlaterDetAux* X, 
		     complex double v[])
{
  int i;
  TwoBodyOperator op_tb_pot = {dim: P->n, opt: 0, par: P, me: tb_pot};

  calcSlaterDetTBMEod(Q, Qp, X, &op_tb_pot, v);

  v[0] = 0.0;
  for (i=1; i<P->n; i++)
    v[0] += v[i];
}


void calcPotentialodrowcol(const Interaction* P,
			   const SlaterDet* Q, const SlaterDet* Qp,
			   const SlaterDetAux* X, complex double v[], 
			   int k, int l)
{
  int i;
  TwoBodyOperator op_tb_pot = {dim: P->n, opt: 0, par: P, me: tb_pot};

  calcSlaterDetTBMEodrowcol(Q, Qp, X, &op_tb_pot, v, k, l);
  
  v[0] = 0.0;
  for (i=1; i<P->n; i++)
    v[0] += v[i];
} 


void calcPotentialsHF(const Interaction* Int,
		      const SlaterDet* Q, const SlaterDetAux* X,
		      void* mes)
{
  TwoBodyOperator op_tb_pot = {dim: Int->n, opt: 0, par: Int, me: tb_pot};
  calcSlaterDetTBHFMEs(Q, X, &op_tb_pot, mes);

  int A=Q->A;
  complex double (*v)[Int->n] = mes;
  int k,m,i;

  for (m=0; m<A; m++)
    for (k=0; k<A; k++) {
      v[k+m*A][0] = 0.0;
      for (i=1; i<Int->n; i++)
       v[k+m*A][0] += v[k+m*A][i];
    }
}      


void calcPotentialHF(const Interaction* Int,
		     const SlaterDet* Q, const SlaterDetAux* X,
		     complex double* vhf)
{
  int A=Q->A;
  complex double v[A*A][Int->n];

  calcPotentialsHF(Int, Q, X, v);

  int k,m;
  for (m=0; m<A; m++)
    for (k=0; k<A; k++)
      vhf[k+m*A] = v[k+m*A][0];
}

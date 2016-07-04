/**

  \file gradPotenial.c

  calculate gradients of interaction matrix elements


  (c) 2003 Thomas Neff

  Changes (HH):
  06/08/03
  + Added factor 0.5 to the gradients of prVTrp implicit in the prescription
    used for the super operator basis (hermitization).
  
  04/12/28
  + Added gradient for {L^2 S12(p,p)}_H
  
  04/09/03
  + Added gradients for S12(p,p) and (p_r v(r) + v(r)p_r)S12(r,p)
  
  04/07/21
  + Added gradients for L2, L2LS, and S12(L,L) interaction terms.
  + Slightly reorganized calculation of auxiliaries.


*/

#include <complex.h>

#include "Gaussian.h"
#include "gradGaussian.h"
#include "SlaterDet.h"
#include "gradSlaterDet.h"
#include "Interaction.h"

#include "gradPotential.h"

#include "numerics/cmath.h"
#include "numerics/coulomb.h"


static void gtb_pot(Interaction* P,
		    const Gaussian* G1, const Gaussian* G2, 
		    const Gaussian* G3, const Gaussian* G4, 
		    const GaussianAux* X13, const GaussianAux* X24, 
		    const gradGaussianAux* dX13,
		    complex double* v, gradGaussian* dv)
{	
  int TT, tautau;

  TT = X13->T * X24->T;
  tautau = ((1-G1->xi*G3->xi)*(1-G2->xi*G4->xi))/4+(G1->xi*G4->xi+G2->xi*G3->xi)/2; 

  if (!(TT || tautau)) return;

  complex double SS; 
  gradSpinor dSS;
  complex double RR;
  gradScalar dRR;
  complex double lambda, alpha, rho[3], rho2, pi[3], pi2;
  complex double dlambda, dalpha;
  gradVector drho, dpi;
  gradScalar drho2, dpi2;
  complex double sigsig;
  gradSpinor dsigsig;

  complex double ltheta, rtheta;
  complex double dltheta, drtheta;
  complex double beta, dbeta, theta, dtheta;
  complex double psi, dpsi;
  
  complex double rhopi;
  gradScalar drhopi;

  complex double rhoxpi[3];
  complex double L2;
  gradScalar dL2;
  
  complex double S[3], LS;
  gradGaussian dLS;
  
  complex double lsigrho, rsigrho, S12;
  gradGaussian dS12;

  complex double lsigpi, rsigpi, lsigrhoxpi, rsigrhoxpi;
  
  complex double S12ll, S12pipi, S12rhopi;
  gradGaussian dS12ll, dS12pipi, dS12rhopi;

  complex double vcoul;
  gradScalar dvcoul;

  int i;
  
  SS = X13->S * X24->S; 
  for (i=0; i<2; i++)
    dSS.chi[i] = dX13->dS.chi[i]* X24->S;
  
  RR = X13->R * X24->R;
  dRR.a = dX13->dR.a* X24->R;
  for (i=0; i<3; i++)
    dRR.b[i] = dX13->dR.b[i]* X24->R;
  
  lambda = X13->lambda + X24->lambda;
  dlambda = dX13->dlambda;

  alpha = X13->alpha + X24->alpha;
  dalpha = dX13->dalpha;
  
  for (i=0; i<3; i++) { 
    rho[i] = X13->rho[i] - X24->rho[i];
    drho.a[i] = dX13->drho.a[i];
  }
  drho.b = dX13->drho.b;

  rho2 = cvec3sqr(rho);
  drho2.a = 2*cvec3mult(drho.a, rho);
  for (i=0; i<3; i++)
    drho2.b[i] = 2*drho.b*rho[i];
  
  for (i=0; i<3; i++) {
    pi[i] = 0.5*(X13->pi[i] - X24->pi[i]);
    dpi.a[i] = 0.5*dX13->dpi.a[i];
  }
  dpi.b = 0.5*dX13->dpi.b;

  sigsig = cvec3mult(X13->sig, X24->sig);
  for (i=0; i<2; i++)
    dsigsig.chi[i] = cvec3mult(dX13->dsig.chi[i], X24->sig);

  if (P->momentump2 || P->momentumpr2 || P->l2 || P->l2ls || P->tll || P->tpp || P->prtrp || P->l2tpp) {
    rtheta = (conj(G1->a)*X13->lambda + conj(G2->a)*X24->lambda);
    ltheta = (G3->a*X13->lambda + G4->a*X24->lambda);
    dltheta = - G3->a* csqr(X13->lambda);
    drtheta = X13->lambda - conj(G1->a)*csqr(X13->lambda);
    
    beta = I*(rtheta-ltheta);
    dbeta = I*(drtheta-dltheta);
    
    theta = ltheta*rtheta;
    dtheta = dltheta*rtheta+ltheta*drtheta;
 
    psi = theta - alpha*lambda;
    dpsi = dtheta - dalpha*lambda - alpha*dlambda;
    
    pi2 = cvec3sqr(pi);
    dpi2.a = 2*cvec3mult(dpi.a, pi);
    for (i=0; i<3; i++)
      dpi2.b[i] = 2*dpi.b*pi[i];    

    rhopi = cvec3mult(rho, pi);
    drhopi.a = cvec3mult(drho.a, pi)+cvec3mult(rho, dpi.a);
    for (i=0; i<3; i++)
      drhopi.b[i] = drho.b*pi[i]+rho[i]*dpi.b;
  }
    
  if (P->l2 || P->spinorbit || P->l2ls || P->tll || P->tpp || P->l2tpp){
    cvec3cross(rho, pi, rhoxpi);
    
    L2 = cvec3sqr(rhoxpi);
    
    dL2.a = drho2.a*pi2 + rho2*dpi2.a - 2*rhopi*drhopi.a;
    for(i=0; i<3; i++)
      dL2.b[i] = drho2.b[i]*pi2 + rho2*dpi2.b[i] - 2*rhopi*drhopi.b[i];
  }
    
  if (P->spinorbit || P->l2ls) {
    for (i=0; i<3; i++)
      S[i] = 0.5*(X13->sig[i]*X24->S+X13->S*X24->sig[i]);
    LS = cvec3mult(rhoxpi, S);

    for (i=0; i<2; i++)
      dLS.chi[i] = 0.5*(cvec3mult(dX13->dsig.chi[i], rhoxpi)*X24->S +
			dX13->dS.chi[i]*cvec3mult(X24->sig, rhoxpi));
    
    dLS.a = cvec3spat(drho.a, pi, S) + cvec3spat(rho, dpi.a, S);
    dLS.b[0] = (drho.b*pi[1]-rho[1]*dpi.b)*S[2] - 
               (drho.b*pi[2]-rho[2]*dpi.b)*S[1];
    dLS.b[1] = (drho.b*pi[2]-rho[2]*dpi.b)*S[0] - 
               (drho.b*pi[0]-rho[0]*dpi.b)*S[2];
    dLS.b[2] = (drho.b*pi[0]-rho[0]*dpi.b)*S[1] - 
               (drho.b*pi[1]-rho[1]*dpi.b)*S[0];
  }

  if (P->tensor || P->tll || P->tpp || P->prtrp || P->l2tpp) {
    lsigrho = cvec3mult(X13->sig, rho); 
    rsigrho = cvec3mult(X24->sig, rho);
    
    S12 = 3*lsigrho*rsigrho - sigsig*rho2;

    for (i=0; i<2; i++)
      dS12.chi[i] = 3*cvec3mult(dX13->dsig.chi[i], rho)*rsigrho - dsigsig.chi[i]*rho2;
    
    dS12.a = 3*(cvec3mult(X13->sig, drho.a)*rsigrho + lsigrho*cvec3mult(X24->sig, drho.a)) - 
      sigsig*drho2.a;
	      
    for (i=0; i<3; i++)
      dS12.b[i] = 3*drho.b*(X13->sig[i]*rsigrho + lsigrho*X24->sig[i]) -
	sigsig* drho2.b[i];
  }
    
  if (P->tll || P->tpp || P->prtrp || P->l2tpp) { 
    lsigpi = cvec3mult(X13->sig, pi);
    rsigpi = cvec3mult(X24->sig, pi);
    		
    // S12pipi
    S12pipi = 3*lsigpi*rsigpi - sigsig*pi2;
    
    dS12pipi.a = 3*(cvec3mult(X13->sig, dpi.a)*rsigpi + lsigpi*cvec3mult(X24->sig, dpi.a)) -
                 sigsig*dpi2.a;
 
    for (i=0; i<3; i++)
      dS12pipi.b[i] = 3*dpi.b*(X13->sig[i]*rsigpi + lsigpi*X24->sig[i]) - sigsig*dpi2.b[i];

    for (i=0; i<2; i++)
      dS12pipi.chi[i] = 3*cvec3mult(dX13->dsig.chi[i], pi)*rsigpi - dsigsig.chi[i]*pi2;
		 
    // S12rhopi    
    S12rhopi = 1.5*(lsigrho*rsigpi + lsigpi*rsigrho) - sigsig*rhopi;
	
    dS12rhopi.a = 1.5*(cvec3mult(X13->sig,drho.a)*rsigpi + 
    		       lsigrho*cvec3mult(X24->sig,dpi.a) +
    		       cvec3mult(X13->sig,dpi.a)*rsigrho + 
		       lsigpi*cvec3mult(X24->sig,drho.a)) -
		  sigsig*drhopi.a;

    for (i=0; i<3; i++)
      dS12rhopi.b[i] = 1.5*(drho.b*X13->sig[i]*rsigpi + lsigrho*dpi.b*X24->sig[i] +
       			    dpi.b*X13->sig[i]*rsigrho + lsigpi*drho.b*X24->sig[i]) -
		       sigsig*drhopi.b[i];
      
    for (i=0; i<2; i++)
      dS12rhopi.chi[i] = 1.5*(cvec3mult(dX13->dsig.chi[i],rho)*rsigpi +
            	              cvec3mult(dX13->dsig.chi[i],pi)*rsigrho) -
			 dsigsig.chi[i]*rhopi;  
    
  }
    
  if (P->tll || P->tpp || P->l2tpp) {
    lsigrhoxpi = cvec3mult(X13->sig, rhoxpi);
    rsigrhoxpi = cvec3mult(X24->sig, rhoxpi);
  
    // S12ll
    S12ll = 3*lsigrhoxpi*rsigrhoxpi - sigsig*L2;
    
    dS12ll.a = 3*((cvec3spat(X13->sig, drho.a, pi) + 
    		   cvec3spat(X13->sig, rho, dpi.a))*rsigrhoxpi +
                  (cvec3spat(X24->sig, drho.a, pi) + 
		   cvec3spat(X24->sig, rho, dpi.a))*lsigrhoxpi) -
	       sigsig*dL2.a;
	       
    dS12ll.b[0] = 3*((X13->sig[1]*(dpi.b*rho[2] - drho.b*pi[2]) -
                      X13->sig[2]*(dpi.b*rho[1] - drho.b*pi[1]))*rsigrhoxpi +
		     (X24->sig[1]*(dpi.b*rho[2] - drho.b*pi[2]) -
		      X24->sig[2]*(dpi.b*rho[1] - drho.b*pi[1]))*lsigrhoxpi) -
		  sigsig*dL2.b[0];

    dS12ll.b[1] = 3*((X13->sig[2]*(dpi.b*rho[0] - drho.b*pi[0]) -
                      X13->sig[0]*(dpi.b*rho[2] - drho.b*pi[2]))*rsigrhoxpi +
		     (X24->sig[2]*(dpi.b*rho[0] - drho.b*pi[0]) -
		      X24->sig[0]*(dpi.b*rho[2] - drho.b*pi[2]))*lsigrhoxpi) -
		  sigsig*dL2.b[1];
    
    dS12ll.b[2] = 3*((X13->sig[0]*(dpi.b*rho[1] - drho.b*pi[1]) -
                      X13->sig[1]*(dpi.b*rho[0] - drho.b*pi[0]))*rsigrhoxpi +
		     (X24->sig[0]*(dpi.b*rho[1] - drho.b*pi[1]) -
		      X24->sig[1]*(dpi.b*rho[0] - drho.b*pi[0]))*lsigrhoxpi) -
		  sigsig*dL2.b[2];
		      		
    for (i=0; i<2; i++)
      dS12ll.chi[i] = 3*cvec3mult(dX13->dsig.chi[i], rhoxpi)*rsigrhoxpi - 
      		      dsigsig.chi[i]*L2;

  }

  // kappa-dependent auxiliaries
  InteractionType type;
  double gamma, kappa;
  complex double iakapi, kapakapi, diakapi, dkapakapi;
  complex double thelakapi, dthelakapi;
  complex double Gi;
  gradScalar dGi;

  
  complex double PiPi, PiPiG, PirPir;
  gradScalar dPiPi, dPiPiG, dPirPir;
  int pipidone=0, pipigdone=0, pirpirdone=0;
  
  complex double LL;
  gradScalar dLL;
  int l2done=0;
  
  complex double L2LS;
  gradGaussian dL2LS;
  int l2lsdone=0;
  
  complex double TLL, TPP, TRP;
  gradGaussian dTLL, dTPP, dTRP;
  int tlldone=0, tppdone=0, trpdone=0;
  
  complex double L2TPP;
  gradGaussian dL2TPP;
  int l2tppdone=0;
  
  complex double tppll, tpppipi, tpprhopi, tpprr;
  complex double trppipi, trprhopi, trprr;
  complex double l2tpppipi, l2tpprhopi, l2tpprr, l2tppll;
  
  gradScalar dl2tpppipi, dl2tpprhopi, dl2tpprr, dl2tppll;
  
  int idx;

  for (idx=0; idx<P->ncomponents; idx++) {

    type = P->c[idx].type;
    gamma = P->c[idx].gamma;
    kappa = P->c[idx].kappa;
    
    if ((idx==0 || kappa != P->c[idx-1].kappa) && type != VC) {
      iakapi = 1.0/(alpha+kappa);
      diakapi = -dalpha*csqr(iakapi);
      kapakapi = kappa*iakapi;
      dkapakapi = kappa*diakapi;
      if (P->l2 || P->l2ls || P->tll || P->tpp || P->prtrp || P->l2tpp) {
        thelakapi = theta*iakapi + lambda*kapakapi;
	dthelakapi = dtheta*iakapi + theta*diakapi + dlambda*kapakapi + lambda*dkapakapi;
      }
      
      Gi = cpow32(kapakapi)*cexp(-0.5*rho2*iakapi);
      dGi.a = (1.5*dkapakapi/kapakapi-0.5*diakapi*rho2-0.5*iakapi*drho2.a)*Gi;
      for (i=0; i<3; i++)
	dGi.b[i] = -0.5*drho2.b[i]*iakapi*Gi;
      
      pipidone = 0;
      pipigdone = 0;
      pirpirdone = 0;
      l2done = 0;
      l2lsdone = 0;
      tlldone = 0;
      tppdone = 0;
      trpdone = 0;
      l2tppdone = 0;
    }

    // p2V
    if (!pipidone && (p2V <= type && type <= tsp2V)) {
      PiPi = pi2-0.5*beta*iakapi*rhopi+
	0.25*theta*csqr(iakapi)*rho2+
	0.75*(lambda-theta*iakapi);
      
      dPiPi.a = dpi2.a - 0.5*dbeta*iakapi*rhopi - 0.5*beta*diakapi*rhopi -
	0.5*beta*iakapi*drhopi.a + 0.25*dtheta*csqr(iakapi)*rho2 +
	0.5*theta*diakapi*iakapi*rho2 + 0.25*theta*csqr(iakapi)*drho2.a +
	0.75*dlambda - 0.75*dtheta*iakapi - 0.75*theta*diakapi;
      for (i=0; i<3; i++)
	dPiPi.b[i] = dpi2.b[i] - 0.5*beta*iakapi*drhopi.b[i] + 
	             0.25*theta*csqr(iakapi)*drho2.b[i];

      pipidone = 1;
    }

    // Vp2
    if (!pipigdone && (Vp2 <= type && type <= tsVp2)) {
      PiPiG = pi2-0.5*beta*iakapi*rhopi+
	(0.25*theta-0.5)*csqr(iakapi)*rho2+
	0.75*lambda-(0.75*theta-1.5)*iakapi;
      
      dPiPiG.a = dpi2.a - 0.5*dbeta*iakapi*rhopi - 0.5*beta*diakapi*rhopi -
	0.5*beta*iakapi*drhopi.a + 0.25*dtheta*csqr(iakapi)*rho2 +
	(0.5*theta-1.0)*diakapi*iakapi*rho2 + (0.25*theta-0.5)*csqr(iakapi)*drho2.a +
	0.75*dlambda - 0.75*dtheta*iakapi - (0.75*theta-1.5)*diakapi;
      for (i=0; i<3; i++)
	dPiPiG.b[i] = dpi2.b[i] - 0.5*beta*iakapi*drhopi.b[i] + 
          (0.25*theta-0.5)*csqr(iakapi)*drho2.b[i];

      pipigdone = 1;
    }

    // pr2V
    if (!pirpirdone && (pr2V <= type && type <= tspr2V)) {
      PirPir = 
	csqr(kapakapi)*(csqr(rhopi)-0.5*beta*iakapi*rhopi*rho2+
			0.25*theta*csqr(iakapi)*csqr(rho2))+
	alpha*kapakapi*(pi2-0.5*beta*iakapi*rhopi+0.25*theta*csqr(iakapi)*rho2) -
	2*csqr(kapakapi)*(-beta*rhopi+theta*iakapi*rho2) +
	0.25*(csqr(kapakapi)*rho2+3*alpha*kapakapi)*(lambda-theta*iakapi) +
	3*csqr(kapakapi)*theta;
      
      dPirPir.a = 
	2*kapakapi*dkapakapi*
	(csqr(rhopi)-0.5*beta*iakapi*rhopi*rho2+0.25*theta*csqr(iakapi)*csqr(rho2))+
	csqr(kapakapi)*
	(2*rhopi*drhopi.a-0.5*dbeta*iakapi*rhopi*rho2-0.5*beta*diakapi*rhopi*rho2-
	 0.5*beta*iakapi*drhopi.a*rho2-0.5*beta*iakapi*rhopi*drho2.a+
	 0.25*dtheta*csqr(iakapi)*csqr(rho2)+0.5*theta*iakapi*diakapi*csqr(rho2)+
	 0.5*theta*csqr(iakapi)*rho2*drho2.a) +
	(dalpha*kapakapi+alpha*dkapakapi)*
	(pi2-0.5*beta*iakapi*rhopi+0.25*theta*csqr(iakapi)*rho2)+
	alpha*kapakapi*
	(dpi2.a-0.5*dbeta*iakapi*rhopi-0.5*beta*diakapi*rhopi-
	 0.5*beta*iakapi*drhopi.a+0.25*dtheta*csqr(iakapi)*rho2+
	 0.5*theta*iakapi*diakapi*rho2+0.25*theta*csqr(iakapi)*drho2.a) -
	4*kapakapi*dkapakapi*(-beta*rhopi+theta*iakapi*rho2) -
	2*csqr(kapakapi)*
	(-dbeta*rhopi-beta*drhopi.a+dtheta*iakapi*rho2+theta*diakapi*rho2+
	 theta*iakapi*drho2.a) +
	0.25*(2*kapakapi*dkapakapi*rho2+csqr(kapakapi)*drho2.a+3*dalpha*kapakapi+
	      3*alpha*dkapakapi)*(lambda-theta*iakapi)+
	0.25*(csqr(kapakapi)*rho2+3*alpha*kapakapi)*(dlambda-dtheta*iakapi-theta*diakapi) +
	6*kapakapi*dkapakapi*theta + 3*csqr(kapakapi)*dtheta;

      for (i=0; i<3; i++)
	dPirPir.b[i] =  
	  csqr(kapakapi)*(2*rhopi*drhopi.b[i]-0.5*beta*iakapi*drhopi.b[i]*rho2-
			  0.5*beta*iakapi*rhopi*drho2.b[i]+
			  0.5*theta*csqr(iakapi)*rho2*drho2.b[i]) +
	  alpha*kapakapi*
	  (dpi2.b[i]-0.5*beta*iakapi*drhopi.b[i]+0.25*theta*csqr(iakapi)*drho2.b[i]) -
	  2*csqr(kapakapi)*(-beta*drhopi.b[i]+theta*iakapi*drho2.b[i]) +
	  0.25*(csqr(kapakapi)*drho2.b[i])*(lambda-theta*iakapi);

      pirpirdone = 1;
    }

    // Vl2
    if (!l2done && (Vl2 <= type && type <= tsVl2)) {
      LL = kapakapi*
             (kapakapi*L2 + 2*alpha*pi2 - beta*rhopi + 0.5*(theta*iakapi+lambda*kapakapi)*rho2 -
    	      1.5*(theta - alpha*lambda));
	     
      dLL.a = dkapakapi*LL/kapakapi +
	      kapakapi*
	        (dkapakapi*L2 + kapakapi*dL2.a + 2*dalpha*pi2 + 2*alpha*dpi2.a -
	         dbeta*rhopi - beta*drhopi.a +
		 0.5*(dtheta*iakapi + theta*diakapi + dlambda*kapakapi + lambda*dkapakapi)*rho2 +
		 0.5*(theta*iakapi + lambda*kapakapi)*drho2.a - 
		 1.5*(dtheta - dalpha*lambda - alpha*dlambda));
	      
      for(i=0; i<3; i++)
        dLL.b[i] = kapakapi*
		     (kapakapi*dL2.b[i] + 2*alpha*dpi2.b[i] - beta*drhopi.b[i] +
	              0.5*(theta*iakapi+lambda*kapakapi)*drho2.b[i]);  
	      
      l2done = 1;
	      
    }


    // Vl2ls
    if (!l2lsdone && (type == Vl2ls || type == tVl2ls)) {
      L2LS = LS*csqr(kapakapi)*
               (kapakapi*L2 + 4*alpha*pi2 - 2*beta*rhopi + (theta*iakapi + lambda*kapakapi)*rho2 -
      		5*(theta - alpha*lambda)) + 
	     2*LS*kapakapi;
		    
      dL2LS.a = (dLS.a*csqr(kapakapi)+2*kapakapi*dkapakapi*LS)*
      		  (kapakapi*L2 + 4*alpha*pi2 - 2*beta*rhopi + (theta*iakapi + lambda*kapakapi)*rho2 - 
		   5*(theta - alpha*lambda)) +
		LS*csqr(kapakapi)*
		  (dkapakapi*L2 + kapakapi*dL2.a + 
		   4*(dalpha*pi2 + alpha*dpi2.a) - 2*(dbeta*rhopi + beta*drhopi.a) +
		   (dtheta*iakapi + theta*diakapi + dlambda*kapakapi + lambda*dkapakapi)*rho2 +
		   (theta*iakapi + lambda*kapakapi)*drho2.a -
		   5*(dtheta - dalpha*lambda - alpha*dlambda)) +
		2*dLS.a*kapakapi + 2*LS*dkapakapi;
		
      for(i=0; i<3; i++)
        dL2LS.b[i] = csqr(kapakapi)*dLS.b[i]*
	               (kapakapi*L2 + 4*alpha*pi2 - 2*beta*rhopi + (theta*iakapi + lambda*kapakapi)*rho2 -
      		        5*(theta - alpha*lambda)) +
		     csqr(kapakapi)*LS*
		       (kapakapi*dL2.b[i] + 4*alpha*dpi2.b[i] - 2*beta*drhopi.b[i] +
		        (theta*iakapi + lambda*kapakapi)*drho2.b[i]) +
		     2*kapakapi*dLS.b[i];
		     
      for(i=0; i<2; i++)
        dL2LS.chi[i] = csqr(kapakapi)*dLS.chi[i]*
	                 (kapakapi*L2 + 4*alpha*pi2 - 2*beta*rhopi + (theta*iakapi + lambda*kapakapi)*rho2 - 
			  5*(theta - alpha*lambda)) + 2*dLS.chi[i]*kapakapi;

      l2lsdone = 1;		    
    }		    
    
    // VTll - S12(l,l)
    if (!tlldone && (type == VTll || type == tVTll)) {
      TLL = kapakapi*(kapakapi*S12ll - alpha*S12pipi -
      	    	      0.25*(theta*iakapi + lambda*kapakapi)*S12 +
	              0.5*beta*S12rhopi);
	    
      dTLL.a = dkapakapi*TLL/kapakapi +
      	       kapakapi*
	         (dkapakapi*S12ll + kapakapi*dS12ll.a - dalpha*S12pipi - alpha*dS12pipi.a -
	       	  0.25*(dtheta*iakapi + theta*diakapi + dlambda*kapakapi + lambda*dkapakapi)*S12 -
		  0.25*(theta*iakapi + lambda*kapakapi)*dS12.a +
		  0.5*dbeta*S12rhopi + 0.5*beta*dS12rhopi.a);
			 
      for(i=0; i<3; i++)
        dTLL.b[i] = kapakapi*(kapakapi*dS12ll.b[i] - alpha*dS12pipi.b[i] - 
	  	    	      0.25*(theta*iakapi + lambda*kapakapi)*dS12.b[i] +
		    	      0.5*beta*dS12rhopi.b[i]);
			      
      for(i=0; i<2; i++)
        dTLL.chi[i] = kapakapi*(kapakapi*dS12ll.chi[i] - alpha*dS12pipi.chi[i] - 
			        0.25*(theta*iakapi + lambda*kapakapi)*dS12.chi[i] +
				0.5*beta*dS12rhopi.chi[i]);	 
      tlldone = 1;
    }
    
    // VTpp - S12(p,p)
    if (!tppdone && (type == VTpp || type == tVTpp)) {
      tppll = 5*alpha + kapakapi*rho2;    
      
      tpppipi = 9*csqr(alpha) + 13*alpha*kapakapi*rho2 + 2*csqr(kapakapi*rho2);
      
      tpprhopi = 4.5*alpha*beta + 16*alpha*kapakapi*rhopi + 2.5*beta*kapakapi*rho2 + 
      		 4*csqr(kapakapi)*rhopi*rho2;
      
      tpprr = 5.25*kapakapi*(theta - alpha*lambda) + 2.25*theta - 4.5 + 
      	      2*csqr(kapakapi*rhopi) + 4*kapakapi*beta*rhopi - 
	      0.75*kapakapi*(theta*iakapi + kapakapi*lambda)*rho2;
      
	      
      TPP = csqr(kapakapi)*
      	      (kapakapi*S12ll*tppll + S12pipi*tpppipi - S12rhopi*tpprhopi + S12*tpprr);
			    
      dTPP.a = 2*TPP/kapakapi*dkapakapi + 
      	       csqr(kapakapi)*
	         ((dkapakapi*S12ll + kapakapi*dS12ll.a)*tppll +
	       	   kapakapi*S12ll*(5*dalpha + dkapakapi*rho2 + kapakapi*drho2.a) +
		   dS12pipi.a*tpppipi +
		   S12pipi*
		     (18*alpha*dalpha + 
		      13*(dalpha*kapakapi*rho2 + alpha*dkapakapi*rho2 + alpha*kapakapi*drho2.a) +
		      4*kapakapi*rho2*(dkapakapi*rho2 + kapakapi*drho2.a)) -
		   dS12rhopi.a*tpprhopi -
		   S12rhopi*
		     (4.5*(dalpha*beta + alpha*dbeta) +
		      16*(dalpha*kapakapi*rhopi + alpha*dkapakapi*rhopi + alpha*kapakapi*drhopi.a) +
		      2.5*(dbeta*kapakapi*rho2 + beta*dkapakapi*rho2 + beta*kapakapi*drho2.a) +
		      4*(2*kapakapi*dkapakapi*rhopi*rho2 + 
		      	 csqr(kapakapi)*(drhopi.a*rho2 + rhopi*drho2.a))) +
 		   dS12.a*tpprr +
	           S12*(5.25*dkapakapi*(theta - alpha*lambda) + 
		        5.25*kapakapi*(dtheta - dalpha*lambda - alpha*dlambda) +
		   	2.25*dtheta + 4*kapakapi*rhopi*(dkapakapi*rhopi + kapakapi*drhopi.a) +
			4*(dkapakapi*beta*rhopi + kapakapi*dbeta*rhopi + kapakapi*beta*drhopi.a) -
			0.75*(dkapakapi*rho2 + kapakapi*drho2.a)*(theta*iakapi + kapakapi*lambda) -
			0.75*kapakapi*
			  (dtheta*iakapi + theta*diakapi + dkapakapi*lambda + kapakapi*dlambda)*rho2));
			  
      for(i=0; i<3; i++)
        dTPP.b[i] = csqr(kapakapi)*
		    (kapakapi*dS12ll.b[i]*tppll + csqr(kapakapi)*S12ll*drho2.b[i] +
		     dS12pipi.b[i]*tpppipi +
		     S12pipi*(13*alpha*kapakapi*drho2.b[i] + 4*csqr(kapakapi)*rho2*drho2.b[i]) -
		     dS12rhopi.b[i]*tpprhopi -
		     S12rhopi*(16*alpha*kapakapi*drhopi.b[i] + 2.5*beta*kapakapi*drho2.b[i] +
		     	       4*csqr(kapakapi)*(drhopi.b[i]*rho2 + rhopi*drho2.b[i])) +
		     dS12.b[i]*tpprr +
		     S12*(4*csqr(kapakapi)*rhopi*drhopi.b[i] + 4*kapakapi*beta*drhopi.b[i] -
		     	  0.75*kapakapi*(theta*iakapi + kapakapi*lambda)*drho2.b[i]));
			    
      for(i=0; i<2; i++)
        dTPP.chi[i] = csqr(kapakapi)*
		      (kapakapi*dS12ll.chi[i]*tppll + dS12pipi.chi[i]*tpppipi - 
		       dS12rhopi.chi[i]*tpprhopi + dS12.chi[i]*tpprr); 
      						
      tppdone = 1;
    }
    
    // prVTrp (p_r v(r) +v(r) p_r) S12(r,p)
    if (!trpdone && (type == prVTrp || type == tprVTrp)) {
      trppipi = 2*alpha*(3*alpha + kapakapi*rho2);
      
      trprhopi = -1.5*alpha*beta*(2 - 7*kapakapi) - 0.5*csqr(kapakapi)*beta*iakapi*csqr(rho2) +
      	         2*kapakapi*(3*alpha + kapakapi*rho2)*rhopi - 
		 0.5*kapakapi*beta*(5 - 12*kapakapi)*rho2;
      
      trprr = 1.5*theta - 2.625*kapakapi*csqr(beta) - 3 + 
              0.5*csqr(kapakapi)*beta*iakapi*rho2*rhopi -
	      2*kapakapi*(alpha*pi2 + kapakapi*csqr(rhopi)) - 
	      0.5*kapakapi*beta*(1 + 9*kapakapi)*rhopi +
	      0.375*kapakapi*iakapi*csqr(beta)*rho2; 
      
      // NOTE: 06/08/03 Added factor 0.5 implied in the hermitization of the super operator basis
      TRP = 0.5*csqr(kapakapi)*(S12pipi*trppipi + S12*trprr + S12rhopi*trprhopi);
      
      // NOTE: 06/08/03 Added factor 0.5 in the (long) 2nd term of dTRP.a - it's included in the 
      //       1st term via TRP ! Corrected dTRP.b and dTRP.chi, too.
      dTRP.a = 2*TRP/kapakapi*dkapakapi +
               0.5*csqr(kapakapi)*
	         (dS12pipi.a*trppipi +
		  S12pipi*2*(6*alpha*dalpha + dalpha*kapakapi*rho2 + 
		  	     alpha*(dkapakapi*rho2 + kapakapi*drho2.a)) +
		  dS12rhopi.a*trprhopi +
		  S12rhopi*
		    (-1.5*(dalpha*beta + alpha*dbeta)*(2-7*kapakapi) + 
		     10.5*alpha*beta*dkapakapi - kapakapi*dkapakapi*beta*iakapi*csqr(rho2) -
		      0.5*csqr(kapakapi)*
		        (dbeta*iakapi*csqr(rho2) + beta*diakapi*csqr(rho2) + 
			 2*beta*iakapi*rho2*drho2.a) +
		      2*dkapakapi*(3*alpha + kapakapi*rho2)*rhopi +
		      2*kapakapi*(3*dalpha + dkapakapi*rho2 + kapakapi*drho2.a)*rhopi +
		      2*kapakapi*(3*alpha + kapakapi*rho2)*drhopi.a -
		      0.5*(dkapakapi*beta + kapakapi*dbeta)*(5 - 12*kapakapi)*rho2 +
		      6*kapakapi*beta*dkapakapi*rho2 -
	              0.5*kapakapi*beta*(5 - 12*kapakapi)*drho2.a) +
		  dS12.a*trprr +
		  S12*(1.5*dtheta - 2.625*dkapakapi*csqr(beta) - 5.25*kapakapi*beta*dbeta +
		       kapakapi*dkapakapi*beta*iakapi*rho2*rhopi +
		       0.5*csqr(kapakapi)*(dbeta*iakapi + beta*diakapi)*rho2*rhopi +
		       0.5*csqr(kapakapi)*beta*iakapi*(drho2.a*rhopi + rho2*drhopi.a) -
		       2*dkapakapi*(alpha*pi2 + kapakapi*csqr(rhopi)) -
		       2*kapakapi*(dalpha*pi2 + alpha*dpi2.a + dkapakapi*csqr(rhopi) + 
		       2*kapakapi*rhopi*drhopi.a) -
		       0.5*(dkapakapi*beta + kapakapi*dbeta)*(1 + 9*kapakapi)*rhopi -
		       4.5*kapakapi*beta*dkapakapi*rhopi -
		       0.5*kapakapi*beta*(1 + 9*kapakapi)*drhopi.a +
		       0.375*(dkapakapi*iakapi + kapakapi*diakapi)*csqr(beta)*rho2 +
		       0.75*kapakapi*iakapi*beta*dbeta*rho2 +
		       0.375*kapakapi*iakapi*csqr(beta)*drho2.a));
      
      for(i=0; i<3; i++)
        dTRP.b[i] = 0.5*csqr(kapakapi)*
	            (dS12pipi.b[i]*trppipi + S12pipi*2*alpha*kapakapi*drho2.b[i] +
		     dS12rhopi.b[i]*trprhopi +
		     S12rhopi*(-csqr(kapakapi)*beta*iakapi*rho2*drho2.b[i] +
		                6*kapakapi*alpha*drhopi.b[i] +
				2*csqr(kapakapi)*(drho2.b[i]*rhopi + rho2*drhopi.b[i]) -
				0.5*kapakapi*beta*(5 - 12*kapakapi)*drho2.b[i]) +
		     dS12.b[i]*trprr +
		     S12*(0.5*csqr(kapakapi)*beta*iakapi*(drho2.b[i]*rhopi +rho2*drhopi.b[i]) -
		          2*kapakapi*(alpha*dpi2.b[i] + 2*kapakapi*rhopi*drhopi.b[i]) -
			  0.5*kapakapi*beta*(1 + 9*kapakapi)*drhopi.b[i] +
			  0.375*kapakapi*iakapi*csqr(beta)*drho2.b[i]));
			 
      for(i=0; i<2; i++)
        dTRP.chi[i] = 0.5*csqr(kapakapi)*
			(dS12pipi.chi[i]*trppipi + dS12rhopi.chi[i]*trprhopi + dS12.chi[i]*trprr); 
	
      trpdone = 1;
    }
    
    // Vl2Tpp - {L2 S12(p,p)}_H
    if (!l2tppdone && (type == Vl2Tpp || type == tVl2Tpp)) {
      
      // l2tpppipi
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

      dl2tpppipi.a = 180*csqr(alpha)*dalpha*kapakapi*pi2 + 60*cpow(alpha,3)*(dkapakapi*pi2 + kapakapi*dpi2.a)
      		   + 234*(alpha*dalpha*csqr(kapakapi) + csqr(alpha)*dkapakapi*kapakapi)*pi2*rho2 
		     + 117*csqr(alpha*kapakapi)*(dpi2.a*rho2 + pi2*drho2.a)
		   - 66*(alpha*dalpha*csqr(kapakapi*rhopi) + kapakapi*dkapakapi*csqr(alpha*rhopi)
		   	 + csqr(alpha*kapakapi)*rhopi*drhopi.a)
		   - 21*(dalpha*cpow(kapakapi,3) + 3*alpha*csqr(kapakapi)*dkapakapi)*rho2*csqr(rhopi)
		     - 21*alpha*cpow(kapakapi,3)*(drho2.a*csqr(rhopi) + 2*rho2*rhopi*drhopi.a)
		   + 33*(dalpha*cpow(kapakapi,3) + 3*alpha*csqr(kapakapi)*dkapakapi)*pi2*csqr(rho2)
		     + 33*alpha*cpow(kapakapi,3)*(dpi2.a*csqr(rho2) + 2*pi2*rho2*drho2.a)
		   + 8*cpow(kapakapi,3)*dkapakapi*L2*csqr(rho2) 
		     + 2*cpow(kapakapi,4)*(dL2.a*csqr(rho2) + 2*L2*rho2*drho2.a)
		   + 9*csqr(kapakapi)*dkapakapi*cpow(rho2,3)*thelakapi
		     + 9*cpow(kapakapi,3)*csqr(rho2)*drho2.a*thelakapi 
		     + 3*cpow(kapakapi*rho2,3)*dthelakapi
		   - 6*dbeta*cpow(kapakapi,3)*csqr(rho2)*rhopi
		     - 18*beta*csqr(kapakapi)*dkapakapi*csqr(rho2)*rhopi
		     - 6*beta*cpow(kapakapi,3)*(2*rho2*drho2.a*rhopi + csqr(rho2)*drhopi.a)
		   + 6*alpha*dalpha*(9 - 35*kapakapi*psi)
		     - 105*csqr(alpha)*(dkapakapi*psi + kapakapi*dpsi)
		   - 42*(dalpha*beta*csqr(kapakapi) + alpha*dbeta*csqr(kapakapi) 
		         + 2*alpha*beta*kapakapi*dkapakapi)*rho2*rhopi
		     - 42*alpha*beta*csqr(kapakapi)*(drho2.a*rhopi + rho2*drhopi.a)
		   - 30*(2*alpha*dalpha*beta*kapakapi + csqr(alpha)*dbeta*kapakapi
		         + csqr(alpha)*beta*dkapakapi)*rhopi
		     - 30*csqr(alpha)*beta*kapakapi*drhopi.a
		   + (dalpha*kapakapi + alpha*dkapakapi)*(74 - 196.5*kapakapi*psi + 15*theta)*rho2
		     + alpha*kapakapi*(-196.5*dkapakapi*psi - 196.5*kapakapi*dpsi + 15*dtheta)*rho2
		     + alpha*kapakapi*(74 - 196.5*kapakapi*psi + 15*theta)*drho2.a
		   + 2*(kapakapi*dkapakapi*csqr(rho2) + csqr(kapakapi)*rho2*drho2.a)
		        *(16 + 21*theta - 52.5*kapakapi*psi)
		     + csqr(kapakapi*rho2)*(21*dtheta - 52.5*dkapakapi*psi - 52.5*kapakapi*dpsi);	   
      
      for(i=0; i<3; i++)
        dl2tpppipi.b[i] = 60*cpow(alpha,3)*kapakapi*dpi2.b[i] 
		   + 117*csqr(alpha*kapakapi)*(dpi2.b[i]*rho2 + pi2*drho2.b[i])
		   - 66*csqr(alpha*kapakapi)*rhopi*drhopi.b[i]
		   - 21*alpha*cpow(kapakapi,3)*(drho2.b[i]*csqr(rhopi) + 2*rho2*rhopi*drhopi.b[i])
		   + 33*alpha*cpow(kapakapi,3)*(dpi2.b[i]*csqr(rho2) + 2*pi2*rho2*drho2.b[i])
		   + 2*cpow(kapakapi,4)*(dL2.b[i]*csqr(rho2) + 2*L2*rho2*drho2.b[i])
		   + 9*cpow(kapakapi,3)*thelakapi*csqr(rho2)*drho2.b[i]
		   - 6*beta*cpow(kapakapi,3)*(2*rho2*drho2.b[i]*rhopi + csqr(rho2)*drhopi.b[i])
		   - 42*alpha*beta*csqr(kapakapi)*(drho2.b[i]*rhopi + rho2*drhopi.b[i])
		   - 30*csqr(alpha)*beta*kapakapi*drhopi.b[i]
		   + alpha*kapakapi*(74 - 196.5*kapakapi*psi + 15*theta)*drho2.b[i]
		   + csqr(kapakapi)*(16 + 21*theta - 52.5*kapakapi*psi)*2*rho2*drho2.b[i];
		  
		   
      // l2tpprhopi       		   
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
		   
      dl2tpprhopi.a = 24*dalpha*cpow(kapakapi*rhopi,3) 
      		     + 72*alpha*csqr(kapakapi*rhopi)*(dkapakapi*rhopi + kapakapi*drhopi.a) 
		   - 192*alpha*kapakapi*(dalpha*kapakapi + alpha*dkapakapi)*pi2*rhopi
		     - 96*csqr(alpha*kapakapi)*(dpi2.a*rhopi + pi2*drhopi.a)
		   - 48*(dalpha*cpow(kapakapi,3) + 3*alpha*csqr(kapakapi)*dkapakapi)*rhopi*rho2*pi2
		     - 48*alpha*cpow(kapakapi,3)*(drhopi.a*rho2*pi2 + rhopi*drho2.a*pi2 + rhopi*rho2*dpi2.a)
		   - 16*cpow(kapakapi,3)*dkapakapi*rho2*rhopi*L2
		     - 4*cpow(kapakapi,4)*(drho2.a*rhopi*L2 + rho2*drhopi.a*L2 + rho2*rhopi*dL2.a)
		   - 6*(3*csqr(kapakapi)*dkapakapi*thelakapi + cpow(kapakapi,3)*dthelakapi)*csqr(rho2)*rhopi
		     - 6*cpow(kapakapi,3)*thelakapi*(2*rho2*drho2.a*rhopi + csqr(rho2)*drhopi.a)
		   - 4.5*(dbeta*cpow(kapakapi,3) + 3*beta*csqr(kapakapi)*dkapakapi)*rho2*L2
		     - 4.5*beta*cpow(kapakapi,3)*(drho2.a*L2 + rho2*dL2.a)
		   + 12*(dbeta*cpow(kapakapi,3) + 3*beta*csqr(kapakapi)*dkapakapi)*rho2*csqr(rhopi)
		     + 12*beta*cpow(kapakapi,3)*(drho2.a*csqr(rhopi) + 2*rho2*rhopi*drhopi.a)
		   - 34.5*(dalpha*beta*csqr(kapakapi) + alpha*dbeta*csqr(kapakapi) 
		   	   + 2*alpha*beta*kapakapi*dkapakapi)*L2
	             - 34.5*alpha*beta*csqr(kapakapi)*dL2.a
		   + 30*(dalpha*beta*csqr(kapakapi) + alpha*dbeta*csqr(kapakapi) 
		   	   + 2*alpha*beta*kapakapi*dkapakapi)*csqr(rhopi)
	             + 60*alpha*beta*csqr(kapakapi)*rhopi*drhopi.a
		   - 30*(2*alpha*dalpha*beta*kapakapi + csqr(alpha)*dbeta*kapakapi 
		   	 + csqr(alpha)*beta*dkapakapi)*pi2
		     - 30*csqr(alpha)*beta*kapakapi*dpi2.a
		   - 1.5*(dalpha*beta + alpha*dbeta)*(9 - 35*kapakapi*psi)
		     + 52.5*alpha*beta*(dkapakapi*psi + kapakapi*dpsi)
		   + ((dalpha*kapakapi + alpha*dkapakapi)*rhopi + alpha*kapakapi*drhopi.a)
		   	*(-138 + 60*theta + 201*kapakapi*psi)
		     + alpha*kapakapi*(60*dtheta + 201*dkapakapi*psi + 201*kapakapi*dpsi)*rhopi
		   + (2*kapakapi*dkapakapi*rho2*rhopi + csqr(kapakapi)*(drho2.a*rhopi + rho2*drhopi.a))
		   	*(-68 + 87*kapakapi*psi + 12*theta)
		     + csqr(kapakapi)*(87*dkapakapi*psi + 87*kapakapi*dpsi + 12*dtheta)*rho2*rhopi
		   - 4.5*(dbeta*csqr(kapakapi)*thelakapi + 2*beta*kapakapi*dkapakapi*thelakapi
		   	  + beta*csqr(kapakapi)*dthelakapi)*csqr(rho2)
		     - 9*beta*csqr(kapakapi)*thelakapi*rho2*drho2.a
		   + ((dbeta*kapakapi + beta*dkapakapi)*rho2 + beta*kapakapi*drho2.a)
		        *(-17.5 - 7.5*theta + 48*kapakapi*psi)
		     + beta*kapakapi*(48*dkapakapi*psi + 48*kapakapi*dpsi -7.5*dtheta)*rho2;  
      
      for(i=0; i<3; i++)
        dl2tpprhopi.b[i] = 72*alpha*cpow(kapakapi,3)*csqr(rhopi)*drhopi.b[i]
		   - 96*csqr(alpha*kapakapi)*(dpi2.b[i]*rhopi + pi2*drhopi.b[i])
		   - 48*alpha*cpow(kapakapi,3)
		   	*(drhopi.b[i]*rho2*pi2 + rhopi*drho2.b[i]*pi2 +rhopi*rho2*dpi2.b[i])
		   - 4*cpow(kapakapi,4)*(drho2.b[i]*rhopi*L2 + rho2*drhopi.b[i]*L2 + rho2*rhopi*dL2.b[i])
		   - 6*cpow(kapakapi,3)*thelakapi*(2*rho2*drho2.b[i]*rhopi + csqr(rho2)*drhopi.b[i])
		   - 4.5*beta*cpow(kapakapi,3)*(drho2.b[i]*L2 + rho2*dL2.b[i])
		   + 12*beta*cpow(kapakapi,3)*(drho2.b[i]*csqr(rhopi) + 2*rho2*rhopi*drhopi.b[i])
		   - 34.5*alpha*beta*csqr(kapakapi)*dL2.b[i]
		   + 60*alpha*beta*csqr(kapakapi)*rhopi*drhopi.b[i]
		   - 30*csqr(alpha)*beta*kapakapi*dpi2.b[i]
		   + alpha*kapakapi*(-138 + 60*theta + 201*kapakapi*psi)*drhopi.b[i]
		   + csqr(kapakapi)*(-68 + 12*theta + 87*kapakapi*psi)*(drho2.b[i]*rhopi + rho2*drhopi.b[i])
		   - 9*beta*csqr(kapakapi)*thelakapi*rho2*drho2.b[i]
		   + beta*kapakapi*(-17.5 - 7.5*theta + 48*kapakapi*psi)*drho2.b[i];
		   
		     
      // l2tppll
      l2tppll = 42*csqr(alpha*kapakapi)*pi2 - 9*alpha*cpow(kapakapi,3)*csqr(rhopi)
      		   + 15*alpha*cpow(kapakapi,3)*pi2*rho2 
		   + cpow(kapakapi,4)*rho2*L2
		   + 1.5*cpow(kapakapi,3)*thelakapi*csqr(rho2)
		   - 3*beta*cpow(kapakapi,3)*rhopi*rho2
		   - 21*alpha*beta*csqr(kapakapi)*rhopi
		   + alpha*kapakapi*(30 - 73.5*kapakapi*psi)
		   + csqr(kapakapi)*(6 - 24*kapakapi*psi + 10.5*theta)*rho2;
		   
      dl2tppll.a = 84*(alpha*dalpha*csqr(kapakapi) + csqr(alpha)*kapakapi*dkapakapi)*pi2
      		     + 42*csqr(alpha*kapakapi)*dpi2.a
		   - 9*(dalpha*cpow(kapakapi,3) + 3*alpha*csqr(kapakapi)*dkapakapi)*csqr(rhopi)
		     - 18*alpha*cpow(kapakapi,3)*rhopi*drhopi.a
		   + 15*(dalpha*cpow(kapakapi,3) + 3*alpha*csqr(kapakapi)*dkapakapi)*pi2*rho2
		     + 15*alpha*cpow(kapakapi,3)*(dpi2.a*rho2 + pi2*drho2.a)
		   + 4*cpow(kapakapi,3)*dkapakapi*rho2*L2
		     + cpow(kapakapi,4)*(drho2.a*L2 + rho2*dL2.a)
		   + 1.5*(3*csqr(kapakapi)*dkapakapi*thelakapi + cpow(kapakapi,3)*dthelakapi)*csqr(rho2)
		     + 3*cpow(kapakapi,3)*thelakapi*rho2*drho2.a
		   - 3*(dbeta*cpow(kapakapi,3) + 3*beta*csqr(kapakapi)*dkapakapi)*rhopi*rho2
		     - 3*beta*cpow(kapakapi,3)*(drhopi.a*rho2 + rhopi*drho2.a)
		   - 21*(dalpha*beta*csqr(kapakapi) + alpha*dbeta*csqr(kapakapi) 
		   	 + 2*alpha*beta*kapakapi*dkapakapi)*rhopi
		     - 21*alpha*beta*csqr(kapakapi)*drhopi.a
		   + (dalpha*kapakapi + alpha*dkapakapi)*(30 - 73.5*kapakapi*psi)
		     -73.5*alpha*kapakapi*(dkapakapi*psi + kapakapi*dpsi)
		   + (2*kapakapi*dkapakapi*rho2 + csqr(kapakapi)*drho2.a)*(6 - 24*kapakapi*psi + 10.5*theta)
		     + csqr(kapakapi)*(-24*dkapakapi*psi - 24*kapakapi*dpsi + 10.5*dtheta)*rho2;
		     
      for(i=0; i<3; i++)
        dl2tppll.b[i] = 42*csqr(alpha*kapakapi)*dpi2.b[i]
		   - 18*alpha*cpow(kapakapi,3)*rhopi*drhopi.b[i]
		   + 15*alpha*cpow(kapakapi,3)*(dpi2.b[i]*rho2 + pi2*drho2.b[i])
		   + cpow(kapakapi,4)*(drho2.b[i]*L2 + rho2*dL2.b[i])
		   + 3*cpow(kapakapi,3)*thelakapi*rho2*drho2.b[i]
		   - 3*beta*cpow(kapakapi,3)*(drhopi.b[i]*rho2 + rhopi*drho2.b[i])
		   - 21*alpha*beta*csqr(kapakapi)*drhopi.b[i]
		   + csqr(kapakapi)*(6 - 24*kapakapi*psi + 10.5*theta)*drho2.b[i]; 	   
	
		   
      // l2tpprr	   
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
		   
      dl2tpprr.a = 12*(dalpha*cpow(kapakapi,3) + 3*alpha*csqr(kapakapi)*dkapakapi)*pi2*csqr(rhopi)
      		     + 12*alpha*cpow(kapakapi,3)*(dpi2.a*csqr(rhopi) + 2*pi2*rhopi*drhopi.a)
		   + 8*cpow(kapakapi,3)*dkapakapi*csqr(rhopi)*L2
		     + 2*cpow(kapakapi,4)*(2*rhopi*drhopi.a*L2 + csqr(rhopi)*dL2.a)
		   + 3.75*(3*csqr(kapakapi)*dkapakapi*thelakapi + cpow(kapakapi,3)*dthelakapi)*rho2*csqr(rhopi)
		     + 3.75*cpow(kapakapi,3)*thelakapi*(drho2.a*csqr(rhopi) + 2*rho2*rhopi*drhopi.a)
		   - 0.75*(3*csqr(kapakapi)*dkapakapi*thelakapi + cpow(kapakapi,3)*dthelakapi)*pi2*csqr(rho2)
		     - 0.75*cpow(kapakapi,3)*thelakapi*(dpi2.a*csqr(rho2) + 2*pi2*rho2*drho2.a)
		   - 1.5*(kapakapi*dkapakapi*csqr(thelakapi*rho2) + csqr(kapakapi*rho2)*thelakapi*dthelakapi)
		     - 1.5*csqr(kapakapi*thelakapi)*rho2*drho2.a
		   + 6*(dbeta*cpow(kapakapi,3) + 3*beta*csqr(kapakapi)*dkapakapi)*rhopi*(L2 - csqr(rhopi))
		     + 6*beta*cpow(kapakapi,3)*(drhopi.a*(L2 - csqr(rhopi)) + rhopi*(dL2.a - 2*rhopi*drhopi.a))
		   + 24*((dalpha*beta + alpha*dbeta)*csqr(kapakapi) + 2*alpha*beta*kapakapi*dkapakapi)*pi2*rhopi
		     + 24*alpha*beta*csqr(kapakapi)*(dpi2.a*rhopi + pi2*drhopi.a)
		   + (dalpha*kapakapi*pi2 + alpha*dkapakapi*pi2 + alpha*kapakapi*dpi2.a)
		   	*(-35 + 34.5*kapakapi*psi + 15*theta)
		     + alpha*kapakapi*(34.5*dkapakapi*psi + 34.5*kapakapi*dpsi + 15*dtheta)*pi2
		   + 2*(kapakapi*dkapakapi*csqr(rhopi) + csqr(kapakapi)*rhopi*drhopi.a)
		   	*(78.5 - 41.25*kapakapi*psi - 56.25*theta)
		     - csqr(kapakapi*rhopi)*(41.25*dkapakapi*psi + 41.25*kapakapi*dpsi + 56.25*dtheta)
		   + (2*kapakapi*dkapakapi*pi2*rho2 + csqr(kapakapi)*(dpi2.a*rho2 + pi2*drho2.a))
		   	*(-14.5 + 12.75*kapakapi*psi + 5.25*theta)
		     + csqr(kapakapi)*pi2*rho2*(12.75*dkapakapi*psi + 12.75*kapakapi*dpsi + 5.25*dtheta)
		   + (dkapakapi*thelakapi*rho2 + kapakapi*dthelakapi*rho2 + kapakapi*thelakapi*drho2.a)
		   	*(-9.75 + 3.75*theta + 13.5*kapakapi*psi)
	             + kapakapi*thelakapi*rho2*(13.5*dkapakapi*psi + 13.5*kapakapi*dpsi + 3.75*dtheta) 
		   + 7.5*(dbeta*csqr(kapakapi)*thelakapi + 2*beta*kapakapi*dkapakapi*thelakapi
		   	+ beta*csqr(kapakapi)*dthelakapi)*rho2*rhopi
		     + 7.5*beta*csqr(kapakapi)*thelakapi*(drho2.a*rhopi + rho2*drhopi.a)
		   + (dbeta*kapakapi*rhopi + beta*dkapakapi*rhopi + beta*kapakapi*drhopi.a)
		   	*(37 - 67.5*kapakapi*psi - 7.5*theta)
		     - beta*kapakapi*rhopi*(67.5*dkapakapi*psi + 67.5*kapakapi*dpsi + 7.5*dtheta)
		   + 6.75*dtheta - 94.5*(kapakapi*dkapakapi*csqr(psi) + csqr(kapakapi)*psi*dpsi)
		   + 5.25*(dkapakapi*psi + kapakapi*dpsi)*(13 - 5*theta)
		     - 26.25*kapakapi*psi*dtheta; 
		     
      for(i=0; i<3; i++)
        dl2tpprr.b[i] = 12*alpha*cpow(kapakapi,3)*(dpi2.b[i]*csqr(rhopi) + 2*pi2*rhopi*drhopi.b[i])
		   + 2*cpow(kapakapi,4)*(2*rhopi*drhopi.b[i]*L2 + csqr(rhopi)*dL2.b[i])
		   + 3.75*cpow(kapakapi,3)*thelakapi*(drho2.b[i]*csqr(rhopi) + 2*rho2*rhopi*drhopi.b[i])
		   - 0.75*cpow(kapakapi,3)*thelakapi*(dpi2.b[i]*csqr(rho2) + 2*pi2*rho2*drho2.b[i])
		   - 1.5*csqr(kapakapi*thelakapi)*rho2*drho2.b[i]
		   + 6*beta*cpow(kapakapi,3)*(drhopi.b[i]*(L2-csqr(rhopi))+rhopi*(dL2.b[i]-2*rhopi*drhopi.b[i]))
		   + 24*alpha*beta*csqr(kapakapi)*(dpi2.b[i]*rhopi + pi2*drhopi.b[i])
		   + alpha*kapakapi*(-35 + 34.5*kapakapi*psi + 15*theta)*dpi2.b[i]
		   + 2*csqr(kapakapi)*(78.5 - 41.25*kapakapi*psi -56.25*theta)*rhopi*drhopi.b[i]
		   + csqr(kapakapi)*(-14.5 + 12.75*kapakapi*psi + 5.25*theta)*(dpi2.b[i]*rho2 + pi2*drho2.b[i])
		   + kapakapi*thelakapi*(-9.75 + 3.75*theta + 13.5*kapakapi*psi)*drho2.b[i]
		   + 7.5*beta*csqr(kapakapi)*thelakapi*(drho2.b[i]*rhopi + rho2*drhopi.b[i])
		   + beta*kapakapi*(37 - 67.5*kapakapi*psi - 7.5*theta)*drhopi.b[i];
		   
		   
      // L2TPP combined	    
      L2TPP = csqr(kapakapi)*(S12pipi*l2tpppipi + S12rhopi*l2tpprhopi + S12*l2tpprr + S12ll*l2tppll);
      
      dL2TPP.a = 2*L2TPP/kapakapi*dkapakapi
      		 + csqr(kapakapi)*(dS12pipi.a*l2tpppipi + S12pipi*dl2tpppipi.a
		                   + dS12rhopi.a*l2tpprhopi + S12rhopi*dl2tpprhopi.a
				   + dS12.a*l2tpprr + S12*dl2tpprr.a
				   + dS12ll.a*l2tppll + S12ll*dl2tppll.a);
				   
      for(i=0; i<3; i++)
        dL2TPP.b[i] = csqr(kapakapi)*(dS12pipi.b[i]*l2tpppipi + S12pipi*dl2tpppipi.b[i]
				   + dS12rhopi.b[i]*l2tpprhopi + S12rhopi*dl2tpprhopi.b[i]
				   + dS12.b[i]*l2tpprr + S12*dl2tpprr.b[i]
				   + dS12ll.b[i]*l2tppll + S12ll*dl2tppll.b[i]);
				   
      for(i=0; i<2; i++)
        dL2TPP.chi[i] = csqr(kapakapi)*(dS12pipi.chi[i]*l2tpppipi + dS12rhopi.chi[i]*l2tpprhopi
				   + dS12.chi[i]*l2tpprr + dS12ll.chi[i]*l2tppll);
      
      l2tppdone = 1;
		    
    }
    
    // Central Potentials
    if (type == V && TT) {
      *v += gamma*TT*SS*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dSS.chi[i]*RR*Gi;
      dv->a += gamma*TT*SS*(dRR.a*Gi+RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*SS*(dRR.b[i]*Gi+RR*dGi.b[i]);
    } 
    else if (type == sV && TT) {
      *v += gamma*TT*sigsig*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dsigsig.chi[i]*RR*Gi;
      dv->a += gamma*TT*sigsig*(dRR.a*Gi+RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*sigsig*(dRR.b[i]*Gi+RR*dGi.b[i]);
    } 
    else if (type == tV && tautau) {
      *v += gamma*tautau*SS*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dSS.chi[i]*RR*Gi;
      dv->a += gamma*tautau*SS*(dRR.a*Gi+RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*SS*(dRR.b[i]*Gi+RR*dGi.b[i]);
    } 
    else if (type == tsV && tautau) {
      *v += gamma*tautau*sigsig*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dsigsig.chi[i]*RR*Gi;
      dv->a += gamma*tautau*sigsig*(dRR.a*Gi+RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*sigsig*(dRR.b[i]*Gi+RR*dGi.b[i]);
    } 

    // p2V Potentials
    else if (type == p2V && TT) {
      *v += gamma*TT*SS*PiPi*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dSS.chi[i]*PiPi*RR*Gi;
      dv->a += gamma*TT*SS*(dPiPi.a*RR*Gi+PiPi*dRR.a*Gi+PiPi*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*SS*(dPiPi.b[i]*RR*Gi+PiPi*dRR.b[i]*Gi+PiPi*RR*dGi.b[i]);
    } 
    else if (type == sp2V && TT) {
      *v += gamma*TT*sigsig*PiPi*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dsigsig.chi[i]*PiPi*RR*Gi;
      dv->a += gamma*TT*sigsig*(dPiPi.a*RR*Gi+PiPi*dRR.a*Gi+PiPi*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*sigsig*(dPiPi.b[i]*RR*Gi+PiPi*dRR.b[i]*Gi+PiPi*RR*dGi.b[i]);
    } 
    else if (type == tp2V && tautau) {
      *v += gamma*tautau*SS*PiPi*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dSS.chi[i]*PiPi*RR*Gi;
      dv->a += gamma*tautau*SS*(dPiPi.a*RR*Gi+PiPi*dRR.a*Gi+PiPi*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*SS*(dPiPi.b[i]*RR*Gi+PiPi*dRR.b[i]*Gi+PiPi*RR*dGi.b[i]);
    } 
    else if (type == tsp2V && tautau) {
      *v += gamma*tautau*sigsig*PiPi*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dsigsig.chi[i]*PiPi*RR*Gi;
      dv->a += gamma*tautau*sigsig*(dPiPi.a*RR*Gi+PiPi*dRR.a*Gi+PiPi*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*sigsig*(dPiPi.b[i]*RR*Gi+PiPi*dRR.b[i]*Gi+PiPi*RR*dGi.b[i]);
    } 

    // Vp2 Potentials
    else if (type == Vp2 && TT) {
      *v += gamma*TT*SS*PiPiG*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dSS.chi[i]*PiPiG*RR*Gi;
      dv->a += gamma*TT*SS*(dPiPiG.a*RR*Gi+PiPiG*dRR.a*Gi+PiPiG*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*SS*(dPiPiG.b[i]*RR*Gi+PiPiG*dRR.b[i]*Gi+PiPiG*RR*dGi.b[i]);
    } 
    else if (type == sVp2 && TT) {
      *v += gamma*TT*sigsig*PiPiG*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dsigsig.chi[i]*PiPiG*RR*Gi;
      dv->a += gamma*TT*sigsig*(dPiPiG.a*RR*Gi+PiPiG*dRR.a*Gi+PiPiG*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*sigsig*(dPiPiG.b[i]*RR*Gi+PiPiG*dRR.b[i]*Gi+PiPiG*RR*dGi.b[i]);
    } 
    else if (type == tVp2 && tautau) {
      *v += gamma*tautau*SS*PiPiG*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dSS.chi[i]*PiPiG*RR*Gi;
      dv->a += gamma*tautau*SS*(dPiPiG.a*RR*Gi+PiPiG*dRR.a*Gi+PiPiG*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*SS*(dPiPiG.b[i]*RR*Gi+PiPiG*dRR.b[i]*Gi+PiPiG*RR*dGi.b[i]);
    } 
    else if (type == tsVp2 && tautau) {
      *v += gamma*tautau*sigsig*PiPiG*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dsigsig.chi[i]*PiPiG*RR*Gi;
      dv->a += gamma*tautau*sigsig*(dPiPiG.a*RR*Gi+PiPiG*dRR.a*Gi+PiPiG*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*sigsig*(dPiPiG.b[i]*RR*Gi+PiPiG*dRR.b[i]*Gi+PiPiG*RR*dGi.b[i]);
    } 

    // pr2 Potentials
    else if (type == pr2V && TT) {
      *v += gamma*TT*SS*PirPir*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dSS.chi[i]*PirPir*RR*Gi;
      dv->a += gamma*TT*SS*(dPirPir.a*RR*Gi+PirPir*dRR.a*Gi+PirPir*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*SS*(dPirPir.b[i]*RR*Gi+PirPir*dRR.b[i]*Gi+
				PirPir*RR*dGi.b[i]);
    } 
    else if (type == spr2V && TT) {
      *v += gamma*TT*sigsig*PirPir*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dsigsig.chi[i]*PirPir*RR*Gi;
      dv->a += gamma*TT*sigsig*(dPirPir.a*RR*Gi+PirPir*dRR.a*Gi+PirPir*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*sigsig*(dPirPir.b[i]*RR*Gi+PirPir*dRR.b[i]*Gi+
				    PirPir*RR*dGi.b[i]);
    } 
    else if (type == tpr2V && tautau) {
      *v += gamma*tautau*SS*PirPir*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dSS.chi[i]*PirPir*RR*Gi;
      dv->a += gamma*tautau*SS*(dPirPir.a*RR*Gi+PirPir*dRR.a*Gi+PirPir*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*SS*(dPirPir.b[i]*RR*Gi+PirPir*dRR.b[i]*Gi+
				    PirPir*RR*dGi.b[i]);
    } 
    else if (type == tspr2V && tautau) {
      *v += gamma*tautau*sigsig*PirPir*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dsigsig.chi[i]*PirPir*RR*Gi;
      dv->a += gamma*tautau*sigsig*(dPirPir.a*RR*Gi+PirPir*dRR.a*Gi+PirPir*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*sigsig*(dPirPir.b[i]*RR*Gi+PirPir*dRR.b[i]*Gi+
					PirPir*RR*dGi.b[i]);
    } 

    // L2 Potentials
    else if (type == Vl2 && TT) {
      *v += gamma*TT*SS*RR*Gi*LL;
        
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*dSS.chi[i]*RR*Gi*LL;
      dv->a += gamma*TT*SS*(dRR.a*Gi*LL + RR*dGi.a*LL + RR*Gi*dLL.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*SS*(dRR.b[i]*Gi*LL + RR*dGi.b[i]*LL + RR*Gi*dLL.b[i]);    
    }
    else if (type == sVl2 && TT) {
      *v += gamma*TT*sigsig*RR*Gi*LL;
        
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*dsigsig.chi[i]*RR*Gi*LL;
      dv->a += gamma*TT*sigsig*(dRR.a*Gi*LL + RR*dGi.a*LL + RR*Gi*dLL.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*sigsig*(dRR.b[i]*Gi*LL + RR*dGi.b[i]*LL + RR*Gi*dLL.b[i]);    
    }
    else if (type == tVl2 && tautau) {
      *v += gamma*tautau*SS*RR*Gi*LL;
        
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*dSS.chi[i]*RR*Gi*LL;
      dv->a += gamma*tautau*SS*(dRR.a*Gi*LL + RR*dGi.a*LL + RR*Gi*dLL.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*SS*(dRR.b[i]*Gi*LL + RR*dGi.b[i]*LL + RR*Gi*dLL.b[i]);    
    }
    else if (type == tsVl2 && tautau) {
      *v += gamma*tautau*sigsig*RR*Gi*LL;
        
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*dsigsig.chi[i]*RR*Gi*LL;
      dv->a += gamma*tautau*sigsig*(dRR.a*Gi*LL + RR*dGi.a*LL + RR*Gi*dLL.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*sigsig*(dRR.b[i]*Gi*LL + RR*dGi.b[i]*LL + RR*Gi*dLL.b[i]);    
    }
    
    // Spin-Orbit Potentials
    else if (type == Vls && TT) {
      *v += gamma*TT*LS*kapakapi*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dLS.chi[i]*kapakapi*RR*Gi;
      dv->a += gamma*TT*(dLS.a*kapakapi*RR*Gi+LS*dkapakapi*RR*Gi+
			LS*kapakapi*dRR.a*Gi+LS*kapakapi*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*kapakapi*(dLS.b[i]*RR*Gi+LS*dRR.b[i]*Gi+LS*RR*dGi.b[i]);
    } 
    else if (type == tVls && tautau) {
      *v += gamma*tautau*LS*kapakapi*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dLS.chi[i]*kapakapi*RR*Gi;
      dv->a += gamma*tautau*(dLS.a*kapakapi*RR*Gi+LS*dkapakapi*RR*Gi+
			    LS*kapakapi*dRR.a*Gi+LS*kapakapi*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*kapakapi*(dLS.b[i]*RR*Gi+LS*dRR.b[i]*Gi+LS*RR*dGi.b[i]);
    } 

    // L2LS Potentials
    else if (type == Vl2ls && TT) {
      *v += gamma*TT*RR*Gi*L2LS;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*RR*Gi*dL2LS.chi[i];
      dv->a += gamma*TT*(dRR.a*Gi*L2LS + RR*dGi.a*L2LS + RR*Gi*dL2LS.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*(dRR.b[i]*Gi*L2LS + RR*dGi.b[i]*L2LS + RR*Gi*dL2LS.b[i]);
    }
    else if (type == tVl2ls && tautau) {
      *v += gamma*tautau*RR*Gi*L2LS;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*RR*Gi*dL2LS.chi[i];
      dv->a += gamma*tautau*(dRR.a*Gi*L2LS + RR*dGi.a*L2LS + RR*Gi*dL2LS.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*(dRR.b[i]*Gi*L2LS + RR*dGi.b[i]*L2LS + RR*Gi*dL2LS.b[i]);
    }

    // Tensor Potentials
    else if (type == VT && TT) {
      *v += gamma*TT*S12*csqr(kapakapi)*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dS12.chi[i]*csqr(kapakapi)*RR*Gi;
      dv->a += gamma*TT*(dS12.a*csqr(kapakapi)*RR*Gi+
			2*S12*kapakapi*dkapakapi*RR*Gi+
			S12*csqr(kapakapi)*dRR.a*Gi+S12*csqr(kapakapi)*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*csqr(kapakapi)*(dS12.b[i]*RR*Gi+S12*dRR.b[i]*Gi+
					     S12*RR*dGi.b[i]);
    } 
    else if (type == tVT && tautau) {
      *v += gamma*tautau*S12*csqr(kapakapi)*RR*Gi;

      for (i=0; i<2; i++)
	dv->chi[i] += gamma*tautau*dS12.chi[i]*csqr(kapakapi)*RR*Gi;
      dv->a += gamma*tautau*(dS12.a*csqr(kapakapi)*RR*Gi+
			    2*S12*kapakapi*dkapakapi*RR*Gi+
			    S12*csqr(kapakapi)*dRR.a*Gi+S12*csqr(kapakapi)*RR*dGi.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*tautau*csqr(kapakapi)*(dS12.b[i]*RR*Gi+S12*dRR.b[i]*Gi+
						 S12*RR*dGi.b[i]);
    } 

    // TLL Potentials
    else if (type == VTll && TT) {
      *v += gamma*TT*RR*Gi*TLL;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*RR*Gi*dTLL.chi[i];
      dv->a += gamma*TT*(dRR.a*Gi*TLL + RR*dGi.a*TLL + RR*Gi*dTLL.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*(dRR.b[i]*Gi*TLL + RR*dGi.b[i]*TLL + RR*Gi*dTLL.b[i]);
    }
    else if (type == tVTll && tautau) {
      *v += gamma*tautau*RR*Gi*TLL;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*RR*Gi*dTLL.chi[i];
      dv->a += gamma*tautau*(dRR.a*Gi*TLL + RR*dGi.a*TLL + RR*Gi*dTLL.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*(dRR.b[i]*Gi*TLL + RR*dGi.b[i]*TLL + RR*Gi*dTLL.b[i]);
    }
             
    // TPP
    else if (type == VTpp && TT) {
      *v += gamma*TT*RR*Gi*TPP;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*RR*Gi*dTPP.chi[i];
      dv->a += gamma*TT*(dRR.a*Gi*TPP + RR*dGi.a*TPP + RR*Gi*dTPP.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*(dRR.b[i]*Gi*TPP + RR*dGi.b[i]*TPP + RR*Gi*dTPP.b[i]);
    }
    else if (type == tVTpp && tautau) {
      *v += gamma*tautau*RR*Gi*TPP;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*RR*Gi*dTPP.chi[i];
      dv->a += gamma*tautau*(dRR.a*Gi*TPP + RR*dGi.a*TPP + RR*Gi*dTPP.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*(dRR.b[i]*Gi*TPP + RR*dGi.b[i]*TPP + RR*Gi*dTPP.b[i]);
    }
    
    // (p_r v(r) + v(r) p_r)S12(r,p)
    else if (type == prVTrp && TT) {
      *v += gamma*TT*RR*Gi*TRP;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*RR*Gi*dTRP.chi[i];
      dv->a += gamma*TT*(dRR.a*Gi*TRP + RR*dGi.a*TRP + RR*Gi*dTRP.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*(dRR.b[i]*Gi*TRP + RR*dGi.b[i]*TRP + RR*Gi*dTRP.b[i]);
    }
    else if (type == tprVTrp && tautau) {
      *v += gamma*tautau*RR*Gi*TRP;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*RR*Gi*dTRP.chi[i];
      dv->a += gamma*tautau*(dRR.a*Gi*TRP + RR*dGi.a*TRP + RR*Gi*dTRP.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*(dRR.b[i]*Gi*TRP + RR*dGi.b[i]*TRP + RR*Gi*dTRP.b[i]);
    }
    
    // {L^2 S12(p,p)}_H
    else if (type == Vl2Tpp && TT) {
      *v += gamma*TT*RR*Gi*L2TPP;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*TT*RR*Gi*dL2TPP.chi[i];
      dv->a += gamma*TT*(dRR.a*Gi*L2TPP + RR*dGi.a*L2TPP + RR*Gi*dL2TPP.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*TT*(dRR.b[i]*Gi*L2TPP + RR*dGi.b[i]*L2TPP + RR*Gi*dL2TPP.b[i]);
    }
    else if (type == tVl2Tpp && tautau) {
      *v += gamma*tautau*RR*Gi*L2TPP;
      
      for (i=0; i<2; i++)
        dv->chi[i] += gamma*tautau*RR*Gi*dL2TPP.chi[i];
      dv->a += gamma*tautau*(dRR.a*Gi*L2TPP + RR*dGi.a*L2TPP + RR*Gi*dL2TPP.a);
      for (i=0; i<3; i++)
        dv->b[i] += gamma*tautau*(dRR.b[i]*Gi*L2TPP + RR*dGi.b[i]*L2TPP + RR*Gi*dL2TPP.b[i]);
    }
    
    
    // Coulomb Potential
    else if (type == VC && TT && G1->xi == 1 && G2->xi == 1) {
      vcoul = 1.0/csqrt(2.0*alpha)*zcoulomb(0.5*rho2/alpha);
      dvcoul.a = -0.5*dalpha/alpha*vcoul +
	1.0/csqrt(2*alpha)*dzcoulomb(0.5*rho2/alpha)*
	(0.5*drho2.a/alpha-0.5*rho2/csqr(alpha)*dalpha);
      for (i=0; i<3; i++)
	dvcoul.b[i] = 1.0/csqrt(2*alpha)*dzcoulomb(0.5*rho2/alpha)*0.5*drho2.b[i]/alpha;

      *v += gamma*TT*SS*vcoul*RR;
 
      for (i=0; i<2; i++)
	dv->chi[i] += gamma*TT*dSS.chi[i]*vcoul*RR;
      dv->a += gamma*TT*SS*(dvcoul.a*RR+vcoul*dRR.a);
      for (i=0; i<3; i++)
	dv->b[i] += gamma*TT*SS*(dvcoul.b[i]*RR+vcoul*dRR.b[i]);
    }
  }
}


void calcgradPotential(const Interaction *P,
		       const SlaterDet* Q, const SlaterDetAux* X, 
		       const gradSlaterDetAux* dX,
		       gradSlaterDet* dv)
{
  gradTwoBodyOperator gop_tb_pot = {opt: 0, par: P, me: gtb_pot};

  calcgradSlaterDetTBME(Q, X, dX, &gop_tb_pot, dv);
}


void calcgradPotentialod(const Interaction *P,
			 const SlaterDet* Q, const SlaterDet* Qp,
			 const SlaterDetAux* X, 
			 const gradSlaterDetAux* dX,
			 gradSlaterDet* dv)
{
  gradTwoBodyOperator gop_tb_pot = {opt: 0, par: P, me: gtb_pot};

  calcgradSlaterDetTBMEod(Q, Qp, X, dX, &gop_tb_pot, dv);
}


void calcgradPotentialrowcol(const Interaction *P,
			     const SlaterDet* Q, const SlaterDetAux* X, 
			     const gradSlaterDetAux* dX,
			     gradSlaterDet* dv,
			     int k, int l)
{
  gradTwoBodyOperator gop_tb_pot = {opt: 0, par: P, me: gtb_pot};

  calcgradSlaterDetTBMErowcol(Q, X, dX, &gop_tb_pot, dv, k, l);
}


void calcgradPotentialodrowcol(const Interaction *P,
			       const SlaterDet* Q, const SlaterDet* Qp,
			       const SlaterDetAux* X, 
			       const gradSlaterDetAux* dX,
			       gradSlaterDet* dv,
			       int k, int l)
{
  gradTwoBodyOperator gop_tb_pot = {opt: 0, par: P, me: gtb_pot};

  calcgradSlaterDetTBMEodrowcol(Q, Qp, X, dX, &gop_tb_pot, dv, k, l);
}

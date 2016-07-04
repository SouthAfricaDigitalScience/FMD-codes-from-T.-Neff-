/**

  \file Interaction.h

  Nucleon-Nucleon interaction

  (c) 2003 Thomas Neff
  
  Changes (HH):
  04/09/02
  + Added support for S12(p,p) and (p_r v(r) + v(r)p_r)S12(r,p) interaction terms.

  04/07/21
  + Added support for L2, L2LS, and S12(L,L) interaction terms.


*/

#ifndef _INTERACTION_H
#define _INTERACTION_H

#define LABELSIZE 10

// maximum number of labeled interaction components
#define MAXINTERACTIONS 32	       

/// implemented interaction types
typedef enum {V, sV, tV, tsV, VT, tVT,
              p2V, sp2V, tp2V, tsp2V, 
              pr2V, spr2V, tpr2V, tspr2V,
              Vp2, sVp2, tVp2, tsVp2,
	      Vl2, sVl2, tVl2, tsVl2,
              Vls, tVls,
	      Vl2ls, tVl2ls,
	      VTll, tVTll,
	      VTpp, tVTpp,
	      Vl2Tpp, tVl2Tpp,
	      prVTrp, tprVTrp,
	      VC} InteractionType;

/// component of interaction
typedef struct {
  InteractionType type;		///< the type
  double gamma;			///< strength gamma [fm^-1]
  double kappa;			///< range kappa [fm^2]
  int ilabel;			///< index for labeling
} InteractionComponent;


typedef struct {
  char name[161];               ///< Interaction name
  char scname[80];              ///< Scaling name
  int ncomponents;
  InteractionComponent* c;
  int n;
  char **label;
  int mescaling;                ///< are matrix elements to be scaled ?
  double *scale;
  int cm;
  int spinorbit;
  int tensor;
  int momentump2;
  int momentumpr2;
  int l2;
  int l2ls;
  int tll;
  int tpp;
  int l2tpp;
  int prtrp;
} Interaction;


int readInteractionfromFile(Interaction* pot, char* filename);

#endif


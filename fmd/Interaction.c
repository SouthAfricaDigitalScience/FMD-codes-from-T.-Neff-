/**

  \file Interaction.c

  Nucleon-Nucleon interaction


  (c) 2003 Thomas Neff
  
  Changes (HH):
  04/12/21
  + Added support for {L2S12(p,p)}_H interaction term.
  
  04/09/02
  + Added support for S12(p,p) and (p_r v(r) + v(r)p_r)S12(r,p) interaction terms.
 
  04/07/21
  + Added support L2, L2LS, and S12(L,L) interaction terms.

  allow to kinds of scaling:
    :: scale the interaction directly
    :  keep the interaction unscaled but read scaling parameters for later scaling of matrix elements

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Interaction.h"
#include "misc/physics.h"
#include "misc/utils.h"


int InteractionComponentCmp(const InteractionComponent* i1, 
			    const InteractionComponent* i2)
{
  return (i1->kappa == i2->kappa ? 0 : (i1->kappa < i2->kappa ? -1 : +1));
}


#define BUFSIZE 255
#define MAXCOMPONENTS 2000

int readInteractionfromFile(Interaction* P, char* fname)
{
  FILE* fp;
  char buf[BUFSIZE];

  InteractionComponent dummyc[MAXCOMPONENTS];
  char dummylabel[MAXCOMPONENTS][LABELSIZE];

  // decompose intname into fname for interaction and possible scaling file
  char *fnintp, fnint[80] = "";
  char *fnscalep, fnscale[80] = "";
  int intscaling = 0;   // scale interaction directly
  int mescaling = 0;    // for scaling matrix elements later

  fnintp = strtok(fname, ":"); 
  strncpy(fnint, fnintp, 80);
  
  fnscalep = strtok(NULL, ":"); 
  if (fnscalep) {
    if (*fnscalep == '*') {
      intscaling = 1;
      strncpy(fnscale, fnscalep+1, 80);
    } else {
      mescaling = 1;
      strncpy(fnscale, fnscalep, 80);
    }
  }

  // reading interaction

  if (!(fp = fopen(fnint, "r"))) {
    fprintf(stderr, "couldn't open %s for reading\n", fnint);
    return -1;
  }

  fprintf(stderr, "... reading Interaction from file %s\n", fnint);

  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Interaction", 12) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Interaction ...>\n");
    return -1;
  }

  char Pname[80];
  sscanf(buf, "<Interaction %s>", Pname);
  strcpy(P->name, stripstr(Pname, ">"));

  char typename[BUFSIZE], label[BUFSIZE];
  double gamma, kappa;

  strcpy(dummylabel[0], "sum\n");

  int i=-1;
  int ilabel=0;

  while (1) {
    fgets(buf, BUFSIZE, fp);
    if (feof(fp)) {
      fprintf(stderr, "did't find </Interaction>\n");
      return -1;
    }
    if (!strncmp(buf, "</Interaction>", 14))
      break;
    if (buf[0] == '#' || buf[0] == '\n')
      continue;

    i++;
    label[0] = '\0';
    if (sscanf(buf, "%s %lf %lf %s", typename, &gamma, &kappa, label) < 3) {
      fprintf(stderr, "malformed line\n>>%s<<\n",
	      stripstr(buf,"\n"));
      return -1;
    }

    if (label[0]) {
      ilabel++;
      strncpy(dummylabel[ilabel], label, LABELSIZE);
    }

    if (i==0 && ilabel==0) {
      fprintf(stderr, "no label defined in first line\n>>%s<<\n",
	      stripstr(buf,"\n"));
      return -1;
    }

    dummyc[i].gamma = gamma;
    dummyc[i].kappa = kappa;
    dummyc[i].ilabel = ilabel;

    if (!strcmp(typename,"V"))
      dummyc[i].type = V;
    else if (!strcmp(typename,"sV"))
      dummyc[i].type = sV;
    else if (!strcmp(typename,"tV"))
      dummyc[i].type = tV;
    else if (!strcmp(typename,"tsV"))
      dummyc[i].type = tsV;
    else if (!strcmp(typename,"p2V"))
      dummyc[i].type = p2V;
    else if (!strcmp(typename,"sp2V"))
      dummyc[i].type = sp2V;
    else if (!strcmp(typename,"tp2V"))
      dummyc[i].type = tp2V;
    else if (!strcmp(typename,"tsp2V"))
      dummyc[i].type = tsp2V;
    else if (!strcmp(typename,"pr2V"))
      dummyc[i].type = pr2V;
    else if (!strcmp(typename,"spr2V"))
      dummyc[i].type = spr2V;
    else if (!strcmp(typename,"tpr2V"))
      dummyc[i].type = tpr2V;
    else if (!strcmp(typename,"tspr2V"))
      dummyc[i].type = tspr2V;
    else if (!strcmp(typename,"Vp2"))
      dummyc[i].type = Vp2;
    else if (!strcmp(typename,"sVp2"))
      dummyc[i].type = sVp2;
    else if (!strcmp(typename,"tVp2"))
      dummyc[i].type = tVp2;
    else if (!strcmp(typename,"tsVp2"))
      dummyc[i].type = tsVp2;
    else if (!strcmp(typename,"Vl2"))
      dummyc[i].type = Vl2;
    else if (!strcmp(typename,"sVl2"))
      dummyc[i].type = sVl2;
    else if (!strcmp(typename,"tVl2"))
      dummyc[i].type = tVl2;
    else if (!strcmp(typename,"tsVl2"))
      dummyc[i].type = tsVl2;
    else if (!strcmp(typename,"Vls"))
      dummyc[i].type = Vls;
    else if (!strcmp(typename,"tVls"))
      dummyc[i].type = tVls;
    else if (!strcmp(typename,"Vl2ls"))
      dummyc[i].type = Vl2ls;
    else if (!strcmp(typename,"tVl2ls"))
      dummyc[i].type = tVl2ls;
    else if (!strcmp(typename,"VTll"))
      dummyc[i].type = VTll;
    else if (!strcmp(typename,"tVTll"))
      dummyc[i].type = tVTll;
    else if (!strcmp(typename,"VTpp"))
      dummyc[i].type = VTpp;
    else if (!strcmp(typename,"tVTpp"))
      dummyc[i].type = tVTpp;
    else if (!strcmp(typename,"Vl2Tpp"))
      dummyc[i].type = Vl2Tpp;
    else if (!strcmp(typename,"tVl2Tpp"))
      dummyc[i].type = tVl2Tpp;
    else if (!strcmp(typename,"VT"))
      dummyc[i].type = VT;
    else if (!strcmp(typename,"tVT"))
      dummyc[i].type = tVT;
    else if (!strcmp(typename,"prVTrp"))
      dummyc[i].type = prVTrp;
    else if (!strcmp(typename,"tprVTrp"))
      dummyc[i].type = tprVTrp;
    else if (!strcmp(typename,"VC")) {
      dummyc[i].type = VC;
      dummyc[i].gamma = alpha;		// finestructure constant
    }
    else {
      fprintf(stderr, "unknown interaction component %s\n",
	      typename);
      return -1;
    }
  }
  fclose(fp);

  qsort(dummyc, i+1, sizeof(InteractionComponent),
	(int (*)(const void*, const void*))InteractionComponentCmp);

  P->ncomponents = i+1;
  P->c = (InteractionComponent*) malloc(P->ncomponents*sizeof(InteractionComponent));
  for (i=0; i<P->ncomponents; i++)
    P->c[i] = dummyc[i];

  P->n = ilabel+1;
  if (P->n > MAXINTERACTIONS) {
    fprintf(stderr, "more then %d labeled interaction components,\nincrease MAXINTERACTIONS in ${FMD}/src/fmd/Interaction.h\n", MAXINTERACTIONS);
    return -1;
  }

  P->label = malloc(P->n*sizeof(char *));
  for (i=0; i<P->n; i++) {
    P->label[i] = (char *) malloc(LABELSIZE*sizeof(char));
    strncpy(P->label[i], dummylabel[i], LABELSIZE);
  }

  P->cm = 0;
  P->spinorbit = 0; P->tensor = 0; P->momentump2 = 0; P->momentumpr2 = 0;
  P->l2 = 0; P->l2ls = 0; P->tll = 0; P->tpp = 0; P->prtrp = 0; P->l2tpp = 0;

  for (i=0; i<P->ncomponents; i++) {
   if (P->c[i].type == Vls || P->c[i].type == tVls)
      P->spinorbit = 1;
   if (P->c[i].type == Vl2ls || P->c[i].type == tVl2ls)
      P->l2ls = 1;
   if (P->c[i].type == VT || P->c[i].type == tVT)
     P->tensor = 1;
   if (P->c[i].type == VTll || P->c[i].type == tVTll)
     P->tll = 1;
   if (P->c[i].type == VTpp || P->c[i].type == tVTpp)
     P->tpp = 1;
   if (P->c[i].type == Vl2Tpp || P->c[i].type == tVl2Tpp)
     P->l2tpp = 1;
   if (P->c[i].type == prVTrp || P->c[i].type == tprVTrp)
     P->prtrp = 1;
   if (P->c[i].type >= Vl2 && P->c[i].type <= tsVl2)
     P->l2 = 1;
   if (P->c[i].type >= p2V && P->c[i].type <= tsp2V || 
       P->c[i].type >= Vp2 && P->c[i].type <= tsVp2)
     P->momentump2 = 1;
   if (P->c[i].type >= pr2V && P->c[i].type <= tspr2V)
     P->momentumpr2 = 1;
  }


  // read scaling parameters

  if (intscaling || mescaling) {

    P->scale = malloc(P->n*sizeof(double));
    for (i=0; i<P->n; i++)
      P->scale[i] = 1.0;

    if (!(fp = fopen(fnscale, "r"))) {
      fprintf(stderr, "couldn't open %s for reading\n", fnscale);
      return -1;
    }

    fprintf(stderr, "... reading Scaling from file %s\n", fnscale);

    do
      fgets(buf, BUFSIZE, fp);
    while (strncmp(buf, "<Scaling", 8) && !feof(fp));
    if (feof(fp)) {
      fprintf(stderr, "did not find <Scaling ...>\n");
      return -1;
    }

    char Pname[80], Sname[80];
    sscanf(buf, "<Scaling %s %s>", Pname, Sname);
    if (strcmp(P->name, Pname)) {
      fprintf(stderr, "Interaction does not match\n");
      return -1;
    }
    strcpy(P->scname, stripstr(Sname, ">"));

    double scale;
    while (1) {
      fgets(buf, BUFSIZE, fp);

      if (feof(fp)) {
        fprintf(stderr, "did't find </Scaling>\n");
        return -1;
      }
      if (!strncmp(buf, "</Scaling>", 10))
        break;
      if (buf[0] == '#' || buf[0] == '\n')
        continue;

      if (sscanf(buf, "%s %lf", label, &scale) < 2) {
        fprintf(stderr, "malformed line\n>>%s<<\n",
                stripstr(buf, "\n"));
        return -1;
      }

      for (i=0; i<P->n; i++) {
        if (!strcmp(P->label[i], label)) {
          P->scale[i] = scale;
          break;
        }
      }


    }

    fprintf(stderr, "Scaling parameters: \n");
    for (i=1; i<P->n; i++) {
      fprintf(stderr, "%s %8.3f\n", P->label[i], P->scale[i]);
    }
  }

  // mescaling: 
  if (mescaling) {
    fprintf(stderr, "scaling parameters for matrix elements\n");
    P->mescaling = 1;
  } else
    P->mescaling = 0;

  // intscaling
  if (intscaling) {
    fprintf(stderr, "scaling interaction\n");
    strcat(P->name, ":");
    strcat(P->name, P->scname);

    // scale interaction components
    for (i=0; i<P->ncomponents; i++) {
      P->c[i].gamma *= P->scale[P->c[i].ilabel];
    }
  }

  return 0;
}


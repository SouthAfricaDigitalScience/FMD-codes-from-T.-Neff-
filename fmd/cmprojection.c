/** 

   \file cmprojection.c

   Gridpoints used for linear momentum projection 

   (c) 2004-2011 Thomas Neff

*/

// This is quite messy - we are using cmintegrationpara not really
// consistently


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "SlaterDet.h"
#include "CenterofMass.h"
#include "Projection.h"

#include "misc/physics.h"

#include "numerics/gaussquad.h"


#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))


double linear[1][3] = {{ 0., 0., 1.}};

double tetrahedron[4][3] = {{0., 0., 1.}, 
			    {0., 0.942809, -0.333333}, 
			    {-0.816497, -0.471405, -0.333333}, 
			    {0.816497, -0.471405, -0.333333}};

double octahedron[6][3] = {{0., 0., 1.},
			   {1., 0., 0.},
			   {0., 1., 0.},
			   {0., 0., -1.},
			   {-1., 0., 0.},
			   {0, -1., 0.}};

double cube[8][3] = {{  0.57735,  0.57735,  0.57735},
		     {  0.57735,  0.57735, -0.57735},
		     {  0.57735, -0.57735,  0.57735},
		     {  0.57735, -0.57735, -0.57735},
		     { -0.57735,  0.57735,  0.57735},
		     { -0.57735,  0.57735, -0.57735},
		     { -0.57735, -0.57735,  0.57735},
		     { -0.57735, -0.57735, -0.57735}};

double (*polyhedron[])[3] = { NULL,
			      linear,
			      NULL,
			      NULL,
			      tetrahedron,
			      NULL,
			      octahedron,
			      NULL,
			      cube };
		     


void initCMintegration(cmintegrationpara* cmpara, const char* projpar)
{
  char projparcpy[strlen(projpar)];
  strcpy(projparcpy, projpar);
  char* c = projparcpy;

  // default is CMNone
  cmpara->type = CMNone;
  cmpara->n = 1;

  c = strtok(c, "-");
  while (c && strncmp(c, "cm", 2)) { 
    c=strtok(NULL, "-"); 
  };
  
  if (c) {
    c=strtok(NULL, "-");
    if (!strncmp(c, "none", 4)) {
      cmpara->type = CMNone;
      cmpara->n = 1;
    } else if (!strncmp(c, "simple", 6)) {
      cmpara->type = CMSimple;
      cmpara->n = 1;
    } else {
      int nr = atoi(c);
      c=strtok(NULL, "-");
      if (!strncmp(c, "tet", 3)) {
	cmpara->type = CMPoly;
	cmpara->poly.nr = nr;
	cmpara->poly.npoly = 4;
	cmpara->n = 4*nr;
      } else if (!strncmp(c, "oct", 3)) {
	cmpara->type = CMPoly;
	cmpara->poly.nr = nr;
	cmpara->poly.npoly = 6;
	cmpara->n = 6*nr;
      } else if (!strncmp(c, "cbe", 3)) {
	cmpara->type = CMPoly;
	cmpara->poly.nr = nr;
	cmpara->poly.npoly = 8;
	cmpara->n = 8*nr;
      } else {
	cmpara->type = CMProd;
	cmpara->prod.nr = nr;
	cmpara->prod.ntheta = atoi(c);
	c = strtok(NULL, "-");
	cmpara->prod.nphi = atoi(c);
	cmpara->n = nr*cmpara->prod.ntheta*cmpara->prod.nphi;
      }
    }
  }

  cmpara->x = malloc(3*cmpara->n*sizeof(double));
  cmpara->w = malloc(cmpara->n*sizeof(double));

  return 0;
}


double _estimateacm(const SlaterDet* Q)
{
  // do we have to initialize workspace first ?
  // if (X.Gaux == NULL)
  //   initSlaterDetAux(Q, &X);

  SlaterDetAux X;
  initSlaterDetAux(Q, &X);

  double tcm;

  calcSlaterDetAux(Q, &X);
  calcTCM(Q, &X, &tcm);
  freeSlaterDetAux(&X);
  
  return (0.75/(tcm*(mproton*Q->Z+mneutron*Q->N)));
}


void _initcmintegration(const Projection* P, 
                        double alpha,
                        cmintegrationpara* par)
{
  if (P->cm == CMNone) {
    par->n = 1;
    par->x = malloc(3*sizeof(double));
    par->w = malloc(1*sizeof(double));
    par->x[0][0] = 0.0; par->x[0][1] = 0.0; par->x[0][2] = 0.0; 
    par->w[0] = 1.0;
  } 

  else if (P->cm == CMSimple) {
    par->n = 1;
    par->x = malloc(3*sizeof(double));
    par->w = malloc(1*sizeof(double));
    par->x[0][0] = 0.0; par->x[0][1] = 0.0; par->x[0][2] = 0.0; 
    par->w[0] = 0.125*pow(M_PI*alpha,-1.5);
  }

  else if (P->cm == CMProd) {
    int nr=P->cmprod.nr; 
    int nx=P->cmprod.ntheta; int nphi=P->cmprod.nphi;
    
    double pointr[nr], weightr[nr];
    GaussSqrHermitePoints(nr, alpha, pointr, weightr);

    double pointx[nx], weightx[nx];
    GaussLegendrePoints(nx, -1.0, 1.0, pointx, weightx);

    double pointphi[nphi], weightphi[nphi];
    ShiftedPeriodicTrapezoidalPoints(nphi, 0, 2*M_PI, 0.00,
				     pointphi, weightphi);

    par->n = nr*nx*nphi;
    par->x = malloc(3*par->n*sizeof(double));
    par->w = malloc(par->n*sizeof(double));

    int i=0;
    int ir, ix, iphi;
    for (ir=0; ir<nr; ir++)
      for (ix=0; ix<nx; ix++)
	for (iphi=0; iphi<nphi; iphi++) {
	  par->x[i][0] = pointr[ir]*sqrt(1-SQR(pointx[ix]))*cos(pointphi[iphi]);
	  par->x[i][1] = pointr[ir]*sqrt(1-SQR(pointx[ix]))*sin(pointphi[iphi]);
	  par->x[i][2] = pointr[ir]*pointx[ix];
	  par->w[i] = 1/CUB(2*M_PI)*
	    weightr[ir]*weightx[ix]*weightphi[iphi];
	  i++;
	}
  }

  else if (P->cm == CMPoly) {
    int nr=P->cmpoly.nr; 
    int npoly=P->cmpoly.npoly;
    
    double pointr[nr], weightr[nr];
    GaussSqrHermitePoints(nr, alpha, pointr, weightr);

    par->n = nr*npoly;
    par->x = malloc(3*par->n*sizeof(double));
    par->w = malloc(par->n*sizeof(double));

    int i=0;	
    int ir, ipoly, j;	
    for (ir=0; ir<nr; ir++)
      for (ipoly=0; ipoly<npoly; ipoly++) {
	for (j=0; j<3; j++)
	  par->x[i][j] = pointr[ir]*polyhedron[npoly][ipoly][j];
	par->w[i] = 1/CUB(2*M_PI)*
	  weightr[ir]*(4*M_PI/npoly);
	i++;
      }
  }
}


void initcmintegration(const Projection* P, 
		       const SlaterDet* Q, const SlaterDet* Qp,
		       cmintegrationpara* par)
{
  alpha = 0.5/(_estimateacm(Q)+_estimateacm(Qp));

  _initcmintegration(P, alpha, par);
}


void reinitcmintegration(const SlaterDet* Q, const SlaterDet* Qp,
			 cmintegrationpara* par)
{
  alpha = 0.5/(_estimateacm(Q)+_estimateacm(Qp));

  if (par->type == CMNone) {
    par->w[0] = 1.0;
    par->x[0][0] = 0.0; par->x[0][1] = 0.0; par->x[0][2] = 0.0; 
  } 

  else if (par->type == CMSimple) {
    par->w[0] = 0.125*pow(M_PI*alpha,-1.5);
    par->x[0][0] = 0.0; par->x[0][1] = 0.0; par->x[0][2] = 0.0; 
  }

  else if (par->type == CMProd) {
    int nr=par->prod.nr; 
    int nx=par->prod.ntheta; int nphi=par->prod.nphi;
    
    double pointr[nr], weightr[nr];
    GaussSqrHermitePoints(nr, alpha, pointr, weightr);

    double pointx[nx], weightx[nx];
    GaussLegendrePoints(nx, -1.0, 1.0, pointx, weightx);

    double pointphi[nphi], weightphi[nphi];
    ShiftedPeriodicTrapezoidalPoints(nphi, 0, 2*M_PI, 0.00,
				     pointphi, weightphi);

    int i=0;
    int ir, ix, iphi;
    for (ir=0; ir<nr; ir++)
      for (ix=0; ix<nx; ix++)
	for (iphi=0; iphi<nphi; iphi++) {
	  par->x[i][0] = pointr[ir]*sqrt(1-SQR(pointx[ix]))*cos(pointphi[iphi]);
	  par->x[i][1] = pointr[ir]*sqrt(1-SQR(pointx[ix]))*sin(pointphi[iphi]);
	  par->x[i][2] = pointr[ir]*pointx[ix];
	  par->w[i] = 1/CUB(2*M_PI)*
	    weightr[ir]*weightx[ix]*weightphi[iphi];
	  i++;
	}
  }

  else if (par->type == CMPoly) {
    int nr=par->poly.nr; 
    int npoly=par->poly.npoly;
    
    double pointr[nr], weightr[nr];
    GaussSqrHermitePoints(nr, alpha, pointr, weightr);

    int i=0;	
    int ir, ipoly, j;	
    for (ir=0; ir<nr; ir++)
      for (ipoly=0; ipoly<npoly; ipoly++) {
	for (j=0; j<3; j++)
	  par->x[i][j] = pointr[ir]*polyhedron[npoly][ipoly][j];
	par->w[i] = 1/CUB(2*M_PI)*
	  weightr[ir]*(4*M_PI/npoly);
	i++;
      }
  }  
}


void freecmintegration(cmintegrationpara* par)
{
  free(par->x);
  free(par->w);
}


void getcmintegrationpoint(int i, const cmintegrationpara* par,
			   double x[3], double* w)
{
  assert(i < par->n);

  int j;
  for (j=0; j<3; j++)
    x[j] = par->x[i][j];
  
  *w = par->w[i];

  // fprintf(stderr, "cm integration point: (%8.5f, %8.5f, %8.5f) weight %8.5f\n", 
  //  	     x[0], x[1], x[2], *w);
}

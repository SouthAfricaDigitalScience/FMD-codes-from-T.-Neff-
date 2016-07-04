/**

  \file Parameterization.h

  Parametrizations of Slater determinants


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATION_H
#define _PARAMETERIZATION_H

#include <stdio.h>

#include "SlaterDet.h"
#include "gradSlaterDet.h"


/// parameters for Parameterization
typedef struct {
  char name[80];       	///< name/description
  int A;		///< number of nucleons
  int Z;		///< number of protons
  int N;		///< number of neutrons

  int n;		///< number of real parameters
  double* x;		///< array of real parameters

  int ngauss;		///< number of gaussians
  void* internals;	///< internals needed for specific parameterization
} Para;
  

/// (virtual) Parameterization of a Slater determinant.
typedef struct {
  char* name;		///< its name

  /// read parameters from file
  int (*Pararead)(FILE* fp, Para* q);
  /// write parameters to file
  int (*Parawrite)(FILE* fp, const Para* q);

  /// clone parameters
  void (*Paraclone)(const Para* q, Para* qp);

  /// init Slater determinant with the parameters
  void (*ParainitSlaterDet)(const Para* q, SlaterDet* Q);
  /// convert paramaters to Slater determinant
  void (*ParatoSlaterDet)(const Para* q, SlaterDet* Q);
  /// project Slater determinant gradient to paramater gradient
  void (*ParaprojectgradSlaterDet)(const Para* q, const gradSlaterDet* dQ, 
				   double* dq); 
} Parameterization;


void shakePara(Para* q, double magnitude);

int readParafromFile(Parameterization* P, Para* q, const char* fname);

#endif

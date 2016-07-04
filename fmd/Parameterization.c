/**

  \file Parameterization.c

  Parametrizations of Slater determinants


  (c) 2003 Thomas Neff

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"

#include "ParameterizationFMD.h"
#include "ParameterizationFMDr.h"
#include "ParameterizationFMDS.h"
#include "ParameterizationFMDreal.h"
#include "ParameterizationAMD.h"
#include "ParameterizationAMDd.h"
#include "ParameterizationAMDS.h"
#include "ParameterizationAMDA.h"
#include "ParameterizationAMDAS.h"
#include "ParameterizationAlpha.h"
#include "ParameterizationAlphaA.h"
#include "ParameterizationAlphaC.h"
#include "ParameterizationCluster.h"
#include "ParameterizationCoreFMD.h"
#include "ParameterizationClusterFMD.h"
#include "ParameterizationClustersFMD.h"
#include "ParameterizationFMDd3h.h"
#include "ParameterizationFMDvxz.h"

#include "misc/utils.h"


// random number
inline static double ranmag(double magnitude)
{
  return magnitude*(double)(rand()-RAND_MAX/2)/(double)(RAND_MAX/2);
}

// primitive shaking 
void shakePara(Para* q, double magnitude)
{
  int i;

  for (i=0; i<q->n; i++)
    q->x[i] = q->x[i]*(1+ranmag(magnitude)) + ranmag(magnitude);
}
  

#define BUFSIZE 255

int readParafromFile(Parameterization* P, Para* q, const char* fname)
{
  FILE* fp;

  if (!(fp = fopen(fname, "r"))) {  
    fprintf(stderr, "couldn't open %s for reading\n", fname);
    return -1;
  }

  fprintf(stderr, "... reading Parameters from file %s\n", fname);

  char buf[BUFSIZE];

  do
    fgets(buf, BUFSIZE, fp);
  while (strncmp(buf, "<Parameterization", 16) && !feof(fp));
  if (feof(fp)) {
    fprintf(stderr, "did't find <Parameterization ...>\n");
    return -1;
  }

  char Pname[80];
  sscanf(buf, "<Parameterization %s>", Pname);
  stripstr(Pname, ">");
  
  if (!strcmp(Pname, "FMD"))
    *P = ParameterizationFMD;
  else if (!strcmp(Pname, "FMDr"))
    *P = ParameterizationFMDr;
  else if (!strcmp(Pname, "FMDS"))
    *P = ParameterizationFMDS;
  else if (!strcmp(Pname, "FMDreal"))
    *P = ParameterizationFMDreal;
  else if (!strcmp(Pname, "AMD"))
    *P = ParameterizationAMD;
  else if (!strcmp(Pname, "AMDd"))
    *P = ParameterizationAMDd;
  else if (!strcmp(Pname, "AMDS"))
    *P = ParameterizationAMDS;
  else if (!strcmp(Pname, "AMDA"))
    *P = ParameterizationAMDA;
  else if (!strcmp(Pname, "AMDAS"))
    *P = ParameterizationAMDAS;
  else if (!strcmp(Pname, "Alpha"))
    *P = ParameterizationAlpha;
  else if (!strcmp(Pname, "AlphaA"))
    *P = ParameterizationAlphaA;
  else if (!strcmp(Pname, "AlphaC"))
    *P = ParameterizationAlphaC;
  else if (!strcmp(Pname, "Cluster"))
    *P = ParameterizationCluster;
  else if (!strcmp(Pname, "CoreFMD"))
    *P = ParameterizationCoreFMD;
  else if (!strcmp(Pname, "ClusterFMD"))
    *P = ParameterizationClusterFMD;
  else if (!strcmp(Pname, "ClustersFMD"))
    *P = ParameterizationClustersFMD;
  else if (!strcmp(Pname, "FMDd3h"))
    *P = ParameterizationFMDd3h;
  else if (!strcmp(Pname, "FMDvxz"))
    *P = ParameterizationFMDvxz;
  else {
    fprintf(stderr, "Parameterization %s not implemented\n", Pname);
    return -1;
  }

  return P->Pararead(fp, q);

}


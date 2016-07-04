/**

  \file ParameterizationAMDA.h

  Parametrization of SlaterDet.

  common fixed width for all single-particle states


  (c) 2003 Thomas Neff

*/


#ifndef _PARAMETERIZATIONAMDA_H
#define _PARAMETERIZATIONAMDA_H


#include "SlaterDet.h"
#include "gradSlaterDet.h"

#include "Parameterization.h"


extern Parameterization ParameterizationAMDA;


int AMDAread(FILE* fp, Para* q);
int AMDAwrite(FILE* fp, const Para* q);

void AMDAclone(const Para* q, Para* qp);

void AMDAinitSlaterDet(const Para* q, SlaterDet* Q);
void AMDAtoSlaterDet(const Para* q, SlaterDet* Q);
void AMDAprojectgradSlaterDet(const Para* q, 
		  	     const gradSlaterDet* dQ, double* dq); 


#endif

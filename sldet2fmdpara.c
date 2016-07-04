/**

  \file sldet2fmdpara.c

  convert SlaterDet into FMD parameterization


  (c) 2003 Thomas Neff

*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "fmd/SlaterDet.h"
#include "fmd/Parameterization.h"
#include "fmd/ParameterizationFMD.h"

#include "misc/utils.h"
#include "misc/physics.h"


#define MAXNG 3

void mirrorSlaterDet(SlaterDet* Q)
{
  int i;
  int Z=Q->Z; int N=Q->N;

  Q->Z=N; Q->N=Z;

  for (i=0; i<Q->ngauss; i++)
    Q->G[i].xi *= -1;
}


void doublegaussianSlaterDet(const SlaterDet* Q, SlaterDet* Qp, int trafo)
{
  int i,j;

  assert(2*Q->ngauss <= MAXNG*Q->A);

  Qp->A = Q->A; Qp->Z = Q->Z; Qp->N = Q->N;
  Qp->ngauss = 2*Q->ngauss;
  Qp->idx = (int*) malloc(Qp->A*sizeof(int));
  Qp->ng = (int*) malloc(Qp->A*sizeof(int));
  Qp->G = (Gaussian*) malloc(MAXNG*Qp->A*sizeof(Gaussian));

  for (i=0; i<Q->A; i++) {
    Qp->idx[i] = 2*Q->idx[i];
    Qp->ng[i] = 2*Q->ng[i];
  }
  for (i=0; i<2*Q->ngauss; i++)
    Qp->G[i] = Q->G[i/2];
  

  // polarization of gaussians
  if (trafo == 0) {
    for (i=1; i<2*Q->ngauss; i=i+2) {
      for (j=0; j<3; j++)
	Qp->G[i].b[j] += 0.5*((double) rand()/RAND_MAX);
      Qp->G[i].chi[0] += 0.1*((double) rand()/RAND_MAX);
      Qp->G[i].chi[1] *= 0.1*((double) rand()/RAND_MAX);  
    }
  }

  // haloize gaussians 
  if (trafo == 1) {
    for (i=1; i<2*Q->ngauss; i=i+2) {
      Qp->G[i].a += 2.0;
      Qp->G[i].chi[0] *= 0.1;
      Qp->G[i].chi[1] *= 0.1;
    }
  }
}


void doublegaussianSlaterDetParity(const SlaterDet* Q, SlaterDet* Qp, int par)
{
  int i;

  assert(2*Q->ngauss <= MAXNG*Q->A);

  Qp->A = Q->A; Qp->Z = Q->Z; Qp->N = Q->N;
  Qp->ngauss = 2*Q->ngauss;
  Qp->idx = (int*) malloc(Qp->A*sizeof(int));
  Qp->ng = (int*) malloc(Qp->A*sizeof(int));
  Qp->G = (Gaussian*) malloc(MAXNG*Qp->A*sizeof(Gaussian));

  for (i=0; i<Q->A; i++) {
    Qp->idx[i] = 2*Q->idx[i];
    Qp->ng[i] = 2*Q->ng[i];
  }
  for (i=0; i<Q->ngauss; i++) {
    Qp->G[0+2*i] = Q->G[i];
    Qp->G[1+2*i] = Q->G[i]; invertGaussian(&Qp->G[1+2*i]);
    if (par == -1) {
      Qp->G[1+2*i].chi[0] *= -1;
      Qp->G[1+2*i].chi[1] *= -1;
    }
  }
}
    


int main(int argc, char *argv[])
{
  createinfo(argc, argv);

  int mirror=0;
  int doublegauss=0;
  int doubleparity=0;
  int halo=0;
  int invert=0;
  double scale=0.0;
  int flip=-1;
  char c;
  
  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [OPTIONS] slaterdetfile fmdparafile\n"
	    "\n   -m         mirror nucleus"
	    "\n   -p         apply parity operator"
	    "\n   -d         double the number of gaussians"
            "\n   -P PAR     double the number of gaussians, second Gaussian is parity partner"
	    "\n   -h         double and haloize"
	    "\n   -s SCALE   scale dimensions in Slater determinant\n",
	    argv[0]);
    exit(-1);
  }
  
  while ((c = getopt(argc, argv, "mf:pdhs:P:")) != -1)
    switch (c) {
    case 'm':
      mirror=1;
      break;
    case 'f':
      flip = atoi(optarg);
      break;
    case 'p':
      invert=1;
      break;
    case 'd':
      doublegauss=1;
      break;
    case 'h':
      doublegauss=1;
      halo=1;
      break;
    case 'P':
      doubleparity=atoi(optarg);
      break;
    case 's':
      scale = atof(optarg);
      break;
    }

  char* slaterdetfile = argv[optind];
  char* fmdparafile = argv[optind+1];

  SlaterDet sldet;
  readSlaterDetfromFile(&sldet, slaterdetfile);

  if (mirror)
    mirrorSlaterDet(&sldet);

  if (flip != -1) {
    int i;
    for (i=0; i<sldet.ng[flip]; i++)
      spinflipGaussian(&sldet.G[sldet.idx[flip]+i]);
  }

  if (doublegauss) {
    SlaterDet sldetd;
    
    // squeeze by default
    doublegaussianSlaterDet(&sldet, &sldetd, halo);
    copySlaterDet(&sldetd, &sldet);
  }

  if (doubleparity) {
    SlaterDet sldetd;

    doublegaussianSlaterDetParity(&sldet, &sldetd, doubleparity);
    copySlaterDet(&sldetd, &sldet);
  }

  if (invert)
    invertSlaterDet(&sldet);

  if (scale != 0.0) {
    scaleSlaterDet(&sldet, scale);
  }

  Para fmdpara;
  SlaterDetinitFMD(&sldet, &fmdpara);

  // write FMD parameters to file
  
  FILE* parafp;

  
  if (!(parafp = fopen(fmdparafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", fmdparafile);
    exit(-1);
  }

  fprintinfo(parafp);

  fprintf(parafp, "\n<Parameterization FMD>\n");
  FMDwrite(parafp, &fmdpara);
  fprintf(parafp, "</Parameterization FMD>\n\n");

  writeSlaterDet(parafp, &sldet);

  fclose(parafp);

  return 0;
}

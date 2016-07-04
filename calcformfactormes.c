/**

  \file calcformfactormes.c

  calculate formfactor matrix elements

  DOES NOT WORK LIKE IT IS NOW, see calcformfactorsme.c


  (c) 2008 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/Projection.h"
#include "fmd/Symmetry.h"
#include "fmd/Formfactors.h"

#include "misc/utils.h"
#include "misc/physics.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/Projectionmpi.h"
#include "fmdmpi/FormfactorSlave.h"
#endif


#define MAXSTATES 100


void cleanup(int ret)
{
#ifdef MPI
  int task=TASKFIN;
  BroadcastTask(&task);

  MPI_Finalize();
#endif

  exit(ret);
}

int main(int argc, char* argv[])
{
  createinfo(argc, argv);

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    FormfactorSlave();

    MPI_Finalize();
  } else {
#endif

  // enough arguments ?

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s [OPTIONS] PROJPARA NUCSFILE\n"
	    "\n     -q QMAX         calculate to max momentum"
	    "\n     -n NPOINTS      number of points"
	    "\n     -p NALPHA-NBETA number of angles"
	    "\n     -l MULTIPOLE    select multipole transition"
	    "\n     -r              recoil\n",
	    argv[0]);
    cleanup(-1);
  }

  int npoints = 50+1;
  int nalpha = 5;
  int ncosb = 4;
  double qmax = 5.0;
  int recoil = 0;
  int l = 0;

  char* projparas;

  /* manage command-line options */

  char c;
  while ((c = getopt(argc, argv, "rq:n:p:l:")) != -1)
    switch (c) {
    case 'r':
      recoil = 1;
      break;
    case 'q':
      qmax = atof(optarg);
      break;
    case 'n':
      npoints = atoi(optarg);
      break;
    case 'p':
      projparas = optarg;
      char* ang = strtok(projparas, "-");
      nalpha = atoi(ang);
      ang = strtok(NULL, "-");
      ncosb = atoi(ang);
      break;
    case 'l':
      l = atoi(optarg);
      break;
    }

  char* projpar = argv[optind];
  char* nucsfile = argv[optind+1];

  char* mbfile[MAXSTATES];
  int n;

  if (readstringsfromfile(nucsfile, &n, mbfile))
    cleanup(-1);

  SlaterDet Q[n];
  Symmetry S[n];

  int i;
  for (i=0; i<n; i++) {
    extractSymmetryfromString(&mbfile[i], &S[i]);
    if (readSlaterDetfromFile(&Q[i], mbfile[i]))
      cleanup(-1);
  }


  // odd numer of nucleons ?
  int odd = Q[0].A % 2;

  // Projection parameters
  Projection P;
  initProjection(&P, odd, projpar);


  FormfactorPara FfP = {
    qmax : qmax,
    npoints : npoints,
    nalpha : nalpha,
    ncosb : ncosb, 
    recoil : recoil
  };


#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastParameters(&FfP, sizeof(FormfactorPara));
  BroadcastA(&Q[0].A);
#endif

  // calculate these multipoles
  int multi[NMULTIPOLES] = {0, 0, 0, 0};
  multi[l] = 1;

  // Project formfactor on angular momentum

  initOpFormfactors(&FfP);

  void* ffactorme[NMULTIPOLES];

  for (i=0; i<NMULTIPOLES; i++)   
    if (multi[i]) {
      ffactorme[i] = initprojectedMBME(&P, &OpMultipoleFormfactor[i]);


  int a,b;

  // read or calculate matrix elements

  for (b=0; b<n; b++)
    for (a=0; a<n; a++)
  
      if (readprojectedMBMEfromFile(mbfile[a], mbfile[b], &P,
				    &OpMultipoleFormfactor[i], S[a], S[b],
				    ffactorme[i])) {
#ifdef MPI				    
	calcprojectedMBMEmpi(&P, &OpMultipoleFormfactor[i], &Q[a], &Q[b],
			  S[a], S[b], ffactorme[i]);
#else
	calcprojectedMBME(&P, &OpMultipoleFormfactor[i], &Q[a], &Q[b],
			  S[a], S[b], ffactorme[i]);
#endif
	writeprojectedMBMEtoFile(mbfile[a], mbfile[b], &P, 
				 &OpMultipoleFormfactor[i], S[a], S[b], 
				 ffactorme[i]);
      } 

      else
	fprintf(stderr, "... matrix elements exist already\n");
    }

  cleanup(0);

#ifdef MPI
  }
#endif

}

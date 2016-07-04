/**

  \file minenergy.c

  minimize energy of Slater determinant


  (c) 2003, 2005 Thomas Neff

*/

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/Parameterization.h"
#include "fmd/ParameterizationFMD.h"
#include "fmd/SlaterDet.h"
#include "fmd/Interaction.h"
#include "fmd/Hamiltonian.h"
#include "fmd/gradSlaterDet.h"
#include "fmd/gradHamiltonian.h"
#include "fmd/gradOscillator.h"
#include "fmd/CenterofMass.h"
#include "fmd/Oscillator.h"
#include "fmd/Observables.h"

#include "MinimizerLBFGS.h"

#include "misc/physics.h"
#include "misc/utils.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/gradHamiltonianmpi.h"
#include "fmdmpi/MinimizerSlave.h"
#endif


volatile sig_atomic_t sigterminate = 0;

void catchterminate(int sig)
{
  sigterminate = 1;
  fprintf(stderr, "... catched signal - terminate minimization\n");
  signal(sig, catchterminate);
}


void cleanup(int ret)
{
  // stop all mpi clients

#ifdef MPI
  int task=TASKFIN;
  BroadcastTask(&task);

  MPI_Finalize();
#endif
 
  // terminate program with return code
  exit(ret);
}


int main(int argc, char *argv[])
{
  createinfo(argc, argv);

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  // fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    MinimizerSlave();

    MPI_Finalize();
  } else {
#endif

  int c;
  int cm=1;
  double omega=0.0;
  int overwrite=0;
  int maxsteps=100;
  double precision=0.0000001;
  double shakemag=0.0;

  // handler for INT and TERM signals
  signal(SIGINT, catchterminate);
  signal(SIGTERM, catchterminate);


  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction slaterdetfile\n"
	    "\n   -h OM           minimize in external oscillator"
	    "\n   -o              use new parameters even if energy has gone up"
	    "\n   -m MAXSTEPS     maximum number of steps (default %d)"
	    "\n   -p PRECISION    desired precision (default %f)"
	    "\n   -s MAGNITUDE    shake parameters before minimization\n",
	    argv[0], maxsteps, precision);
    cleanup(-1);
  }
  
  while((c = getopt(argc, argv, "h:om:p:s:")) != -1)
    switch (c) {
    case 'h':
      omega = atof(optarg)/hbc;
      break;
    case 'o':
      overwrite = 1;
      break;
    case 'm':
      maxsteps = atoi(optarg);
      break;
    case 'p':
      precision = atof(optarg);
      break;
    case 's':
      shakemag = atof(optarg);
      break;
    }
  
  char* interactionfile = argv[optind];
  char* parafile = argv[optind+1];

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    cleanup(-1);
  Int.cm = cm;

  Parameterization P;
  Para qinitial;
  if (readParafromFile(&P, &qinitial, parafile))
    cleanup(-1);

  SlaterDet Q;
  SlaterDetAux X;

  P.ParainitSlaterDet(&qinitial, &Q);
  initSlaterDetAux(&Q, &X);
  P.ParatoSlaterDet(&qinitial, &Q);

  calcSlaterDetAux(&Q, &X);
  double einitial;
  calcHamiltonian(&Int, &Q, &X, &einitial); 
  fprintf(stderr, "\ninitial:\tE = %8.3f MeV\n\n", hbc*einitial);

#ifdef MPI
  BroadcastInteraction(&Int);
  BroadcastA(&Q.A);
#endif

  Para q;
  P.Paraclone(&qinitial, &q);

  // shaking the parameters ?
  if (shakemag)
    shakePara(&q, shakemag);
  
  gradSlaterDetAux dX;
  gradSlaterDet dH;

  initgradSlaterDetAux(&Q, &dX);
  initgradSlaterDet(&Q, &dH);
  
  int steps = -1;
  double e, h;
  double dh[q.n];

  Minimizer mini;
  initMinimizer(&mini, q.n, precision);

  do {

    steps++;

    P.ParatoSlaterDet(&q, &Q);
    
    calcSlaterDetAux(&Q, &X);
    calcgradSlaterDetAux(&Q, &X, &dX);

#ifdef MPI
    calcgradHamiltonianmpi(&Int, &Q, &X, &dX, &dH);
#else
    calcgradHamiltonian(&Int, &Q, &X, &dX, &dH);
#endif
    e = dH.val;

    // in oscillator potential ?
    if (omega)
      calcgradOsci2(&Q, &X, &dX, omega, &dH);

    h = dH.val;
    P.ParaprojectgradSlaterDet(&q, &dH, dh);

    if (omega)
      fprintf(stderr, "step: %3d\tE = %8.3f MeV,\t Vosci = %8.3f MeV\n", 
	      steps, e*hbc, (h-e)*hbc);
    else
      fprintf(stderr, "step: %3d\tE = %8.3f MeV\n", 
	      steps, e*hbc);

  } while (!sigterminate && 
	   MinimizerStep(&mini, q.x, h, dh) && steps < maxsteps);

  // no improvement ? then exit
  if (!overwrite && einitial < e) 
    cleanup(0);

  P.ParatoSlaterDet(&q, &Q);

  // move SlaterDet to origin in phase space
  // in FMD parameterization we copy the relocated SlaterDet

  moveboostorientSlaterDet(&Q, &X);

  // normalize SlaterDet
  normalizeSlaterDet(&Q, &X);

  if (!strcmp(P.name, "FMD"))
    SlaterDetinitFMD(&Q, &q);

  // calculate the observables
  Observables Obs;
  
  calcObservables(&Int, &Q, &Obs);

  double vosci;
  if (omega) {
    calcSlaterDetAux(&Q, &X);
    calcOsci2(&Q, &X, omega, &vosci);
  }

  backup(parafile);

  fprintf(stderr, "... writing Parameters to file %s\n", parafile);

  FILE* outfp;
  if (!(outfp = fopen(parafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", parafile);
    cleanup(-1);
  }

  fprintinfo(outfp);

  fprintf(outfp, "\n# minimized %s for %s in %s parameterization\n"
	  "# using %s interaction\n", 
	  cm ? "< Hintr >" : "< H >", q.name, P.name, Int.name);

  if (sigterminate)
    fprintf(outfp, "\n# minimization was terminated prematurely\n");

  if (omega)
    fprintf(outfp, "\n# minimized in external oscillator,\tVosci = %8.3f MeV\n",
	    hbc*vosci);

  fprintObservables(outfp, &Int, &Q, &Obs);

  fprintf(outfp, "\n# Parameterization\n");
  fprintf(outfp, "<Parameterization %s>\n", P.name);
  P.Parawrite(outfp, &q);
 
  fprintf(outfp, "\n# SlaterDet\n");
  writeSlaterDet(outfp, &Q);

  fclose(outfp);

  fprintf(stderr, "... %4.2f minutes computing time used\n", usertime()/60.0);

  cleanup(0);

#ifdef MPI
  }
#endif
}

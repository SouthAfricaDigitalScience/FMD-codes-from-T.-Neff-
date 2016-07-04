/**

  \file minenergyp.c

  minimize energy of parity projected Slater determinant


  (c) 2003,2004 Thomas Neff

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


int main(int argc, char *argv[])
{
  createinfo(argc, argv);

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

  fprintf(stderr, "... [%2d] %s\n", mpirank, hostname());

  if (mpirank != 0) {
    MinimizerSlaveod();
    MPI_Finalize();
    exit(0);
  } else {
#endif

  int c;
  int cm=0;
  int par=0;
  double omega=0.0;
  int overwrite=0;
  int maxsteps=250;
  double precision=0.00001;
  double shakemag=0.0;

  // handler for INT and TERM signals
  signal(SIGINT, catchterminate);
  signal(SIGTERM, catchterminate);


  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction slaterdetfile\n"
	    "\n   -c              minimize H-Tcm"
	    "\n   -p PARITY       minimize parity projected energy"
	    "\n   -h OM           minimize in external oscillator"
	    "\n   -o              use new parameters even if energy has gone up"
	    "\n   -m MAXSTEPS     maximum number of steps (default %d)"
	    "\n   -a PRECISION    desired precision (default %f)"
	    "\n   -s MAGNITUDE    shake parameters before minimization\n",
	    argv[0], maxsteps, precision);
    exit(-1);
  }
  
  while((c = getopt(argc, argv, "a:ch:om:p:s:")) != -1)
    switch (c) {
    case 'c':
      cm = 1;
      break;
    case 'p':
      par = atoi(optarg);
      break;
    case 'h':
      omega = atof(optarg)/hbc;
      break;
    case 'o':
      overwrite = 1;
      break;
    case 'm':
      maxsteps = atoi(optarg);
      break;
    case 'a':
      precision = atof(optarg);
      break;
    case 's':
      shakemag = atof(optarg);
      break;
    }
  
  char* interactionfile = argv[optind];
  char* parafile = argv[optind+1];

  Interaction Int;
  readInteractionfromFile(&Int, interactionfile);
  Int.cm = cm;

#ifdef MPI
  BroadcastInteraction(&Int);
#endif

  Parameterization P;
  Para qinitial;
  readParafromFile(&P, &qinitial, parafile); 

  SlaterDet Q, Qp;
  SlaterDetAux X;

  P.ParainitSlaterDet(&qinitial, &Q);
  P.ParainitSlaterDet(&qinitial, &Qp);
  initSlaterDetAux(&Q, &X);
  P.ParatoSlaterDet(&qinitial, &Q);

#ifdef MPI
  BroadcastA(&Q.A);
#endif

  complex double hd, hp, nd, np;
  double einitial, eintr;

  if (!par) {
    par = (Q.A % 2) ? -1 : 1;
  }
   
  calcSlaterDetAuxod(&Q, &Q, &X);
  nd = X.ovlap;
  calcHamiltonianod(&Int, &Q, &Q, &X, &hd);
  eintr = hd/nd;

  copySlaterDet(&Q, &Qp);
  invertcmSlaterDet(&Qp, &X);
  calcSlaterDetAuxod(&Q, &Qp, &X);
  np = X.ovlap;
  calcHamiltonianod(&Int, &Q, &Qp, &X, &hp);
  
  einitial = (hd+par*hp)/(nd+par*np);
  
  fprintf(stderr, "\ninitial:\tE = %8.3f MeV,   Eintr = %8.3f MeV\n\n", 
	  hbc*einitial, hbc*eintr);

  Para q;
  P.Paraclone(&qinitial, &q);

  // shaking the parameters ?
  if (shakemag)
    shakePara(&q, shakemag);
  
  gradSlaterDetAux dX;
  gradSlaterDet dH, dhd, dhp, dnd, dnp;
  double H, N;
    
  initgradSlaterDetAux(&Q, &dX);
  initgradSlaterDet(&Q, &dhd);
  initgradSlaterDet(&Q, &dnd);
  initgradSlaterDet(&Q, &dhp);
  initgradSlaterDet(&Q, &dnp);
  initgradSlaterDet(&Q, &dH);
  
  int steps = -1;
  double hn;
  double dhn[q.n];

  Minimizer mini;
  initMinimizer(&mini, q.n, precision);

  do {

    steps++;

    P.ParatoSlaterDet(&q, &Q);

    calcSlaterDetAuxod(&Q, &Q, &X);
    calcgradSlaterDetAuxod(&Q, &Q, &X, &dX);
    calcgradOvlapod(&Q, &Q, &X, &dX, &dnd);
#ifdef MPI
    calcgradHamiltonianodmpi(&Int, &Q, &Q, &X, &dX, &dhd);
#else
    calcgradHamiltonianod(&Int, &Q, &Q, &X, &dX, &dhd);
#endif
    hd = dhd.val;
    nd = X.ovlap;
    
    copySlaterDet(&Q, &Qp);
    invertcmSlaterDet(&Qp, &X);
    calcSlaterDetAuxod(&Q, &Qp, &X);
    calcgradSlaterDetAuxod(&Q, &Qp, &X, &dX);
    calcgradOvlapod(&Q, &Qp, &X, &dX, &dnp);
#ifdef MPI
    calcgradHamiltonianodmpi(&Int, &Q, &Qp, &X, &dX, &dhp);
#else
    calcgradHamiltonianod(&Int, &Q, &Qp, &X, &dX, &dhp);
#endif
    hp = dhp.val;
    np = X.ovlap;
    
    H = hd+par*hp;
    N = nd+par*np;
    
    zerogradSlaterDet(&dH);
    addmulttogradSlaterDet(&dH, &dhd, 1.0/N);
    addmulttogradSlaterDet(&dH, &dhp, par* 1.0/N);
    addmulttogradSlaterDet(&dH, &dnd, -H/(N*N));
    addmulttogradSlaterDet(&dH, &dnp, -par*H/(N*N));
    
    hn = H/N;

    // in oscillator potential ?
    // if (omega)
    //  calcgradOsci2(&Q, &X, &dX, omega, &dH);

    P.ParaprojectgradSlaterDet(&q, &dH, dhn);

    // if (omega)
    //  fprintf(stderr, "step: %3d\tE = %8.3f MeV,\t Vosci = %8.3f MeV\n", 
    //          steps, e*hbc, (h-e)*hbc);
    // else
    fprintf(stderr, "step: %3d\tE = %8.3f MeV,   Eintr = %8.3f MeV\n", 
	    steps, hbc*hn, hbc*creal(hd/nd));

  } while (!sigterminate && 
	   MinimizerStep(&mini, q.x, hn, dhn) && steps < maxsteps);

  // no improvement ? then write initial parameters
  if (!overwrite && einitial < hn) 
    P.Paraclone(&qinitial, &q);

  P.ParatoSlaterDet(&q, &Q);

  // orient along the main axes
  // in FMD parameterization we copy the relocated SlaterDet

  moveboostorientSlaterDet(&Q, &X);

  if (!strcmp(P.name, "FMD"))
    SlaterDetinitFMD(&Q, &q);

  // calculate the observables
  Observables Obs;
  
  calcObservablesParity(&Int, par, &Q, &Obs);

  // double vosci;
  // if (omega)
  //  calcOsci2(&Q, &X, omega, &vosci);

  backup(parafile);

  fprintf(stderr, "... writing Parameters to file %s\n", parafile);

  FILE* outfp;
  if (!(outfp = fopen(parafile, "w"))) {
    fprintf(stderr, "couldn't open %s for writing\n", parafile);
    exit(-1);
  }

  fprintinfo(outfp);

  fprintf(outfp, "\n# minimized %s for %s in %s parameterization\n"
	  "# parity projected to %+d\n"
	  "# using %s interaction\n", 
	  cm ? "< Hintr >" : "< H >", q.name, P.name, par, Int.name);

  if (sigterminate)
    fprintf(outfp, "\n# minimization was terminated prematurely\n");

  // if (omega)
  //  fprintf(outfp, "\n# minimized in external oscillator,\tVosci = %8.3f MeV\n",
  //	    hbc*vosci);

  fprintObservables(outfp, &Int, &Q, &Obs);

  fprintf(outfp, "\n# Parameterization\n");
  fprintf(outfp, "<Parameterization %s>\n", P.name);
  P.Parawrite(outfp, &q);
 
  fprintf(outfp, "\n# SlaterDet\n");
  writeSlaterDet(outfp, &Q);

  fclose(outfp);

  fprintf(stderr, "... %4.2f minutes computing time used\n", usertime()/60.0);

#ifdef MPI
  int task=TASKFIN;
  BroadcastTask(&task);
  }

  MPI_Finalize();
#endif

  return 0;
}

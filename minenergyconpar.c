/**

  \file minenergyconpar.c

  minimize energy of constrained and parity projected Slater determinant
  constraints are evaluated for parity projected state


  (c) 2003,2004 Thomas Neff

*/

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
#include "fmd/Oscillator.h"
#include "fmd/gradOscillator.h"

#include "fmd/Constraint.h"
#include "fmd/ConstraintCM.h"
// #include "fmd/ConstraintCMXP.h"
#include "fmd/ConstraintT2.h"
#include "fmd/ConstraintS2.h"
#include "fmd/Constraintl2.h"
// #include "fmd/Constraintls.h"
#include "fmd/ConstraintJ2.h"
#include "fmd/ConstraintLS.h"
#include "fmd/ConstraintR2.h"
#include "fmd/ConstraintDipole.h"
#include "fmd/ConstraintQuadrupole.h"
#include "fmd/ConstraintOctupole.h"

#include "fmd/CenterofMass.h"
#include "fmd/SpatialOrientation.h"
#include "fmd/Observables.h"

#include "MinimizerDONLP2par.h" 

#include "misc/physics.h"
#include "misc/utils.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/MinimizerSlave.h"
#endif


#define SQR(x) (x)*(x)

#define MAXCONSTRAINT 10

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
    MinimizerSlaveod();
    MPI_Finalize();
    exit(0);
  } else {
#endif

  int c;
  int cm=1;
  int par=0;
  int overwrite=0;
  int maxsteps=250;
  double shakemag=0.0;
  int constT2=0, constS2=0, constJ2=0, constLS=0;
  int constmradius=0, constmquadrupole=0;
  int consteradius=0, constequadrupole=0;
  double constT2val=0.0, constS2val=0.0, constJ2val=0.0, constLSval=0.0;
  double constmrval=0.0, constmqval=0.0;
  double consterval=0.0, consteqval=0.0;

  /*
  int constl2=0, constls=0, constj2=0;
  int constmradius=0, constmquadrupole=0, constmoctupole=0;
  int consteradius=0, constedipole=0, constequadrupole=0, consteoctupole=0;
  int constbeta=0, constgamma=0;
  double constt2val=0.0, consts2val=0.0;
  double constl2val=0.0, constlsval=0.0, constj2val=0.0;
  double constmrval=0.0, constmqval=0.0, constmoval=0.0;
  double consterval=0.0, constedval=0.0, consteqval=0.0, consteoval=0.0;
  double constbval=0.0, constgval=0.0;
  */

  // first constraint is ConstraintCM
  // special handling as it is of type Constraint and not Constraintod
  Constraintod Const[MAXCONSTRAINT];
  int nconst = 1;

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction slaterdetfile\n"
	    "\n   -p PARITY       project parity"
	    "\n   -o              overwrite in all cases"
	    "\n   -m MAXSTEPS     maximum number of steps (default %d)"
	    "\n   -s MAGNITUDE    shake parameters before minimization"
	    "\n   -C              constrain center of mass"
	    "\n   -T T2           constrain isospin"
	    "\n   -S S2           constrain spin"
	    "\n   -L LS           constrain single-particle l2"
	    "\n   -J J2           constrain angular momentum"
	    "\n   -R RADIUS       constrain radius"
	    "\n   -D DIPOLE       constrain dipole moments"
	    "\n   -Q QUADRUPOLE   constrain quadrupole moments"
	    "\n   -O OCTUPOLE     constrain octupole moments"
	    "\n   -B BETA         constrain quadrupole deformation"
	    "\n   -G GAMMA        constrain quadrupole deformation\n"
	    // "\n   -M              constrain main axes of quadrupole tensor"
	    // "\n   -I Jz           cranking constrain"
	    , argv[0], maxsteps);
    cleanup(-1);
  }
  
  while((c = getopt(argc, argv, "p:om:s:CT:L:S:J:R:D:Q:O:B:G:")) != -1)
    switch (c) {
    case 'p':
      par = atoi(optarg);
      break;
    case 'o':
      overwrite = 1;
      break;
    case 'm':
      maxsteps = atoi(optarg);
      break;
    case 's':
      shakemag = atof(optarg);
      break;
    case 'T':
      constT2 = 1;
      constT2val = atof(optarg);
      break;
    case 'L':
      constLS = 1;
      constLSval = atof(optarg);
      break;
    case 'S':
      constS2 = 1;
      constS2val = atof(optarg);
      break;
    case 'J':
      constJ2 = 1;
      constJ2val = atof(optarg);
      break;
    case 'R':
      if (optarg[0] == 'E') {
	consteradius = 1;
	consterval = atof(optarg+2);
      } else {
	constmradius = 1;
	constmrval = atof(optarg);
      }
      break;
      /*
    case 'D':
      constedipole = 1;
      constedval = atof(optarg);
      break;
      */
    case 'Q':
      if (optarg[0] == 'E') {
	constequadrupole = 1;
	consteqval = atof(optarg+2);
      } else {
	constmquadrupole = 1;
	constmqval = atof(optarg);
      }
      break;
      /*
    case 'O':
      if (optarg[0] == 'E') {
	consteoctupole = 1;
	consteoval = atof(optarg+2);
      } else {
	constmoctupole = 1;
	constmoval = atof(optarg);
      }
      break;
    case 'B':
      constbeta = 1;
      constbval = atof(optarg);
      break;
    case 'G':
      constgamma = 1;
      constgval = atof(optarg);
      break;
      */
    }
  
  if (constT2) {
    Const[nconst] = ConstraintT2od;
    Const[nconst].val = constT2val;
    nconst++;
  }
  if (constS2) {
    Const[nconst] = ConstraintS2od;
    Const[nconst].val = constS2val;
    nconst++;
  }
  /*
  if (constl2) {
    Const[nconst] = Constraintl2;
    Const[nconst].val = constl2val;
    nconst++;
  }
  if (constls) {
    Const[nconst] = Constraintls;
    Const[nconst].val = constlsval;
    nconst++;
  }
  */
  if (constJ2) {
    Const[nconst] = ConstraintJ2od;
    Const[nconst].val = constJ2val;
    nconst++;
  }
  if (constLS) {
    Const[nconst] = ConstraintLSod;
    Const[nconst].val = constLSval;
    nconst++;
  }
  if (constmradius) {
    Const[nconst] = ConstraintR2od;
    Const[nconst].val = SQR(constmrval);
    nconst++;
  }
  if (consteradius) {
    Const[nconst] = ConstraintER2od;
    Const[nconst].val = SQR(consterval);
    nconst++;
  }
  /*
  if (constedipole) {
    Const[nconst] = ConstraintED2;
    Const[nconst].val = constedval;
    nconst++;
  }
  */
  if (constmquadrupole) {
    Const[nconst] = ConstraintQ2od;
    Const[nconst].val = constmqval;
    nconst++;
  }
  if (constequadrupole) {
    Const[nconst] = ConstraintEQ2od;
    Const[nconst].val = consteqval;
    nconst++;
  }
  /*
  if (constmoctupole) {
    Const[nconst] = ConstraintO2;
    Const[nconst].val = constmoval;
    nconst++;
  }
  if (consteoctupole) {
    Const[nconst] = ConstraintEO2;
    Const[nconst].val = consteoval;
    nconst++;
  }
  if (constbeta) {
    Const[nconst] = ConstraintQ2;
    Const[nconst].val = sqrt(24*M_PI/5)*constbval;
    nconst++;
  }
  if (constgamma) {
    if (!constbeta) {
      fprintf(stderr, "gamma constraint only in combination with beta constraint\n");
      exit(-1);
    }
    Const[nconst] = ConstraintDetQ;
    Const[nconst].val = sqrt(4*M_PI/5)*constbval*cbrt(2*cos(3*constgval*M_PI/180));
    nconst++;
  }
  */

  char* interactionfile = argv[optind];
  char* parafile = argv[optind+1];

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    cleanup(-1);
  Int.cm = cm;

#ifdef MPI
  BroadcastInteraction(&Int);
#endif

  Parameterization P;
  Para qinitial;
  if (readParafromFile(&P, &qinitial, parafile))
    cleanup(-1);

  SlaterDet Q, Qp;
  SlaterDetAux X;

  P.ParainitSlaterDet(&qinitial, &Q);
  P.ParainitSlaterDet(&qinitial, &Qp);
  if (!par)
    par = (Q.A % 2) ? -1 : +1;

  initSlaterDetAux(&Q, &X);
  P.ParatoSlaterDet(&qinitial, &Q);
  P.ParatoSlaterDet(&qinitial, &Q);

#ifdef MPI
  BroadcastA(&Q.A);
#endif

  complex double hd, hp, nd, np;
  double einitial, eintr;

  moveboostorientSlaterDet(&Q, &X);

  calcSlaterDetAuxod(&Q, &Q, &X);
  nd = X.ovlap;
  calcHamiltonianod(&Int, &Q, &Q, &X, &hd);
  eintr = hd/nd;

  copySlaterDet(&Q, &Qp);
  invertSlaterDet(&Qp);
  calcSlaterDetAuxod(&Q, &Qp, &X);
  np = X.ovlap;
  calcHamiltonianod(&Int, &Q, &Qp, &X, &hp);
  
  einitial = (hd+par*hp)/(nd+par*np);

  int i;
  for (i=1; i<nconst; i++)
    fprintf(stderr, "# constraining %4s to %8.3f\n", 
	    Const[i].label, Const[i].output(Const[i].val)); 

  fprintf(stderr, "\ninitial:\tE = %8.3f MeV,   Eintr = %8.3f MeV\n\n", 
	  hbc*einitial, hbc*eintr);

  Para q;
  P.Paraclone(&qinitial, &q);

  // shaking the parameters ?
  if (shakemag) {
    shakePara(&q, shakemag);
  }

  // in FMD parameterization move to origin in phase space
  // if (!strcmp(P.name, "FMD")) {
  //   P.ParatoSlaterDet(&q, &Q);
  //   moveboostorientSlaterDet(&Q, &X);
  //   SlaterDetinitFMD(&Q, &q);
  // }

  // minimize !
  MinimizeDONLP2par(&Int, par, Const, nconst, &P, &q, maxsteps); 

  double e;
  P.ParatoSlaterDet(&q, &Q);

  calcSlaterDetAuxod(&Q, &Q, &X);
  nd = X.ovlap;
  calcHamiltonianod(&Int, &Q, &Q, &X, &hd);
  eintr = hd/nd;

  copySlaterDet(&Q, &Qp);
  invertSlaterDet(&Qp);
  calcSlaterDetAuxod(&Q, &Qp, &X);
  np = X.ovlap;
  calcHamiltonianod(&Int, &Q, &Qp, &X, &hp);
  
  e = (hd+par*hp)/(nd+par*np);

  fprintf(stderr, "\nfinal:  \tE = %8.3f MeV\n\n", hbc*e);

  // output files

  // backup of parafile

  backup(parafile);

  // save minimization result in MIN directory
  {
    P.ParatoSlaterDet(&q, &Q);

    // orient SlaterDet
    moveboostorientSlaterDet(&Q, &X);

    // in FMD parameterization we copy the relocated SlaterDet
    if (!strcmp(P.name, "FMD"))
      SlaterDetinitFMD(&Q, &q);

    // calculate the observables
    Observables Obs;
  
    calcObservablesParity(&Int, par, &Q, &Obs);

    ensuredir("MIN");
    char outfile[1024];
    sprintf(outfile, "MIN/%s-minpar.%d", parafile, (int) time(NULL));

    fprintf(stderr, "... writing Parameters to file %s\n", outfile);

    FILE* outfp;
    if (!(outfp = fopen(outfile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", outfile);
      cleanup(-1);
    }

    fprintinfo(outfp);

    fprintf(outfp, "\n# minimized %s for %s in %s parameterization\n"
            "# parity projected to %+d\n"
            "# using %s interaction\n", 
            cm ? "< Hintr >" : "< H >", q.name, P.name, par, Int.name);
  
    fprintf(outfp, "# einitial: %8.3f MeV\n", hbc*einitial);
    fprintf(outfp, "# efinal:   %8.3f MeV\n", hbc*e);

    fprintObservables(outfp, &Int, &Q, &Obs);

    fprintf(outfp, "\n# Parameterization\n");
    fprintf(outfp, "<Parameterization %s>\n", P.name);
    P.Parawrite(outfp, &q);
 
    fprintf(outfp, "\n# SlaterDet\n");
    writeSlaterDet(outfp, &Q);

    fclose(outfp);  
  }

  // save minimization result in parafile
  {
    // no improvement ? then write initial parameters
    if (!overwrite && einitial < e) 
      P.Paraclone(&qinitial, &q);

    P.ParatoSlaterDet(&q, &Q);

    // orient SlaterDet
    moveboostorientSlaterDet(&Q, &X);

    // in FMD parameterization we copy the relocated SlaterDet
    if (!strcmp(P.name, "FMD"))
      SlaterDetinitFMD(&Q, &q);

    // calculate the observables
    Observables Obs;
  
    calcObservablesParity(&Int, par, &Q, &Obs);

    fprintf(stderr, "... writing Parameters to file %s\n", parafile);

    FILE* outfp;
    if (!(outfp = fopen(parafile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", parafile);
      cleanup(-1);
    }

    fprintinfo(outfp);

    fprintf(outfp, "\n# minimized %s for %s in %s parameterization\n"
            "# parity projected to %+d\n"
            "# using %s interaction\n", 
            cm ? "< Hintr >" : "< H >", q.name, P.name, par, Int.name);
  
    fprintf(outfp, "\n# einitial: %8.3f MeV\n", hbc*einitial);
    fprintf(outfp, "# efinal:   %8.3f MeV\n", hbc*e);

    if (overwrite)
      fprintf(outfp, "\n# overwrite flag: use new parameters\n\n");
    if (!overwrite && einitial < e)
      fprintf(outfp, "\n# no improvement by minimization: use initial parameters\n\n");

    fprintObservables(outfp, &Int, &Q, &Obs);

    fprintf(outfp, "\n# Parameterization\n");
    fprintf(outfp, "<Parameterization %s>\n", P.name);
    P.Parawrite(outfp, &q);
 
    fprintf(outfp, "\n# SlaterDet\n");
    writeSlaterDet(outfp, &Q);

    fclose(outfp);  
  }

  fprintf(stderr, "... %4.2f minutes computing time used\n", usertime()/60.0);

  cleanup(0);

#ifdef MPI
  }
#endif
}

/**

  \file minenergyconp.c

  minimize energy of constrained and parity projected Slater determinant


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
// #include "fmd/Constraintl2.h"
// #include "fmd/Constraintls.h"
#include "fmd/ConstraintJ2.h"
#include "fmd/ConstraintJ2.h"
#include "fmd/ConstraintR2.h"
#include "fmd/ConstraintNOsci.h"
#include "fmd/ConstraintDipole.h"
#include "fmd/ConstraintQuadrupole.h"
#include "fmd/ConstraintOctupole.h"

#include "fmd/ConstraintSymmetries.h"

#include "fmd/CenterofMass.h"
#include "fmd/SpatialOrientation.h"
#include "fmd/Observables.h"

#include "MinimizerDONLP2p.h" 

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
  int log=0;
  int maxsteps=250;
  double shakemag=0.0;
  int constcm=1, constT2=0, constL2=0, constLS=0, constJ2=0;
  int constS2=0, constpS2=0, constnS2=0;
  int constj2=0;
  int constnosci=0, constpnosci=0, constnnosci=0;
  int constmradius=0, constmquadrupole=0, constmoctupole=0;
  int consteradius=0, constedipole=0, constequadrupole=0, consteoctupole=0;
  int constnradius=0, constnquadrupole=0, constnoctupole=0;
  int constbeta=0, constgamma=0;
  int constpbeta=0, constpgamma=0, constnbeta=0, constngamma=0;
  int constparp=0, constparn=0, constparch=0, constd3h=0, constsxz=0, constsyz=0;

  double constT2val=0.0;
  double constS2val=0.0, constpS2val=0.0, constnS2val=0.0;
  double constL2val=0.0, constLSval=0.0, constJ2val=0.0;
  double constj2val=0.0;
  double constnoscival=0.0, constpnoscival=0.0, constnnoscival=0.0;
  double constmrval=0.0, constmqval=0.0, constmoval=0.0;
  double consterval=0.0, constedval=0.0, consteqval=0.0, consteoval=0.0;
  double constnrval=0.0, constnqval=0.0, constnoval=0.0;
  double constbval=0.0, constgval=0.0;
  double constpbval=0.0, constpgval=0.0, constnbval=0.0, constngval=0.0;

  Constraint Const[MAXCONSTRAINT];
  int nconst = 0;

  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction slaterdetfile\n"
	    "\n   -p PARITY       project parity"
	    "\n   -o              overwrite in all cases"
	    "\n   -l EVERY        log EVERY step"
	    "\n   -m MAXSTEPS     maximum number of steps (default %d)"
	    "\n   -s MAGNITUDE    shake parameters before minimization"
	    "\n   -C              constrain center of mass"
	    "\n   -S S2           constrain spin"
	    "\n   -T T2           constrain isospin"
	    "\n   -L l2           constrain single-particle l2"
	    "\n   -J J2           constrain angular momentum"
	    "\n   -j j2           constrain single-particle angular momentum"
	    "\n   -N NOSCI        constrain number of oscillator quanta"
	    "\n   -R [P|N:]RADIUS       constrain radius"
	    "\n   -D DIPOLE       constrain dipole moments"
	    "\n   -Q [P|N:]QUADRUPOLE   constrain quadrupole moments"
	    "\n   -O [P|N:]OCTUPOLE     constrain octupole moments"
	    "\n   -B [P|N:]BETA         constrain quadrupole deformation"
	    "\n   -G [P|N:]GAMMA        constrain quadrupole deformation"
            "\n   -P[+1|-1]       constrain parity of intrinsic state"
            "\n   -C              constrain on symmetry under parity and charge"
            "\n   -V [xz|yz]      constrain reflection symmetry of intrinsic state"
            "\n   -3              constrain intrinsic state to d3h symmetry\n"
	    // "\n   -M              constrain main axes of quadrupole tensor"
	    // "\n   -I Jz           cranking constrain"
	    , argv[0], maxsteps);
    cleanup(-1);
  }
  
  while((c = getopt(argc, argv, "p:ol:m:s:T:L:S:J:j:N:R:D:Q:O:B:G:P:CV:3")) != -1)
    switch (c) {
    case 'p':
      par = atoi(optarg);
      break;
    case 'o':
      overwrite = 1;
      break;
    case 'l':
      log = atoi(optarg);
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
      constL2 = 1;
      constL2val = atof(optarg);
      break;
    case 'S':
      if (optarg[0] == 'P') {
        constpS2 = 1;
        constpS2val = atof(optarg+2);
      } else if (optarg[0] == 'N') {
        constnS2 = 1;
        constnS2val = atof(optarg+2);
      } else {
        constS2 = 1;
        constS2val = atof(optarg);
      }
      break;
    case 'J':
      constJ2 = 1;
      constJ2val = atof(optarg);
      break;
    case 'j':
      constj2 = 1;
      constj2val = atof(optarg);
      break;
    case 'N':
      if (optarg[0] == 'P') {
        constpnosci = 1;
        constpnoscival = atof(optarg+2);
      } else if (optarg[0] == 'N') {
        constnnosci = 1;
        constnnoscival = atof(optarg+2);
      } else {
        constnosci = 1;
        constnoscival = atof(optarg);
      }
      break;
    case 'R':
      if (optarg[0] == 'E' || optarg[0] == 'P') {
	consteradius = 1;
	consterval = atof(optarg+2);
      } else if (optarg[0] == 'N') {
	constnradius = 1;
	constnrval = atof(optarg+2);
      } else {
	constmradius = 1;
	constmrval = atof(optarg);
      }
      break;
    case 'D':
      constedipole = 1;
      constedval = atof(optarg);
      break;
    case 'Q':
      if (optarg[0] == 'E' || optarg[0] == 'P') {
	constequadrupole = 1;
	consteqval = atof(optarg+2);
      } else if (optarg[0] == 'N') {
	constnquadrupole = 1;
	constnqval = atof(optarg+2);
      } else {
	constmquadrupole = 1;
	constmqval = atof(optarg);
      }
      break;
    case 'O':
      if (optarg[0] == 'E') {
	consteoctupole = 1;
	consteoval = atof(optarg+2);
      } else if (optarg[0] == 'N') {
	constnoctupole = 1;
	constnoval = atof(optarg+2);
      } else {
	constmoctupole = 1;
	constmoval = atof(optarg);
      }
      break;
    case 'B':
      if (optarg[0] == 'P') {
	constpbeta = 1;
	constpbval = atof(optarg+2);
      } else if (optarg[0] == 'N') {
	constnbeta = 1;
	constnbval = atof(optarg+2);
      } else {
	constbeta = 1;
	constbval = atof(optarg);
      }
      break;
    case 'G':
      if (optarg[0] == 'P') {
	constpgamma = 1;
	constpgval = atof(optarg+2);
      } else if (optarg[0] == 'N') {
	constngamma = 1;
	constngval = atof(optarg+2);
      } else {
	constgamma = 1;
	constgval = atof(optarg);
      }
      break;
    case 'P':
      if (atoi(optarg) == +1) 
        constparp = 1;
      else
        constparn = 1;
      break;
    case 'C':
      constparch = 1;
      break;
    case 'V':
      if (!strcmp(optarg, "xz"))
        constsxz = 1;
      else if (!strcmp(optarg, "yz"))
        constsyz = 1;
      break;
    case '3':
      constd3h = 1;
      break;
    }
  
  if (constcm) {
    Const[nconst] = ConstraintCM;
    nconst++;
    // Const[nconst  ] = ConstraintX;
    // Const[nconst+1] = ConstraintY;
    // Const[nconst+2] = ConstraintZ;
    // Const[nconst+3] = ConstraintPX;
    // Const[nconst+4] = ConstraintPY;
    // Const[nconst+5] = ConstraintPZ;
    // nconst += 6;
  }
  if (constT2) {
    Const[nconst] = ConstraintT2;
    Const[nconst].val = constT2val;
    nconst++;
  }
  if (constS2) {
    Const[nconst] = ConstraintS2;
    Const[nconst].val = constS2val;
    nconst++;
  }
  if (constpS2) {
    Const[nconst] = ConstraintPS2;
    Const[nconst].val = constpS2val;
    nconst++;
  }
  if (constnS2) {
    Const[nconst] = ConstraintNS2;
    Const[nconst].val = constnS2val;
    nconst++;
  }
/*   if (constL2) { */
/*     Const[nconst] = ConstraintL2; */
/*     Const[nconst].val = constL2val; */
/*     nconst++; */
/*   } */
/*   if (constLS) { */
/*     Const[nconst] = ConstraintLS; */
/*     Const[nconst].val = constLSval; */
/*     nconst++; */
/*   } */
  if (constJ2) {
    Const[nconst] = ConstraintJ2;
    Const[nconst].val = constJ2val;
    nconst++;
  }
  if (constj2) {
    Const[nconst] = ConstraintJ2;
    Const[nconst].val = constj2val;
    nconst++;
  }
  if (constnosci) {
    Const[nconst] = ConstraintNOsci;
    Const[nconst].val = SQR(constnoscival);
    nconst++;
  }
  if (constpnosci) {
    Const[nconst] = ConstraintPNOsci;
    Const[nconst].val = SQR(constpnoscival);
    nconst++;
  }
  if (constnnosci) {
    Const[nconst] = ConstraintNNOsci;
    Const[nconst].val = SQR(constnnoscival);
    nconst++;
  }
  if (constmradius) {
    Const[nconst] = ConstraintR2;
    Const[nconst].val = SQR(constmrval);
    nconst++;
  }
  if (consteradius) {
    Const[nconst] = ConstraintER2;
    Const[nconst].val = SQR(consterval);
    nconst++;
  }
  if (constnradius) {
    Const[nconst] = ConstraintNR2;
    Const[nconst].val = SQR(constnrval);
    nconst++;
  }
  if (constedipole) {
    Const[nconst] = ConstraintED2;
    Const[nconst].val = constedval;
    nconst++;
  }
  if (constmquadrupole) {
    Const[nconst] = ConstraintQ2;
    Const[nconst].val = constmqval;
    nconst++;
  }
  if (constequadrupole) {
    Const[nconst] = ConstraintEQ2;
    Const[nconst].val = consteqval;
    nconst++;
  }
  if (constnquadrupole) {
    Const[nconst] = ConstraintNQ2;
    Const[nconst].val = constnqval;
    nconst++;
  }
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
  if (constnoctupole) {
    Const[nconst] = ConstraintNO2;
    Const[nconst].val = constnoval;
    nconst++;
  }
  if (constbeta) {
    Const[nconst] = ConstraintQ2;
    Const[nconst].val = sqrt(24*M_PI/5)*constbval;
    nconst++;
  }
  if (constpbeta) {
    Const[nconst] = ConstraintEQ2;
    Const[nconst].val = sqrt(24*M_PI/5)*constpbval;
    nconst++;
  }
  if (constnbeta) {
    Const[nconst] = ConstraintNQ2;
    Const[nconst].val = sqrt(24*M_PI/5)*constnbval;
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
  if (constpgamma || constngamma) {
    fprintf(stderr, "proton/neutron gamma constraint not implemented yet !\n");
    exit(-1);
  }
  if (constparp) {
    Const[nconst] = ConstraintParityP;
    nconst++;
  }
  if (constparn) {
    Const[nconst] = ConstraintParityN;
    nconst++;
  }
  if (constparch) {
    Const[nconst] = ConstraintParityCharge;
    nconst++;
  }
  if (constsxz) {
    Const[nconst] = ConstraintSxz;
    nconst++;
  }
  if (constsyz) {
    Const[nconst] = ConstraintSyz;
    nconst++;
  }
  if (constd3h) {
    Const[nconst] = ConstraintD3H;
    nconst++;
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
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastInteraction(&Int);
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
  for (i=0; i<nconst; i++)
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

  // in FMD parameterization move SlaterDet to origin in phase space
  // if (!strcmp(P.name, "FMD")) {
  //   P.ParatoSlaterDet(&q, &Q);
  //   moveboostorientSlaterDet(&Q, &X);
  //   SlaterDetinitFMD(&Q, &q);
  // }

  // minimize !
  MinimizeDONLP2p(&Int, par, Const, nconst, &P, &q, maxsteps, log, parafile); 

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

    // normalize SlaterDet
    normalizeSlaterDet(&Q, &X);

    // in FMD parameterization we copy the relocated SlaterDet
    if (!strcmp(P.name, "FMD"))
      SlaterDetinitFMD(&Q, &q);

    // calculate the observables
    Observables Obs;
  
    calcObservablesParity(&Int, par, &Q, &Obs);

    ensuredir("MIN");
    char outfile[1024];
    sprintf(outfile, "MIN/%s-minp.%d", parafile, (int) time(NULL));

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
  
    fprintf(outfp, "\n# einitial: %8.3f MeV\n", hbc*einitial);
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

    // normalize SlaterDet
    normalizeSlaterDet(&Q, &X);

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

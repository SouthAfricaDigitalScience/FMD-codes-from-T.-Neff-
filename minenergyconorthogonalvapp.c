/**

  \file minenergyconorthogonalvapp.c

  minimize projected energy orthogonal to given multiconfig state(s) 


  (c) 2003-2009 Thomas Neff

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
#include "fmd/ConstraintNOsci.h"
#include "fmd/ConstraintR2.h"
#include "fmd/ConstraintDipole.h"
#include "fmd/ConstraintQuadrupole.h"
#include "fmd/ConstraintOctupole.h"

#include "fmd/CenterofMass.h"
#include "fmd/SpatialOrientation.h"
#include "fmd/Observables.h"

#include "MinimizerDONLP2orthogonalvapp.h" 

#include "misc/physics.h"
#include "misc/utils.h"

#include "numerics/zcw.h"

#ifdef MPI
#include <mpi.h>
#include "fmdmpi/Communication.h"
#include "fmdmpi/MinimizerprojSlave.h"
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
    MinimizerprojSlave();

    MPI_Finalize();
  } else {
#endif

  int c;
  int cm=1;
  int ifix=1;
  int j=0;
  int par=0;
  double alpha=0.1;
  int overwrite=0;
  int log=0;
  int maxsteps=250;
  double shakemag=0.0;
  int constcm=1, constT2=0, constS2=0, constL2=0, constLS=0, constJ2=0;
  int constj2=0;
  int constnosci=0, constpnosci=0, constnnosci=0;
  int constmradius=0, constmquadrupole=0, constmoctupole=0;
  int consteradius=0, constedipole=0, constequadrupole=0, consteoctupole=0;
  int constnradius=0, constnquadrupole=0, constnoctupole=0;
  double constT2val=0.0, constS2val=0.0;
  double constL2val=0.0, constLSval=0.0, constJ2val=0.0;
  double constj2val=0.0;
  double constnoscival=0.0, constpnoscival=0.0, constnnoscival=0.0;
  double constmrval=0.0, constmqval=0.0, constmoval=0.0;
  double consterval=0.0, constedval=0.0, consteqval=0.0, consteoval=0.0;
  double constnrval=0.0, constnqval=0.0, constnoval=0.0;

  Constraint Const[MAXCONSTRAINT];
  int nconst = 0;

  if (argc < 5) {
    fprintf(stderr, "\nusage: %s PROJPARA interaction fixslaterdetfile slaterdetfile\n"
            "\n   -i I            optimize orthogonal to first I eigenstates"
	    "\n   -j J            project onto angular momentum"
	    "\n   -p[+1|-1]       project onto parity"
            "\n   -a ALPHA        minimize (emulti + ALPHA eintr)"
	    "\n   -l EVERY        log EVERY step"
	    "\n   -o              overwrite in all cases"
	    "\n   -m MAXSTEPS     maximum number of steps (default %d)"
	    "\n   -s MAGNITUDE    shake parameters before minimization"
	    "\n   -T T2           constrain isospin"
	    "\n   -L l2           constrain single-particle l2"
	    "\n   -S ls           constrain single-particle ls"
	    "\n   -J J2           constrain angular momentum"
	    "\n   -j j2           constrain single-particle angular momentum"
	    "\n   -N NOSCI        constrain number of oscillator quanta"	
            "\n   -R RADIUS       constrain radius"
	    "\n   -D DIPOLE       constrain dipole moments"
	    "\n   -Q QUADRUPOLE   constrain quadrupole moments"
	    "\n   -O OCTUPOLE     constrain octupole moments\n"
	    // "\n   -I Jz           cranking constrain"
	    , argv[0], maxsteps);
    cleanup(-1);
  }
  
  while((c = getopt(argc, argv, "i:j:p:a:Lol:m:s:T:L:S:J:N:R:D:Q:O:")) != -1)
    switch (c) {
    case 'i':
      ifix = atoi(optarg);
      break;
    case 'j':
      j = atoi(optarg);
      break;
    case 'p':
      par = atoi(optarg);
      break;
    case 'a':
      alpha = atof(optarg);
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
      constS2 = 1;
      constS2val = atof(optarg);
      break;
    case 'J':
      constJ2 = 1;
      constJ2val = atof(optarg);
      break;
      /*
    case 'j':
      constj2 = 1;
      constj2val = atof(optarg);
      break;
      */
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
      if (optarg[0] == 'E' || optarg[0] == 'P') {
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
    }
  
  if (constcm) {
    Const[nconst] = ConstraintCM;
    nconst++;
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

  int i;
  char* projpar = argv[optind];
  char* interactionfile = argv[optind+1];

  int nfix = argc-optind-3;
  char* fixmbfile[nfix];
  for (i=0; i<nfix; i++)
    fixmbfile[i] = argv[optind+2+i];

  char* parafile = argv[optind+2+nfix];

  // set up angular integration
  angintegrationpara AngPar;
  initAngintegration(&AngPar, projpar);

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    cleanup(-1);
  Int.cm = cm;

  // consider single fixed configuration for now
  SlaterDet Qfix[nfix];
  for (i=0; i<nfix; i++) 
    if (readSlaterDetfromFile(&Qfix[i], fixmbfile[i])) {
      fprintf(stderr, "couldn't read %s\n", fixmbfile[i]);
      cleanup(-1);
    }

  // configuration to be varied
  Parameterization P;
  Para qinitial;
  if (readParafromFile(&P, &qinitial, parafile)) {
    fprintf(stderr, "couldn't read %s\n", parafile);
    cleanup(-1);
  }

  SlaterDet Q;
  P.ParainitSlaterDet(&qinitial, &Q);
  P.ParatoSlaterDet(&qinitial, &Q);

  if (!par)
    par = (Q.A % 2) ? -1 : +1;

#ifdef MPI
  int task=TASKSTART;
  BroadcastTask(&task);

  BroadcastInteraction(&Int);
  BroadcastA(&Q.A);
#endif

  initWork(nfix, Qfix, j, par, &AngPar, &Int, &Q);
  
  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);

  moveboostorientSlaterDet(&Q, &X);

  double eintri, eorthoi;
  calcprojectedHamiltonian(nfix, Qfix, &Q, &Int, j, par, ifix, &AngPar,
                           &eintri, &eorthoi);

  for (int i=0; i<nconst; i++)
    fprintf(stderr, "# constraining %4s to %8.3f\n", 
	    Const[i].label, Const[i].output(Const[i].val)); 

  double einitial = eorthoi + alpha*eintri;

  fprintf(stderr, "\ninitial:\tE = %8.3f MeV, Eorthoi = %8.3f MeV, Eintr = %8.3f MeV\n\n", 
	  hbc*einitial, hbc*eorthoi, hbc*eintri);

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
  MinimizeDONLP2orthogonalvapp(&Int, j, par, ifix, alpha,
			       &AngPar, 
			       Const, nconst,
			       nfix, Qfix,
			       &P, &q, maxsteps,
			       log, parafile); 

  P.ParatoSlaterDet(&q, &Q);

  double eintrf, eorthof;

  calcprojectedHamiltonian(nfix, Qfix, &Q, &Int, j, par, ifix, &AngPar,
                           &eintrf, &eorthof);

  double e = eorthof + alpha*eintrf;

  fprintf(stderr, "\nfinal:\tE = %8.3f MeV, Eortho = %8.3f, Eintr = %8.3f MeV\n\n", 
	  hbc*e, hbc*eorthof, hbc*eintrf);

  moveboostorientSlaterDet(&Q, &X);

  calcprojectedHamiltonian(nfix, Qfix, &Q, &Int, j, par, ifix, &AngPar,
                           &eintrf, &eorthof);

  e = eorthof + alpha*eintrf;

  fprintf(stderr, "\nboosted:\tE = %8.3f MeV, Eortho = %8.3f, Eintr = %8.3f MeV\n\n", 
	  hbc*e, hbc*eorthof, hbc*eintrf);

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
  
    calcObservables(&Int, &Q, &Obs);

    ensuredir("MIN");
    char outfile[1024];
    sprintf(outfile, "MIN/%s-minortogonalvapp.%d", parafile, (int) time(NULL));

    fprintf(stderr, "... writing Parameters to file %s\n", outfile);

    FILE* outfp;
    if (!(outfp = fopen(outfile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", outfile);
      cleanup(-1);
    }

    fprintinfo(outfp);

    fprintf(outfp, "\n# minimized %s for %s in %s parameterization\n"
            "# projected on J^pi = %s\n",
            cm ? "< Hintr >" : "< H >", q.name, P.name, 
            AngmomtoStr(j, par == +1 ? 0 : 1));
    for (int ifix=0; ifix<nfix; ifix++)
      fprintf(outfp, "# fixed config: %s\n", fixmbfile[ifix]);
    fprintf(outfp, "# orthogonal to first %d eigenstates\n"
            "# using %s interaction\n", 
            ifix, Int.name);
  
    fprintf(outfp, "\n# initial: e = %8.3f MeV, eortho = %8.3f MeV, eintr = %8.3f MeV\n", 
            hbc*einitial, hbc*eorthoi, hbc*eintri);
    fprintf(outfp, "# final:   e = %8.3f MeV, eortho = %8.3f MeV, eintr = %8.3f MeV\n", 
            hbc*e, hbc*eorthof, hbc*eintrf);

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
  
    calcObservables(&Int, &Q, &Obs);

    fprintf(stderr, "... writing Parameters to file %s\n", parafile);

    FILE* outfp;
    if (!(outfp = fopen(parafile, "w"))) {
      fprintf(stderr, "couldn't open %s for writing\n", parafile);
      cleanup(-1);
    }

    fprintinfo(outfp);

    fprintf(outfp, "\n# minimized %s for %s in %s parameterization\n"
            "# projected on J^pi = %s\n",
            cm ? "< Hintr >" : "< H >", q.name, P.name, 
            AngmomtoStr(j, par == +1 ? 0 : 1));
    for (int ifix=0; ifix<nfix; ifix++)
      fprintf(outfp, "# fixed config: %s\n", fixmbfile[ifix]);
    fprintf(outfp, "# orthogonal to first %d eigenstates\n"
            "# using %s interaction\n", 
            ifix, Int.name);
  
    fprintf(outfp, "\n# initial: e = %8.3f MeV, eortho = %8.3f MeV, eintr = %8.3f MeV\n", 
            hbc*einitial, hbc*eorthoi, hbc*eintri);
    fprintf(outfp, "# final:   e = %8.3f MeV, eortho = %8.3f MeV, eintr = %8.3f MeV\n", 
            hbc*e, hbc*eorthof, hbc*eintrf);

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


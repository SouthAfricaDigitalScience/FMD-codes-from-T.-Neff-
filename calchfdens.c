/**

  \file calchfdens.c

  calculate the HF single-particle states and densities of a FMD Slater determinant and plot the densities


  (c) 2005 Thomas Neff

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>

#include "fmd/SlaterDet.h"
#include "fmd/CenterofMass.h"
#include "fmd/Potential.h"
#include "fmd/Hamiltonian.h"
#include "fmd/AngularMomenta.h"
#include "fmd/Isospin.h"
#include "fmd/NOscillator.h"

#include "numerics/cmat.h"
#include "misc/utils.h"
#include "misc/physics.h"

#define SQR(x)	(x)*(x)

void expect(int n, complex double* A, complex double* V, double* a)
{
  int i;
  int k,m;

  for (i=0; i<n; i++) {
    a[i] = 0.0;
    for (m=0; m<n; m++)
      for (k=0; k<n; k++)
	a[i] += conj(V[k+i*n])*A[k+m*n]*V[m+i*n];
  }
}


double angroot(double j2)
{
  return (-0.5+sqrt(j2+0.25));
}


typedef struct {
  int idx;
  double rank;
  double e;
  double j2;
  double l2;
  double hosci;
  double t3;
} SingleParticleState;


int cmpSingleParticleState(SingleParticleState* spa, SingleParticleState* spb)
{
  if (spa->rank < spb->rank)
    return +1;
  else
    return -1;
}


int main(int argc, char *argv[])
{
  createinfo(argc, argv);
  
  if (argc < 3) {
    fprintf(stderr, "\nusage: %s interaction slaterdetfile\n"
	    "\n   -s           sort sp states according to n, l, j"
	    "\n   -V           use only potential"
	    "\n   -O OMEGA     use oscillator constant [MeV]"
	    "\n   -T           T-Tcm"
				"\n   -p            no frame, no annotation"
				"\n   -g            debug"
				"\n   -q            quiet"
				"\n   -C LOGO       show logo"
				"\n   -o            reorient"
				"\n   -n            plain"
				"\n   -l LABEL	use label"
				"\n   -v {xy|yz|xz} select view"
				"\n   -x RANGE      coordinate range"
				"\n   -m RANGE	momentum range\n",
		argv[0]);
    exit(-1);
  }

  int c;
  int nljorder=0;
  int potentialonly=0;
  int cm=0;
  double omega=0.0;
  int debug=0, quiet=0;
  int coordinate=0; int momentum=0;
  int orient=0;
  int view[3] = {1, 1, 1};
  int nplots[2] = {2, 3};
  int plain=0;
  
  double xmax=5.5;
  double pmax=3.5;
  const int npoints=40;
  char* logo = NULL;
  char* label = NULL;
  
    
  while((c = getopt(argc, argv, "sVTO:v:c:oc:m:nl:ogC:qv:")) != -1)
    switch (c) {
    case 's':
      nljorder = 1;
      break;
    case 'V':
      potentialonly=1;
      break;
    case 'T':
      cm = 1;
      break;
    case 'O':
      omega = atof(optarg)/hbc;
      break;
		case 'g':
			debug=1;
			break;
		case 'q':
			quiet=1;
			break;
		case 'n':
			plain=1;
			break;
		case 'o':
			orient=1;
			break;
		case 'C':
			logo = optarg;
			break;
		case 'l':
			label=optarg;
			break;
		case 'c':
			coordinate = 1;
			xmax = atof(optarg);
			break;
		case 'm':
			momentum = 1;
			pmax = atof(optarg);
			break;
		case 'v':
			view[0]=0; view[1]=0; view[2]=0;
			nplots[1] = 1;
			if (!strcmp(optarg, "xy")) view[2]=1;
			if (!strcmp(optarg, "xz")) view[1]=1;
			if (!strcmp(optarg, "yz")) view[0]=1;
			break;
    }
	
	if (!(coordinate || momentum)) {
		coordinate = 1;
		momentum = 1;
	}

	nplots[0] = (coordinate && momentum) ? 2 : 1;
  
  char* interactionfile = argv[optind];
  char* slaterdetfile = argv[optind+1];

  Interaction Int;
  if (readInteractionfromFile(&Int, interactionfile))
    exit(-1);
  Int.cm = cm;

  SlaterDet Q;
  if (readSlaterDetfromFile(&Q, slaterdetfile))
    exit(-1);

  // number of nucleons
  int A = Q.A;

  // single-particle overlap matrix
  SlaterDetAux X;
  initSlaterDetAux(&Q, &X);
  calcSlaterDetAux(&Q, &X);

  if (orient) {
	  orientSlaterDet(&Q, &X);
	  calcSlaterDetAux(&Q, &X);
  }
  
  // Tcm
  double tcm;
  calcTCM(&Q, &X, &tcm);

  // derive oscillator constant from center of mass motion
  double omegacm = 4.0/3.0*tcm;

  if (omega == 0.0)
    omega = omegacm;

  complex double n[A*A];
  copycmat(A, X.n, n);

  // hartree-fock matrix in Gaussian single-particle basis

  complex double hhf[A*A];
  if (potentialonly)	calcPotentialHF(&Int, &Q, &X, hhf);
  else			calcHamiltonianHF(&Int, &Q, &X, hhf);
  
  complex double l2hf[A*A];
  calcl2HF(&Q, &X, l2hf);
  
  complex double j2hf[A*A];
  calcj2HF(&Q, &X, j2hf);

  complex double hoscihf[A*A];
  calcHOsciHF(&Q, &X, omega, hoscihf);

  complex double t3hf[A*A];
  calct3HF(&Q, &X, t3hf);
  
  // solve eigenvalue problem
  complex double lambda[A];
  complex double V[A*A];
  int dim;

  if (nljorder) {
    complex double nljhf[A*A];
    int i;
    for (i=0; i<A*A; i++)
      nljhf[i] = 10000*hoscihf[i]/omega - 100*j2hf[i] + l2hf[i];

      generalizedeigensystem(nljhf, n, A, 0.0, lambda, V, &dim);
  } else
    generalizedeigensystem(hhf, n, A, 0.0, lambda, V, &dim);

  // calculate expectation values
  double norm[A];
  expect(A, n, V, norm);

  double h[A];
  expect(A, hhf, V, h);

  double l2[A];
  expect(A, l2hf, V, l2);
  
  double j2[A];
  expect(A, j2hf, V, j2);

  double hosci[A];
  expect(A, hoscihf, V, hosci);

  double t3[A];
  expect(A, t3hf, V, t3);

  fprintinfo(stdout);

  // sort single-particle states by energy
  SingleParticleState sp[A];

  int i;
  for (i=0; i<A; i++) {
    sp[i].idx = i;
    sp[i].rank = lambda[i];
    sp[i].e = h[i]/norm[i];
    sp[i].j2 = j2[i]/norm[i];
    sp[i].l2 = l2[i]/norm[i];
    sp[i].hosci = hosci[i]/norm[i];
    sp[i].t3 = t3[i]/norm[i];
  }	
  
  qsort(sp, A, sizeof(SingleParticleState), cmpSingleParticleState);

	if (!label)
		label = nucleusIDLformat(nucleusname(Q.A, Q.Z));
  
	double* dens = malloc(3*SQR(npoints)*A*sizeof(double));
  
	// write densities to data files

	char datafile[255];
	FILE *datafp;
	int v;
	int j,k,m;

	int idx;
	for (idx=0; idx<A; idx++) {
		snprintf(datafile, 255, "%s.%d.hfdens", slaterdetfile, idx);
		if (!(datafp = fopen(datafile, "w"))) {
			fprintf(stderr, "couldn't open %s for writing\n", datafile);
			exit(-1);
		}

		// coordinate space densities
		if (coordinate) {
			for (v=0; v<3; v++)
				if (view[v]) {
				complex double* densxhf = malloc(3*SQR(npoints)*A*A*sizeof(complex double));
				calcDensitiesCoordinateHF(&Q, &X, v, npoints, xmax, densxhf);
				
				for (j=0; j<SQR(npoints); j++) {
					for (i=0; i<A; i++) {
						dens[j+i*3*SQR(npoints)] = 0.0;
						for (m=0; m<A; m++)
							for (k=0; k<A; k++)
								dens[j+i*3*SQR(npoints)] += conj(V[k+i*A])*densxhf[j+k*3*SQR(npoints)+m*A*3*SQR(npoints)]*V[m+i*A];
					}
				}
				free(densxhf);
				
				for (j=0; j<npoints; j++) {
					for (k=0; k<npoints; k++)
						fprintf(datafp, "%f  ", dens[k+j*npoints+3*idx*SQR(npoints)]/rho0/norm[idx]);
					fprintf(datafp, "\n");
				}
				}
		}

	// momentum space densities
				if (momentum) {
					for (v=0; v<3; v++)
						if (view[v]) {
						complex double* densphf = malloc(3*SQR(npoints)*A*A*sizeof(complex double));
						calcDensitiesMomentumHF(&Q, &X, v, npoints, pmax, densphf);
						
						for (j=0; j<SQR(npoints); j++) {
							for (i=0; i<A; i++) {
								dens[j+i*3*SQR(npoints)] = 0.0;
								for (m=0; m<A; m++)
									for (k=0; k<A; k++)
										dens[j+i*3*SQR(npoints)] += conj(V[k+i*A])*densphf[j+k*3*SQR(npoints)+m*A*3*SQR(npoints)]*V[m+i*A];
							}
						}
						free(densphf);
						
						for (j=0; j<npoints; j++) {
							for (k=0; k<npoints; k++)
								fprintf(datafp, "%f  ", dens[k+j*npoints+3*idx*SQR(npoints)]/norm[idx]);
							fprintf(datafp, "\n");
						}
						}
				}
		
	fclose(datafp);
	
	// writing IDL script

	char scriptfile[255];
	FILE* scriptfp;

	char* fmdhome;

	if (!(fmdhome = getenv("FMD"))) {
		fprintf(stderr, "environment variable FMD not defined\n");
		exit(-1);
	}

	snprintf(scriptfile, 255, "%s.%d.hfdens.script", slaterdetfile, idx);
	if (!(scriptfp = fopen(scriptfile, "w"))) {
		fprintf(stderr, "couldn't open %s for writing\n", scriptfile);
		exit(-1);
	}

	fprintf(scriptfp, 
			"; plot densities of nuclei\n; written by calchfdens\n\n");
	fprintf(scriptfp, "dens=fltarr(%d,%d,%d)\n", 
			npoints, npoints, nplots[0]*nplots[1]);

	fprintf(scriptfp, "!path = '%s/lib:' + !path\n", fmdhome);
	fprintf(scriptfp, "openr, unit, '%s', /get_lun\n", datafile);
	fprintf(scriptfp, "readf, unit, dens\n");
	fprintf(scriptfp, "free_lun, unit\n");
  
	fprintf(scriptfp, "\n.run multipost, densityplot\n");
	fprintf(scriptfp, 
			"pos = initpost('%s.%d.hfdens.eps',%d,%d, gapx=1.5, gapy=1.5, /color)\n",
			slaterdetfile, idx, nplots[0], nplots[1]);
	fprintf(scriptfp, "!x.minor = 4\n!y.minor = 4\n");
	fprintf(scriptfp, "loadct, 3\n");

	i = -1;

	if (coordinate) {
		for (v=0; v<3; v++) {
			if (view[v]) {
				i++;
				fprintf(scriptfp, "!p.position = pos(%d,*)\n", i);
				fprintf(scriptfp, "densityplot, dens(*,*,%d), %d, %f, %f, %f, %f, $\n",
						i, v, -xmax, xmax, -xmax, xmax);
				fprintf(scriptfp, "\t/cut, /coordinate, /cont, /dens");
				if (plain) fprintf(scriptfp, ", /noframe, /noannot");
				else {
					fprintf(scriptfp, ", key = '%s'", label);
					if (logo && nplots[0] == 1) fprintf(scriptfp, ", logo = '%s'", logo);
				} 
				fprintf(scriptfp, "\n");
			}
		}
	}

	if (momentum) {
		for (v=0; v<3; v++) {
			if (view[v]) {
				i++;
				fprintf(scriptfp, "!p.position = pos(%d,*)\n", i);
				fprintf(scriptfp, "densityplot, dens(*,*,%d), %d, %f, %f, %f, %f, $\n",
						i, v, -pmax, pmax, -pmax, pmax);
				fprintf(scriptfp, "\t/cut, /momentum, /cont, /dens");
				if (plain) fprintf(scriptfp, ", /noframe, /noannot");
				else {
					fprintf(scriptfp, ", key = '%s'", label);
					if (logo) fprintf(scriptfp, ", logo = '%s'", logo);
				}
				fprintf(scriptfp, "\n");
			}
		}
	}

	fclose(scriptfp);

  // calling IDL

	char call[255];

	snprintf(call, 255, "idl < %s", scriptfile);
	system(call);

	char epsfile[255];
	snprintf(epsfile, 255, "%s.%d.hfdens.eps", slaterdetfile, idx);

  // calling gv

	if (!quiet) {
		snprintf(call, 255, "gv %s &", epsfile);
		system(call);
	}

  // clean up

	if (!debug) {
		remove(datafile); 
		remove(scriptfile);
	}
	}
	
	  // print proton and neutron levels

	fprintf(stdout, "\nusing oscillator constant: hbar Omega = %8.3f MeV\n",
			omega*hbc);

	fprintf(stdout, "\nproton levels:\n");
	for (i=0; i<A; i++)
		if (sp[i].t3 > 0.0) {
		fprintf(stdout, "[%d]\te: %8.3f MeV, j: %5.3f, l: %5.3f, nosci: %5.3f\n",
				sp[i].idx, hbc*sp[i].e, angroot(sp[i].j2), angroot(sp[i].l2),
				sp[i].hosci/omega-1.5);
		}

		fprintf(stdout, "\nneutron levels:\n");
		for (i=0; i<A; i++)
			if (sp[i].t3 < 0.0) {
			fprintf(stdout, "[%d]\te: %8.3f MeV, j: %5.3f, l: %5.3f, nosci: %5.3f\n",
					sp[i].idx , hbc*sp[i].e, angroot(sp[i].j2), angroot(sp[i].l2),
					sp[i].hosci/omega-1.5);
			}
  return 0;
}

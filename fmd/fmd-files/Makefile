include ../Makefile.inc

CFLAGS := $(CFLAGS) -I..

OBJLIBS = ../libfmd.a
OBJS 	= Gaussian.o SlaterDet.o Ovlap.o \
	  CenterofMass.o KineticEnergy.o Interaction.o Potential.o \
	  Oscillator.o NOscillator.o \
	  Radii.o SDRadii.o RadiiAll.o RadiiLS.o Quadrupole.o \
	  AngularMomenta.o Isospin.o Parity.o TimeReversal.o NOsci.o \
	  Observables.o SpatialOrientation.o \
	  Densities.o Densities3d.o SpinDensities.o \
	  gradGaussian.o gradSlaterDet.o gradCMSlaterDet.o \
	  gradCenterofMass.o gradKineticEnergy.o gradPotential.o \
	  ConstraintCMXP.o ConstraintT2.o ConstraintS2.o \
	  ConstraintCM.o ConstraintJ.o ConstraintJ2.o \
	  Constraintl2.o ConstraintJ2.o \
	  ConstraintR2.o ConstraintLS.o \
	  ConstraintNOsci.o \
	  ConstraintDipole.o ConstraintQuadrupole.o ConstraintOctupole.o ConstraintQuadrupole-detEQ.o\
	  ConstraintSymmetries.o \
          Energies.o \
	  Hamiltonian.o gradHamiltonian.o gradOscillator.o \
	  Parameterization.o ParameterizationFMD.o ParameterizationFMDS.o \
	  ParameterizationFMDr.o ParameterizationFMDreal.o \
	  ParameterizationAMD.o ParameterizationAMDd.o ParameterizationAMDS.o \
	  ParameterizationAMDA.o ParameterizationAMDAS.o \
	  ParameterizationAlpha.o ParameterizationAlphaA.o \
	  ParameterizationAlphaC.o \
	  ParameterizationCluster.o \
	  ParameterizationCoreFMD.o \
	  ParameterizationClusterFMD.o ParameterizationClustersFMD.o \
	  ParameterizationFMDd3h.o ParameterizationFMDvxz.o \
	  Projection.o Symmetry.o cmprojection.o angprojection.o \
	  ProjectedObservables.o ElectroMagneticMultipole.o \
	  GamovTeller.o \
	  Formfactors.o \
	  OneNucleonOvlaps.o TwoNucleonOvlaps.o \
	  TwoBodyDensity.o \
	  MultiSlaterDet.o SimpleSlaterDet.o SymmetricSlaterDet.o \
	  SymmetricMultiSlaterDet.o DiClusterProjSlaterDet.o \
	  DiClusterMultiProjSlaterDet.o DiClusterMulticonfig.o \
	  HOBasis.o DensityMatrixHO.o ProjectedDensityMatrixHO.o \
	  MeanOscillatorQuanta.o ProjectedShellOccupations.o


all:	$(OBJLIBS)


# create dependencies

include $(OBJS:.o=.d)

../libfmd.a:	$(OBJS)
	$(AR) $(ARFLAGS) $@ $? 


clean:
	$(RM) $(OBJS)

veryclean:
	$(RM) $(OBJS) $(OBJS:.o=.d)

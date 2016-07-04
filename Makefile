include Makefile.inc

CFLAGS := $(CFLAGS) -I.

LIBS = 		git_version.o -lfmd -lnumerics -lmisc $(LAPACKLIBS) $(SYSLIBS)
LIBSMPI =	git_version.o -lfmdmpi -lfmd -lnumerics -lmisc $(LAPACKLIBS) $(MPILIBS) $(SYSLIBS)

DIRS	= 	fmd fmdmpi numerics misc
OBJLIBS = 	git_version.o libfmd.a libnumerics.a libmisc.a
OBJLIBSMPI =	git_version.o libfmdmpi.a libfmd.a libnumerics.a libmisc.a

OBJS =		gnd2fmdpara.o sldet2fmdpara.o transformsldet.o joinsldets.o \
		minenergy.o minenergyp.o MinimizerLBFGS.o \
		minenergycon.o MinimizerDONLP2.o minenergycon-detQ.o minenergycon-detEQ.o\
          	minenergyconp.o MinimizerDONLP2p.o \
		minenergyconpar.o MinimizerDONLP2par.o \
		minenergyconproj.o MinimizerDONLP2proj.o\
                minenergyconproj-detEQ.o \
		minenergyconvap.o MinimizerDONLP2vap.o \
                minenergyconvap-detEQ.o \
		minenergyconvapp.o MinimizerDONLP2vapp.o \
		minenergyconvappc.o MinimizerDONLP2vappc.o \
		minenergyconvappcm.o MinimizerDONLP2vappcm.o \
		minenergyconvappiso.o MinimizerDONLP2vappiso.o \
		minenergyconmultivapp.o MinimizerDONLP2multivapp.o \
		minenergyconmultivappcm.o MinimizerDONLP2multivappcm.o \
		minenergyconorthogonalvapp.o MinimizerDONLP2orthogonalvapp.o \
		minenergyconorthogonalproj.o MinimizerDONLP2orthogonalproj.o \
		calcenergy.o calcenergyp.o calcnorm.o calcnormp.o \
		calcquadrupole.o calcorientation.o calcconstraints.o \
		recalcenergy.o recalcenergyp.o \
		calchflevels.o \
		calcovlapproj.o calcenergyprojme.o calcovlapprojme.o \
		calcenergyproj.o calcenergyprojks.o calcenergymultiproj.o \
		recalcenergymultiproj.o \
		calcenergymultiselproj.o calcenergymultiprojsel.o \
		calcbasisovlap.o \
		calcdensities.o calcdensitiespn.o calcdensitiespnp.o \
		calcdensitiespmn.o calcdensitiesspin.o \
		calcdensities3d.o calcdens.o calchfdens.o \
		calcradiiall.o calcradiiallme.o calcradiils.o \
		calcsdradii.o \
		calctransitions.o calctransitionsme.o calctransitionsgt.o \
		calcchargeformfactors.o calctransitionchargeformfactors.o \
		calcpointformfactors.o calctransitionpointformfactors.o \
		calcformfactorsme.o \
		calcpointdensities.o calctransitionpointdensities.o \
		calcmatterdensities.o \
		calcchargedensities.o calctransitionchargedensities.o \
		calconenucleonovlaps.o \
		calctwonucleonovlapst.o calctwonucleonovlapsy.o \
		calcpairs.o \
		calctwobodydensitiescoord.o calctwobodydensitiescoordl.o \
		calctwobodydensitiesmom.o calctwobodydensitiesmoml.o \
		calctwobodydensitiescoordme.o calctwobodydensitiescoordlme.o \
		calctwobodydensitiesmomme.o calctwobodydensitiesmomlme.o \
		calcenergyprojmulti.o calcenergymultiprojmulti.o \
		calcenergyprojmemulti.o calcovlapprojmemulti.o \
		calctransitionsmulti.o calctransitionsmemulti.o \
		calcbasisovlapmulti.o \
		calcoccupationnumbersho.o \
		calcoccupationnumbershoproj.o calctransitionnumbershoproj.o \
		calcoccupationnumbershoprojme.o calcoccupationnumbershoprojmes.o \
		calcshelloccupations.o calcshelloccupationsme.o \
		calctimereversal.o \


OBJSMPI = 	minenergy.mpi.o minenergyp.mpi.o \
		minenergycon.mpi.o MinimizerDONLP2.mpi.o \
		minenergyconp.mpi.o MinimizerDONLP2p.mpi.o \
		minenergyconproj.mpi.o MinimizerDONLP2proj.mpi.o\
                minenergyconproj.mpi-detEQ.o \
		minenergyconvap.mpi.o MinimizerDONLP2vap.mpi.o \
                minenergyconvap-detEQ.mpi.o \
		minenergyconvapp.mpi.o MinimizerDONLP2vapp.mpi.o \
		minenergyconvappc.mpi.o MinimizerDONLP2vappc.mpi.o \
		minenergyconvappcm.mpi.o MinimizerDONLP2vappcm.mpi.o \
		minenergyconvappiso.mpi.o MinimizerDONLP2vappiso.mpi.o \
		minenergyconmultivapp.mpi.o MinimizerDONLP2multivapp.mpi.o \
		minenergyconmultivappcm.mpi.o MinimizerDONLP2multivappcm.mpi.o \
		minenergyconorthogonalvapp.mpi.o MinimizerDONLP2orthogonalvapp.mpi.o \
		minenergyconorthogonalproj.mpi.o MinimizerDONLP2orthogonalproj.mpi.o \
		minenergyconpar.mpi.o MinimizerDONLP2par.mpi.o \
                calcenergy.mpi.o \
		calcenergyproj.mpi.o calcenergymultiproj.mpi.o \
		calcenergymultiprojsel.mpi.o \
		calcenergyprojmulti.mpi.o calcenergymultiprojmulti.mpi.o \
		calconenucleonovlaps.mpi.o \
		calctwonucleonovlapst.mpi.o calctwonucleonovlapsy.mpi.o

BINARIES = 	gnd2fmdpara sldet2fmdpara transformsldet joinsldets \
	  	minenergy minenergyp \
		minenergycon minenergycon-detQ minenergyconp minenergyconpar \
		minenergyconproj minenergyconorthogonalproj \
                minenergyconproj-detEQ minenergycon-detEQ \
		minenergyconvap minenergyconvapp minenergyconvappiso minenergyconvap-detEQ \
		minenergyconvappc minenergyconvappcm \
		minenergyconmultivapp minenergyconmultivappcm \
		minenergyconorthogonalvapp \
	  	calcenergy calcenergyp calcnorm calcnormp \
		calcquadrupole calcorientation calcconstraints \
		recalcenergy recalcenergyp \
		calchflevels \
		calcenergyprojme calcovlapproj calcovlapprojme \
		calcenergyproj calcenergyprojks \
		calcenergymultiproj recalcenergymultiproj \
		calcenergymultiselproj calcenergymultiprojsel \
		calcbasisovlap \
	  	calcdensities calcdensitiespn calcdensitiespnp \
		calcdensitiespmn calcdensitiesspin \
		calcdensities3d calcdens calchfdens \
		calcradiiall calcsdradii calcradiiallme calcradiils \
		calctransitions calctransitionsme calctransitionsgt \
		calcchargeformfactors calctransitionchargeformfactors \
		calcpointformfactors calctransitionpointformfactors \
		calcformfactorsme \
		calcpointdensities calctransitionpointdensities \
		calcmatterdensities \
		calcchargedensities calctransitionchargedensities \
		calconenucleonovlaps \
		calctwonucleonovlapst calctwonucleonovlapsy \
		calcpairs \
		calctwobodydensitiescoord calctwobodydensitiescoordl \
		calctwobodydensitiesmom calctwobodydensitiesmoml \
		calctwobodydensitiescoordme calctwobodydensitiescoordlme \
		calctwobodydensitiesmomme calctwobodydensitiesmomlme \
		calcenergyprojmulti calcenergymultiprojmulti \
		calcenergyprojmemulti calcovlapprojmemulti \
		calctransitionsmulti calctransitionsmemulti \
		calcoccupationnumbersho \
		calcoccupationnumbershoproj calctransitionnumbershoproj \
		calcoccupationnumbershoprojme calcoccupationnumbershoprojmes \
		calcshelloccupations calcshelloccupationsme \
		calctimereversal \


BINARIESMPI =	mpiminenergy mpiminenergyp \
		mpiminenergycon mpiminenergyconp \
		mpiminenergyconproj mpiminenergyconorthogonalproj\
                mpiminenergyconproj-detEQ \
		mpiminenergyconvap mpiminenergyconvapp mpiminenergyconvappiso mpiminenergyconvap-detEQ \
		mpiminenergyconvappc mpiminenergyconvappcm \
		mpiminenergyconmultivapp mpiminenergyconmultivappcm \
		mpiminenergyconorthogonalvapp \
		mpiminenergyconpar \
                mpicalcenergy \
		mpicalcenergyproj mpicalcenergymultiproj \
		mpicalcenergymultiprojsel \
		mpicalcenergyprojmulti \
		mpicalcoccupationnumbershoprojmes \
		mpicalconenucleonovlaps \
		mpicalctwonucleonovlapst mpicalctwonucleonovlapsy


all: $(BINARIES) $(BINARIESMPI)


gnd2fmdpara:	gnd2fmdpara.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ gnd2fmdpara.o $(LIBS)

sldet2fmdpara:	sldet2fmdpara.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ sldet2fmdpara.o $(LIBS)

transformsldet:	transformsldet.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ transformsldet.o $(LIBS)

joinsldets:	joinsldets.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ joinsldets.o $(LIBS)

minenergy:	minenergy.o MinimizerLBFGS.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergy.o MinimizerLBFGS.o $(LIBS)

mpiminenergy:	minenergy.mpi.o MinimizerLBFGS.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergy.mpi.o MinimizerLBFGS.o $(LIBSMPI)

minenergyp:	minenergyp.o MinimizerLBFGS.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyp.o MinimizerLBFGS.o $(LIBS)

mpiminenergyp:	minenergyp.mpi.o MinimizerLBFGS.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyp.mpi.o MinimizerLBFGS.o $(LIBSMPI)

minenergycon:	minenergycon.o MinimizerDONLP2.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergycon.o MinimizerDONLP2.o $(LIBS)

minenergycon-detQ:	minenergycon-detQ.o MinimizerDONLP2.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergycon-detQ.o MinimizerDONLP2.o $(LIBS)

minenergycon-detEQ:	minenergycon-detEQ.o MinimizerDONLP2.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergycon-detEQ.o MinimizerDONLP2.o $(LIBS)

minenergycon-detEQ-gradient:	minenergycon-detEQ-gradient.o MinimizerDONLP2.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergycon-detEQ-gradient.o MinimizerDONLP2.o $(LIBS)

mpiminenergycon:	minenergycon.mpi.o MinimizerDONLP2.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergycon.mpi.o MinimizerDONLP2.mpi.o $(LIBSMPI)

minenergyconp:	minenergyconp.o MinimizerDONLP2p.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconp.o MinimizerDONLP2p.o $(LIBS)

mpiminenergyconp:	minenergyconp.mpi.o MinimizerDONLP2p.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconp.mpi.o MinimizerDONLP2p.mpi.o $(LIBSMPI)

minenergyconproj:	minenergyconproj.o MinimizerDONLP2proj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconproj.o MinimizerDONLP2proj.o $(LIBS)

minenergyconproj-detEQ:	minenergyconproj-detEQ.o MinimizerDONLP2proj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconproj-detEQ.o MinimizerDONLP2proj.o $(LIBS)

mpiminenergyconproj:	minenergyconproj.mpi.o MinimizerDONLP2proj.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconproj.mpi.o MinimizerDONLP2proj.mpi.o $(LIBSMPI)

mpiminenergyconproj-detEQ:	minenergyconproj-detEQ.mpi.o MinimizerDONLP2proj.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconproj-detEQ.mpi.o MinimizerDONLP2proj.mpi.o $(LIBSMPI)

minenergyconvap:	minenergyconvap.o MinimizerDONLP2vap.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvap.o MinimizerDONLP2vap.o $(LIBS)

minenergyconvap-detEQ:	minenergyconvap-detEQ.o MinimizerDONLP2vap.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvap-detEQ.o MinimizerDONLP2vap.o $(LIBS)

mpiminenergyconvap:	minenergyconvap.mpi.o MinimizerDONLP2vap.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvap.mpi.o MinimizerDONLP2vap.mpi.o $(LIBSMPI)

mpiminenergyconvap-detEQ:	minenergyconvap-detEQ.mpi.o MinimizerDONLP2vap.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvap-detEQ.mpi.o MinimizerDONLP2vap.mpi.o $(LIBSMPI)

minenergyconvapp:	minenergyconvapp.o MinimizerDONLP2vapp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvapp.o MinimizerDONLP2vapp.o $(LIBS)

mpiminenergyconvapp:	minenergyconvapp.mpi.o MinimizerDONLP2vapp.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvapp.mpi.o MinimizerDONLP2vapp.mpi.o $(LIBSMPI)

minenergyconvappc:	minenergyconvappc.o MinimizerDONLP2vappc.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvappc.o MinimizerDONLP2vappc.o $(LIBS)

mpiminenergyconvappc:	minenergyconvappc.mpi.o MinimizerDONLP2vappc.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvappc.mpi.o MinimizerDONLP2vappc.mpi.o $(LIBSMPI)

minenergyconvappcm:	minenergyconvappcm.o MinimizerDONLP2vappcm.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvappcm.o MinimizerDONLP2vappcm.o $(LIBS)

mpiminenergyconvappcm:	minenergyconvappcm.mpi.o MinimizerDONLP2vappcm.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvappcm.mpi.o MinimizerDONLP2vappcm.mpi.o $(LIBSMPI)

minenergyconvappiso:	minenergyconvappiso.o MinimizerDONLP2vappiso.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvappiso.o MinimizerDONLP2vappiso.o $(LIBS)

mpiminenergyconvappiso:	minenergyconvappiso.mpi.o MinimizerDONLP2vappiso.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvappiso.mpi.o MinimizerDONLP2vappiso.mpi.o $(LIBSMPI)

minenergyconvapn:	minenergyconvapn.o MinimizerDONLP2vapn.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconvapn.o MinimizerDONLP2vapn.o $(LIBS)

mpiminenergyconvapn:	minenergyconvapn.mpi.o MinimizerDONLP2vapn.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconvapn.mpi.o MinimizerDONLP2vapn.mpi.o $(LIBSMPI)

minenergyconmultivapp:	minenergyconmultivapp.o MinimizerDONLP2multivapp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconmultivapp.o MinimizerDONLP2multivapp.o $(LIBS)

mpiminenergyconmultivapp:	minenergyconmultivapp.mpi.o MinimizerDONLP2multivapp.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconmultivapp.mpi.o MinimizerDONLP2multivapp.mpi.o $(LIBSMPI)

minenergyconmultivappcm:	minenergyconmultivappcm.o MinimizerDONLP2multivappcm.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconmultivappcm.o MinimizerDONLP2multivappcm.o $(LIBS)

mpiminenergyconmultivappcm:	minenergyconmultivappcm.mpi.o MinimizerDONLP2multivappcm.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconmultivappcm.mpi.o MinimizerDONLP2multivappcm.mpi.o $(LIBSMPI)

minenergyconorthogonalvapp:	minenergyconorthogonalvapp.o MinimizerDONLP2orthogonalvapp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconorthogonalvapp.o MinimizerDONLP2orthogonalvapp.o $(LIBS)

mpiminenergyconorthogonalvapp:	minenergyconorthogonalvapp.mpi.o MinimizerDONLP2orthogonalvapp.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconorthogonalvapp.mpi.o MinimizerDONLP2orthogonalvapp.mpi.o $(LIBSMPI)

minenergyconorthogonalproj:	minenergyconorthogonalproj.o MinimizerDONLP2orthogonalproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconorthogonalproj.o MinimizerDONLP2orthogonalproj.o $(LIBS)

mpiminenergyconorthogonalproj:	minenergyconorthogonalproj.mpi.o MinimizerDONLP2orthogonalproj.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconorthogonalproj.mpi.o MinimizerDONLP2orthogonalproj.mpi.o $(LIBSMPI)

minenergyconpar:	minenergyconpar.o MinimizerDONLP2par.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ minenergyconpar.o MinimizerDONLP2par.o $(LIBS)

mpiminenergyconpar:	minenergyconpar.mpi.o MinimizerDONLP2par.mpi.o $(OBJLIBS)
	$(MPILD) $(MPILDFLAGS) -o $@ minenergyconpar.mpi.o MinimizerDONLP2par.mpi.o $(LIBSMPI)

calcenergy:	calcenergy.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergy.o $(LIBS)

mpicalcenergy:	calcenergy.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcenergy.mpi.o $(LIBSMPI)

calcenergyp:	calcenergyp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergyp.o $(LIBS)

calcnorm:	calcnorm.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcnorm.o $(LIBS)

calcnormp:	calcnormp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcnormp.o $(LIBS)

calcquadrupole:	calcquadrupole.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcquadrupole.o $(LIBS)

calcorientation:	calcorientation.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcorientation.o $(LIBS)

calcconstraints:	calcconstraints.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcconstraints.o $(LIBS)

recalcenergy:	recalcenergy.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ recalcenergy.o $(LIBS)

recalcenergyp:	recalcenergyp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ recalcenergyp.o $(LIBS)

calchflevels:	calchflevels.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calchflevels.o $(LIBS)

calcenergyprojme:	calcenergyprojme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergyprojme.o $(LIBS)

calcovlapproj:	calcovlapproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcovlapproj.o $(LIBS)

calcovlapprojme:	calcovlapprojme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcovlapprojme.o $(LIBS)

calcenergyproj:	calcenergyproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergyproj.o $(LIBS)

calcenergyprojks:	calcenergyprojks.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergyprojks.o $(LIBS)

mpicalcenergyproj:	calcenergyproj.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcenergyproj.mpi.o $(LIBSMPI)

calcenergymultiproj:	calcenergymultiproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergymultiproj.o $(LIBS)

recalcenergymultiproj:	recalcenergymultiproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ recalcenergymultiproj.o $(LIBS)

mpicalcenergymultiproj:	calcenergymultiproj.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcenergymultiproj.mpi.o $(LIBSMPI)

calcenergymultiselproj:	calcenergymultiselproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergymultiselproj.o $(LIBS)

calcenergymultiprojsel:	calcenergymultiprojsel.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergymultiprojsel.o $(LIBS)

mpicalcenergymultiprojsel:	calcenergymultiprojsel.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcenergymultiprojsel.mpi.o $(LIBSMPI)

calcbasisovlap:		calcbasisovlap.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcbasisovlap.o $(LIBS)

calcradiiall:		calcradiiall.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcradiiall.o $(LIBS)

calcradiiallme:		calcradiiallme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcradiiallme.o $(LIBS)

calcsdradii:		calcsdradii.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcsdradii.o $(LIBS)

calcradiils:		calcradiils.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcradiils.o $(LIBS)

calcdensities:		calcdensities.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdensities.o $(LIBS)

calcdensitiespn:	calcdensitiespn.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdensitiespn.o $(LIBS)

calcdensitiespmn:	calcdensitiespmn.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdensitiespmn.o $(LIBS)

calcdensitiespnp:	calcdensitiespnp.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdensitiespnp.o $(LIBS)

calcdensitiesspin:	calcdensitiesspin.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdensitiesspin.o $(LIBS)

calcdensities3d:	calcdensities3d.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdensities3d.o $(LIBS)

calcdens:		calcdens.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcdens.o $(LIBS)

calchfdens:		calchfdens.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calchfdens.o $(LIBS)

calctransitions:	calctransitions.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitions.o $(LIBS)

calctransitionsme:	calctransitionsme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionsme.o $(LIBS)

calctransitionsgt:	calctransitionsgt.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionsgt.o $(LIBS)

calcchargeformfactors:	calcchargeformfactors.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcchargeformfactors.o $(LIBS)

calctransitionchargeformfactors:	calctransitionchargeformfactors.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionchargeformfactors.o $(LIBS)

calctransitionchargeformfactorsod:	calctransitionchargeformfactorsod.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionchargeformfactorsod.o $(LIBS)

calcpointformfactors:	calcpointformfactors.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcpointformfactors.o $(LIBS)

calctransitionpointformfactors:	calctransitionpointformfactors.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionpointformfactors.o $(LIBS)

calcformfactormes:	calcformfactormes.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcformfactormes.o $(LIBS)

mpicalcformfactormes:	calcformfactormes.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcformfactormes.mpi.o $(LIBSMPI)

calcformfactorsme:	calcformfactorsme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcformfactorsme.o $(LIBS)

calcpointdensities:	calcpointdensities.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcpointdensities.o $(LIBS)

calctransitionpointdensities:	calctransitionpointdensities.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionpointdensities.o $(LIBS)

calcmatterdensities:	calcmatterdensities.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcmatterdensities.o $(LIBS)

calcchargedensities:	calcchargedensities.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcchargedensities.o $(LIBS)

calctransitionchargedensities:	calctransitionchargedensities.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionchargedensities.o $(LIBS)

calcspectroscopicamplitudes:	calcspectroscopicamplitudes.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcspectroscopicamplitudes.o $(LIBS)

calconenucleonovlaps:	calconenucleonovlaps.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calconenucleonovlaps.o $(LIBS)

mpicalconenucleonovlaps:	calconenucleonovlaps.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calconenucleonovlaps.mpi.o $(LIBSMPI)

calctwonucleonovlapst:	calctwonucleonovlapst.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwonucleonovlapst.o $(LIBS)

mpicalctwonucleonovlapst:	calctwonucleonovlapst.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calctwonucleonovlapst.mpi.o $(LIBSMPI)

calctwonucleonovlapsy:	calctwonucleonovlapsy.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwonucleonovlapsy.o $(LIBS)

mpicalctwonucleonovlapsy:	calctwonucleonovlapsy.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calctwonucleonovlapsy.mpi.o $(LIBSMPI)

calcpairs:	calcpairs.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcpairs.o $(LIBS)

calctwobodydensitiescoord:	calctwobodydensitiescoord.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiescoord.o $(LIBS)

calctwobodydensitiescoordl:	calctwobodydensitiescoordl.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiescoordl.o $(LIBS)

calctwobodydensitiesmom:	calctwobodydensitiesmom.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiesmom.o $(LIBS)

calctwobodydensitiesmoml:	calctwobodydensitiesmoml.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiesmoml.o $(LIBS)

calctwobodydensitiescoordme:	calctwobodydensitiescoordme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiescoordme.o $(LIBS)

calctwobodydensitiescoordlme:	calctwobodydensitiescoordlme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiescoordlme.o $(LIBS)

calctwobodydensitiesmomme:	calctwobodydensitiesmomme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiesmomme.o $(LIBS)

calctwobodydensitiesmomlme:	calctwobodydensitiesmomlme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctwobodydensitiesmomlme.o $(LIBS)

calcenergyprojmulti:	calcenergyprojmulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergyprojmulti.o $(LIBS)

calcenergyprojmemulti:	calcenergyprojmemulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergyprojmemulti.o $(LIBS)

calcenergymultiprojmulti:	calcenergymultiprojmulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcenergymultiprojmulti.o $(LIBS)

calctransitionsmulti:	calctransitionsmulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionsmulti.o $(LIBS)

calctransitionsmemulti:	calctransitionsmemulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionsmemulti.o $(LIBS)

calcovlapprojmemulti:	calcovlapprojmemulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcovlapprojmemulti.o $(LIBS)

calcbasisovlapmulti:	calcbasisovlapmulti.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcbasisovlapmulti.o $(LIBS)

mpicalcenergyprojmulti:	calcenergyprojmulti.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcenergyprojmulti.mpi.o $(LIBSMPI)

mpicalcenergymultiprojmulti:	calcenergymultiprojmulti.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcenergymultiprojmulti.mpi.o $(LIBSMPI)

createSymmetricSlaterDet:	createSymmetricSlaterDet.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ createSymmetricSlaterDet.o $(LIBS)

createjoinedSymmetricSlaterDet:	createjoinedSymmetricSlaterDet.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ createjoinedSymmetricSlaterDet.o $(LIBS)

calcoccupationnumbersho:	calcoccupationnumbersho.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcoccupationnumbersho.o $(LIBS)

calcoccupationnumbershoproj:	calcoccupationnumbershoproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcoccupationnumbershoproj.o $(LIBS)

calctransitionnumbershoproj:	calctransitionnumbershoproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctransitionnumbershoproj.o $(LIBS)

calcoccupationnumbershoprojme:	calcoccupationnumbershoprojme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcoccupationnumbershoprojme.o $(LIBS)

calcoccupationnumbershoprojmes:	calcoccupationnumbershoprojmes.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcoccupationnumbershoprojmes.o $(LIBS)

mpicalcoccupationnumbershoprojmes:	calcoccupationnumbershoprojmes.mpi.o $(OBJLIBSMPI)
	$(MPILD) $(MPILDFLAGS) -o $@ calcoccupationnumbershoprojmes.mpi.o $(LIBSMPI)

calctimereversal:	calctimereversal.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calctimereversal.o $(LIBS)

calcmeanoscillatorquanta:	calcmeanoscillatorquanta.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcmeanoscillatorquanta.o $(LIBS)

calcshelloccupations:	calcshelloccupations.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcshelloccupations.o $(LIBS)

calcshelloccupationsme:	calcshelloccupationsme.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ calcshelloccupationsme.o $(LIBS)


libfmd.a:	force
	cd fmd; $(MAKE) $(MFLAGS)

libfmdmpi.a:	force
	cd fmdmpi; $(MAKE) $(MFLAGS)

libnumerics.a:	force
	cd numerics; $(MAKE) $(MFLAGS)

libmisc.a:	force
	cd misc; $(MAKE) $(MFLAGS)

# on every build, record the working copy revision string

git_version.c:	force
	@echo -n 'const char* git_version(void) { const char* GIT_Version = "' > git_version.c
	git show -s --pretty=format:"%h" >> git_version.c
	@echo ' ("__DATE__" "__TIME__")"; return GIT_Version; }' >> git_version.c

# create .gitignore

.gitignore:	force
	@echo "# created by Makefile" > $@
	@echo ".gitignore" >> $@
	@echo "*~\n*.swp\n*.d\n*.o\n*.a" >> $@
	@echo "Makefile.inc" >> $@
	@echo "git_version.c" >> $@
	@echo "doc" >> $@
	@echo $(BINARIES) | sed 's/ /\n/g' >> $@
	@echo $(BINARIESMPI) | sed 's/ /\n/g' >> $@


# targets depending on force will be made always

force:
	true


# create dependencies

include $(OBJS:.o=.d)



install:	$(BINARIES) $(BINARIESMPI)
ifndef FMD
		$(ECHO) "Environment Variable FMD not defined"
else
		mkdir -p $(FMD)/bin $(FMD)/data $(FMD)/lib
		$(INSTALL) $(BINARIES) $(FMD)/bin
		$(INSTALL) $(BINARIESMPI) $(FMD)/bin
		$(INSTALL) data/* $(FMD)/data
		$(INSTALL) lib/* $(FMD)/lib
endif

doc:
	git log -v > doc/CHANGES
	doxygen doc/Doxyfile

clean:
	$(RM) $(OBJS) $(OBJSMPI) $(OBJLIBS) $(OBJLIBSMPI)
	for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

veryclean:
	$(RM) $(OBJS) $(OBJSMPI) $(OBJLIBS) $(OBJLIBSMPI)
	$(RM) $(OBJS:.o=.d)
	$(RM) $(BINARIES) $(BINARIESMPI)
	for d in $(DIRS); do (cd $$d; $(MAKE) veryclean ); done

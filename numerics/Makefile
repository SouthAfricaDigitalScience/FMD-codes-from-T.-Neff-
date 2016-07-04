include ../Makefile.inc

OBJLIBS = ../libnumerics.a
COBJS 	= cmat.o rotationmatrices.o coulomb.o clebsch.o \
		legendrep.o sphericalharmonics.o sphericalbessel.o \
		zcw.o gaussquad.o interpol.o
FOBJS	= zdet.o djmnb.o lbfgs.o iqd.o dcsint.o coulcc.o
OBJS	= $(COBJS) $(FOBJS) donlp2.o

all:	$(OBJLIBS)


# special treatment for donlp2

donlp2.o:	donlp2/donlp2-all.f
	$(FC) $(FFLAGS) -Idonlp2 -c donlp2/donlp2-all.f -o donlp2.o


# create dependencies automatically

include $(COBJS:.o=.d)

../libnumerics.a:	$(COBJS) $(FOBJS) donlp2.o
	$(AR) $(ARFLAGS) $@ $? 


clean:
	$(RM) $(OBJS)

veryclean:
	$(RM) $(OBJS)
	$(RM) $(COBJS:.o=.d)

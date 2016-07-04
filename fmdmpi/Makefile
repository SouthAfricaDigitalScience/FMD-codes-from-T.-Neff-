include ../Makefile.inc

CFLAGS := $(MPICFLAGS) -I.. -DMPI
CC = $(MPICC)

OBJLIBS = ../libfmdmpi.a
OBJS 	= Communication.o \
	Hamiltonianmpi.o gradHamiltonianmpi.o \
	MinimizerSlave.o MinimizerprojSlave.o \
	Projectionmpi.o ProjectionMultimpi.o \
	ProjectionSlave.o FormfactorSlave.o DiagonalDensityHOSlave.o \
	OneNucleonOvlapsSlave.o TwoNucleonOvlapsSlave.o


all:	$(OBJLIBS)


# create dependencies

include $(OBJS:.o=.d)

../libfmdmpi.a:	$(OBJS)
	$(AR) $(ARFLAGS) $@ $? 


clean:
	$(RM) $(OBJS)

veryclean:
	$(RM) $(OBJS) $(OBJS:.o=.d)
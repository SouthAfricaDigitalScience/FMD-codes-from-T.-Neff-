.PHONY: all doc clean veryclean install

.SUFFIXES:
.SUFFIXES: .c .f .o .mpi.o .a


SHELL = /bin/sh
ECHO = echo


CC = gcc
CLANGFLAGS = --std=gnu99 -Wall
COPTFLAGS = -O2 -march=native -ffast-math
CFLAGS = $(CLANGFLAGS) $(COPTFLAGS) -I$(NCURSES_DIR)/include -I$(OPENMPI_DIR)/include -I$(LAPACK_DIR)/include

FC = gfortran
FLANGFLAGS = -Wall
FOPTFLAGS = -O2 -march=native -ffast-math
FFLAGS = $(FLANGFLAGS) $(FOPTFLAGS)

LD = gcc
LDFLAGS = -L. -L$(LAPACK_DIR)/lib -L$(NCURSES_DIR)/lib -L$(OPENMPI_DIR)/lib -L$(GCC_DIR)/lib64
SYSLIBS = -lgfortran -lm -lz

LAPACKLIBS = -llapack -lblas

# compare with mpicc -showme:compile, mpicc -showme:link

MPIPATH = ${OPENMPI_DIR}

MPICC = mpicc
MPICFLAGS = -I$(MPIPATH)/include -I$(MPIPATH)/include/openmpi -pthread $(CFLAGS)

MPILD = mpicc
MPILDFLAGS = -L. -pthread -L$(MPI_DIR)/lib -L$(LAPACK_DIR)/lib -L$(NCURSES_DIR)/lib -L$(OPENMPI_DIR)/lib -L$(GCC_DIR)/lib64
MPILIBS = -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl


AR = ar

INSTALL = install -m 750


# implicit rules

%.mpi.o: %.c
	$(MPICC) -c -DMPI $(CPPFLAGS) $(MPICFLAGS) $< -o $@

# automatic dependencies

%.d: %.c
	$(SHELL) -ec '$(CC) -MM $(CFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

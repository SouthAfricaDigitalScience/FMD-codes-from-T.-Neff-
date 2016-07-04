SRC = ..

include $(SRC)/Makefile.inc

CFLAGS := 	$(CFLAGS) -I$(SRC)
LIBS =		-L$(SRC) -lfmd -lnumerics -lmisc $(LAPACKLIBS) $(SRC)/svn_version.o $(SYSLIBS)

OBJLIBS =	$(SRC)/libfmd.a $(SRC)/libnumerics.a $(SRC)/libmisc.a
OBJS =		printobsmeproj.o printovlmeproj.o \
		printobsmeproj-fmd-frozen.o printovlmeproj-fmd-frozen.o


all: 	printobsmeproj printovlmeproj \
	printobsmeproj-fmd-frozen printovlmeproj-fmd-frozen


printobsmeproj-fmd-frozen:	printobsmeproj-fmd-frozen.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ printobsmeproj-fmd-frozen.o $(LIBS)

printovlmeproj-fmd-frozen:	printovlmeproj-fmd-frozen.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ printovlmeproj-fmd-frozen.o $(LIBS)

printobsmeproj:	printobsmeproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ printobsmeproj.o $(LIBS)

printovlmeproj:	printovlmeproj.o $(OBJLIBS)
	$(LD) $(LDFLAGS) -o $@ printovlmeproj.o $(LIBS)


# create dependencies

include $(OBJS:.o=.d)


clean:
	$(RM) $(OBJS)

veryclean:
	$(RM) $(OBJS) $(OBJS:.o=.d)

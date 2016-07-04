include ../Makefile.inc

OBJLIBS = ../libmisc.a
OBJS 	= md5.o utils.o physics.o

all:	$(OBJLIBS)


# create dependencies automatically

include $(OBJS:.o=.d)

../libmisc.a:	$(OBJS)
	$(AR) $(ARFLAGS) $@ $? 


clean:
	$(RM) $(OBJS)

veryclean:
	$(RM) $(OBJS)
	$(RM) $(OBJS:.o=.d)
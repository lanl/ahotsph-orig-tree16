TREEHOME=..
treedir_sed=\.\.
##### Application-specific stuff goes here

ifndef prog
src= grav.c interp_cylgrid.c mac.c physics.c physics_sph.c \
	sphcofm.c sphinit.c sph.c sphprint.c sphplus.c ghosts.c com.c \
	ranlib.c wvt.c polint.c initial.c

programname=wvt
else
src=$(prog).c
fsrc=$(prog).f
programname=$(prog)
endif

treedir=$(TREEHOME)

appexcludes:=-name data

##### End of application-specific setup

include $(treedir)/Make-common/Make.$(ARCH)

include $(treedir)/Make-common/Make.generic

$(objdir)/grav$(objsuf) : grav.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c grav.c
	-@mv grav$(objsuf) $(objdir)

$(objdir)/sph$(objsuf) : sph.c
	$(CC) $(CFLAGS) $(AGGRESSIVE_OPT) -c sph.c
	-@mv sph$(objsuf) $(objdir)

# DO NOT DELETE THIS LINE -- make depend depends on it.

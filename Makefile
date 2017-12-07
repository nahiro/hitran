CC			= gcc
FC			= gfortran
LD			= g++
CFLAGS			= -g -O2
#LIBS			= -lfrtbegin -lg2c -lm -lgcc_s
LIBS			= -lg2c -lm
CERNDIR			= /cern/pro
MINUITDIR		= $(CERNDIR)/c-minuit
INCDIR			= -I$(MINUITDIR)
GCCVER			:= $(shell gcc -dumpversion)
LIBDIR			= -L/usr/lib/gcc-lib/i386-redhat-linux/$(GCCVER) -L/usr/lib
DEPEND			= $(MINUITDIR)/minuit/code/libminuit.a
DEFS			= -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"c-minuit\" -DVERSION=\"CVS-20031201\" -DHAVE_LIBM=1 -Df2cFortran=1 -DSTDC_HEADERS=1 -DTIME_WITH_SYS_TIME=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_UNISTD_H=1 -DHAVE_GETOPT_H=1 -DHAVE_GETLINE=1
ROOTDIR			= /cern/pro/root/bin
ISROOT			:= $(shell stat -c "%F" $(ROOTDIR) 2>/dev/null | tr A-Z a-z)
ifeq ($(ISROOT),directory)
ROOTCLIB0S		:= $(shell root-config --cflags)
ROOTLIBS		:= $(shell root-config --libs) -lMinuit
ROOTGLIBS		:= $(shell root-config --glibs)
endif
ROOTINC			= $(CERNDIR)/root/include
BINDIR			= ../bin
LIB0			= -lm
LIB1			= -lcfitsio
LIB2			= -ljpeg
LIB3			= -lgsl -lgslcblas
LIB4			= -lfftw3
LIB5			= -L/usr/lib/gcc/i386-redhat-linux/4.0.2 -lstdc++
OPT			= -O
WARNING			= -Wall
OBJ			= $(BINDIR)/strutil.o
OBJ1			= $(BINDIR)/strutil.o $(BINDIR)/numutil.o
OBJ2			= $(BINDIR)/novas.o $(BINDIR)/novascon.o $(BINDIR)/solsys2.o $(BINDIR)/readeph0.o $(BINDIR)/jplint.o $(BINDIR)/jplsubs.o
OBJ3			= $(BINDIR)/astro_common.o
OBJ4			= $(BINDIR)/BD_TIPS_2004.o $(BINDIR)/BD_TIPS_2008.o $(BINDIR)/BD_TIPS_2012.o $(BINDIR)/BD_TIPS_2016.o
OBJ5			= $(BINDIR)/fftvoigt.o
OBJ6			= $(BINDIR)/faddeeva.o
OBJS			= $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)
PROG			= direct_sun hitran_absorption hitran_intensity hitran_output read_hitran hitran_tips
#----------------------------------------------------
all			: clean $(PROG)
			rm -f $(OBJS)
clean			:
			rm -f $(OBJS)
hitran_absorption	: hitran_absorption.c $(OBJ1) $(OBJ4) $(OBJ6)
			$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
			$(FC) $(BINDIR)/$@.o $(OBJ1) $(OBJ4) $(OBJ6) $(LIB0) $(LIB3) $(LIB4) -o $(BINDIR)/$@ $(WARNING)
			rm -f $(BINDIR)/$@.o
hitran_intensity	: hitran_intensity.c $(OBJ1) $(OBJ4) $(OBJ6)
			$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
			$(FC) $(BINDIR)/$@.o $(OBJ1) $(OBJ4) $(OBJ6) $(LIB0) $(LIB3) -o $(BINDIR)/$@ $(WARNING)
			rm -f $(BINDIR)/$@.o
hitran_output		: hitran_output.c $(OBJ1) $(OBJ4) $(OBJ6)
			$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
			$(FC) $(BINDIR)/$@.o $(OBJ1) $(OBJ4) $(OBJ6) $(LIB0) $(LIB3) -o $(BINDIR)/$@ $(WARNING)
			rm -f $(BINDIR)/$@.o
hitran_tips		: hitran_tips.c $(OBJ) $(OBJ4) $(OBJ6)
			$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
			$(FC) $(BINDIR)/$@.o $(OBJ) $(OBJ4) $(OBJ6) $(LIB0) $(LIB3) -o $(BINDIR)/$@ $(WARNING)
			rm -f $(BINDIR)/$@.o
direct_sun		: direct_sun.c $(OBJ1) $(OBJ4) $(OBJ6)
			$(CC) -c $@.c -o $(BINDIR)/$@.o $(WARNING)
			$(FC) $(BINDIR)/$@.o $(OBJ1) $(OBJ4) $(OBJ6) $(LIB0) $(LIB3) -o $(BINDIR)/$@ $(WARNING)
			rm -f $(BINDIR)/$@.o
read_hitran		: read_hitran.c
			$(CC) $@.c -o $(BINDIR)/$@ $(LIB0) $(WARNING)
$(BINDIR)/strutil.o	: strutil.c
			$(CC) -c strutil.c   -o $@ $(WARNING)
$(BINDIR)/numutil.o	: numutil.c
			$(CC) -c numutil.c   -o $@ $(WARNING)
$(BINDIR)/tutil.o	: tutil.c
			$(CC) -c tutil.c     -o $@ $(WARNING)
$(BINDIR)/BD_TIPS_2004.o: BD_ISO_2004.for BD_MOL_2004.for BD_TIPS_2004.for Isotops_2004.cmn Molec_2004.cmn Species_2004.cmn
			$(FC) -c BD_TIPS_2004.for -o $@
$(BINDIR)/BD_TIPS_2008.o: BD_ISO_2008.for BD_MOL_2008.for BD_TIPS_2008.for Isotops_2008.cmn Molec_2008.cmn Species_2008.cmn
			$(FC) -c BD_TIPS_2008.for -o $@
$(BINDIR)/BD_TIPS_2012.o: BD_ISO_2012.for BD_MOL_2012.for BD_TIPS_2012.for Isotops_2012.cmn Molec_2012.cmn Species_2012.cmn
			$(FC) -c BD_TIPS_2012.for -o $@
$(BINDIR)/BD_TIPS_2016.o: BD_ISO_2016.for BD_MOL_2016.for BD_TIPS_2016.for Isotops_2016.cmn Molec_2016.cmn Species_2016.cmn
			$(FC) -c BD_TIPS_2016.for -ffixed-line-length-132 -o $@
$(BINDIR)/fftvoigt.o	: fftvoigt.c
			$(CC) -c fftvoigt.c  -o $@ $(WARNING)
$(BINDIR)/faddeeva.o	: faddeeva.f
			$(FC) -c faddeeva.f -o $@

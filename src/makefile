# Makefile for dns
#
#

UNAME = $(shell uname)
ARCH = $(shell uname -m)


ifeq ($(UNAME),Linux)
   #F90 = lf95
   #FFLAGS = -O -DLINUX $(MPI)
   #LDFLAGS = --staticlink 

   F90 = pgf90
   FFLAGS = -fast -DLINUX $(MPI)

   HOST = $(shell hostname)
   ifeq ($(HOST),darkstar)
      MPI = -DUSE_MPI -I/contrib/lam/include
      MPILIB =  -L/contrib/lam/lib -lpmpi -llamf77mpi -lmpi -llam 
   endif
endif
ifeq ($(UNAME),SunOS)
   F90 = f95
   #MPILIB = -L/usr/local/src/mpich/lib/solaris/ch_p4  -lmpi -lsocket -lnsl
   #MPI = -DUSE_MPI -I/usr/local/src/mpich/include
   # default -fast is -xvector=yes, but this introduces bugs in test.F90!
   FFLAGS = -fast -xvector=no -xarch=v9  -DSUNOS $(MPI)
endif
ifeq ($(UNAME),IRIX64)
   F90 = f90
   MPILIB = -lmpi 
   LDFLAGS =  -l scs 
   FFLAGS = -64 -O3 -macro_expand -DUSE_MPI  -DIRIX64  -DSGIFFT
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) 



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o fft_interface.o mpi.o
#modules/routines which use OBJS1
OBJS2 = transpose.o
#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o 
OPS    = fftops.o test.o integrals.o
FILEIO = fileio.o
OBJS3 = $(INITS) $(OPS) $(FILEIO)

MODEL_GRID  = ns3d.o bc.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_SPEC  = ns3dspectral.o bc.o $(OBJS1) $(OBJS2) $(OBJS3)



#rules
%.o: %.f90
	$(COMPILE)  $< 
%.o: %.F90
	$(COMPILE)  $< 
%.o: %.f
	$(COMPILE) -o $@ $< 

dns:  $(MODEL_SPEC) dns.o 
	$(LINK) -o dns dns.o $(MODEL_SPEC) $(LDFLAGS) $(MPILIB)

dnsgrid:  $(MODEL_GRID) dns.o 
	$(LINK) -o dns dns.o $(MODEL_GRID) $(LDFLAGS) $(MPILIB)


clean: 
	rm -f *.o *.mod dns  

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends

ifeq ($(UNAME),IRIX64)
fft_interface.o:  fft99_interface.o
	cp -f fft99_interface.o fft_interface.o
else
fft_interface.o:  fft_james_interface.o
	cp -f fft_james_interface.o fft_interface.o
endif

params.h:  
	echo "Run gridsetup.py to generage a params.h file!"
	exit 1

transpose.h:
	echo "Run gridsetup.py to generage a transpose.h file!"
	exit 1








 
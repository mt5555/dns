# Makefile for dns
#
# FFT OPTIONS:
#     USE_STKFFT
#     USE_FFT99
#     USE_SGIFFT
#
#

UNAME = $(shell uname)
ARCH = $(shell uname -m)



ifeq ($(UNAME),Linux)
   HOST = $(shell hostname)
   ifeq ($(HOST),darkstar)
      MPI = -DUSE_MPI -I/contrib/lam/include
      MPILIB =  -L/contrib/lam/lib -lpmpi -llamf77mpi -lmpi -llam 
   endif

   #F90 = lf95
   #FFLAGS = -O -DLINUX $(MPI) -DUSE_STKFFT
   #FFLAGS = -O -DLINUX $(MPI) -DUSE_FFT99
   #LDFLAGS = --staticlink  -Lstk -lstk
   #CC = gcc

   F90 = pgf90
   #FFLAGS = -fast -DLINUX $(MPI)  -DUSE_FFT99 
   FFLAGS = -g -fast -DLINUX $(MPI)  -DUSE_STKFFT
   LDFLAGS = -Lstk -lstk
endif
ifeq ($(UNAME),SunOS)
   F90 = f95
   #MPILIB = -L/usr/local/src/mpich/lib/solaris/ch_p4  -lmpi -lsocket -lnsl
   #MPI = -DUSE_MPI -I/usr/local/src/mpich/include
   # default -fast is -xvector=yes, but this introduces bugs in test.F90!
   FFLAGS = -fast -xvector=no -xarch=v9  -DSUNOS $(MPI) -DUSE_STKFFT
   LDFLAGS = -Lstk -lstk
   CC = cc -xarch=v9 -c
endif
ifeq ($(UNAME),IRIX64)
   F90 = f90
   MPI = -DUSE_MPI
   MPILIB = -lmpi 
   LDFLAGS =  -l scs -Lstk -lstk
   FFLAGS = -64 -O3 -macro_expand  -DIRIX64  -DUSE_SGIFFT $(MPI)
   CC = cc -64
endif
ifeq ($(UNAME),OSF1)
   F90 = f90
   #MPI = -DUSE_MPI -I/usr/opt/MPI190/include
   #MPILIB = -L/usr/opt/MPI190/lib -lfmpi -lmpi 
   LDFLAGS =  $(MPILIB) -Lstk -lstk -lcxml

   #FFLAGS = -O3  -DOSF1 $(MPI) -DUSE_FFT99
   #FFLAGS = -O3  -DOSF1 $(MPI) -DUSE_STKFFT  
   FFLAGS = -O3  -DOSF1 $(MPI) -DUSE_CPQFFT  
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) 



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o  mpi.o fft_interface.o
#modules/routines which use OBJS1
OBJS2 = transpose.o 
#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o 
OPS    = fftops.o test.o integrals.o
FILEIO = fileio.o fileioc.o
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
%.o: %.c
	$(CC) -c -o $@ $< 

dns:  $(MODEL_SPEC) dns.o 
	$(LINK) -o dns dns.o $(MODEL_SPEC) $(LDFLAGS) $(MPILIB)

dnsgrid:  $(MODEL_GRID) dns.o 
	$(LINK) -o dns dns.o $(MODEL_GRID) $(LDFLAGS) $(MPILIB)


clean: 
	rm -f *.o *.mod dns  

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends


params.h:  
	echo "Run gridsetup.py to generage a params.h file!"
	exit 1

transpose.h:
	echo "Run gridsetup.py to generage a transpose.h file!"
	exit 1


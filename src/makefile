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
      MPI = 
      MPILIB =  
   else
      MPI = -DUSE_MPI -I/contrib/lam/include
      MPILIB =  -L/contrib/lam/lib -lpmpi -llamf77mpi -lmpi -llam 
   endif

   ## Lahey
   #F90 = lf95
   ##FFLAGS = -g -O -DLINUX $(MPI) -DUSE_STKFFT
   #FFLAGS = -g -O -DLINUX $(MPI) -DUSE_FFT99
   #LDFLAGS = --staticlink  #-Lstk -lstk
   #CC = gcc

   ## PGI
   #F90 = pgf90
   #CC=pgcc
   #FFLAGS = -Mrecursive -fast -DLINUX $(MPI)  -DUSE_FFT99 
   ##LDFLAGS = -Lstk -lstk


   #INTEL	
   F90 = ifc 
   CC=gcc
   FFLAGS =  -auto -w90  -Vaxlib -O1 -tpp7  -DLINUX $(MPI)  -DUSE_FFT99 
   #LDFLAGS = -Lstk -lstk

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
   LDFLAGS =  -l scs #-Lstk -lstk
   FFLAGS = -64 -Ofast -macro_expand  -DIRIX64  -DNDEBUG -DUSE_SGIFFT $(MPI)
   CC = cc -64
endif
ifeq ($(UNAME),OSF1)
   F90 = f90
   ifeq ($(HOST),mahi)
      MPI = 
      MPILIB =  
   else
   ifeq ($(HOST),truchas)
      MPI = -DUSE_MPI -I/usr/opt/MPI190/include
      MPILIB = -L/usr/opt/MPI190/lib -lfmpi -lmpi 
   else	
      MPI = -DUSE_MPI -I/usr/opt/MPI/include
      MPILIB = -L/usr/opt/MPI/lib -lfmpi -lmpi -lelan
   endif
   endif

   LDFLAGS =  $(MPILIB) -lcxml #-Lstk -lstk 

   #FFLAGS = -recursive -fast -DOSF1 $(MPI) -DUSE_FFT99
   #FFLAGS = -recursive -fast -DOSF1  $(MPI) -DUSE_STKFFT  
   FFLAGS = -fast  -DOSF1 $(MPI) -DUSE_CPQFFT  
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) 



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o  mpi.o fft_interface.o 
#modules/routines which use OBJS1
OBJS2 = transpose.o structf.o sforcing.o solve.o
#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o 
OPS    = fftops.o test.o integrals.o
FILEIO = fileio.o fileioc.o
OBJS3 = $(INITS) $(OPS) $(FILEIO)

MODEL_GRID       = ns_grid.o  dns.o bc.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_IMPULSE    = ns_impulse.o  dns.o bc.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_SPEC_SLOW  = ns_decomp3.o dns.o bc.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_SPEC       = ns.o  dns.o bc.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_ALPHA      = ns_alpha.o dns.o  bc.o $(OBJS1) $(OBJS2) $(OBJS3)


#rules
%.o: %.f90
	$(COMPILE)  $< 
%.o: %.F90
	$(COMPILE)  $< 
%.o: %.f
	$(COMPILE) -o $@ $< 
%.o: %.c
	$(CC) -c -o $@ $< 

dns:  $(MODEL_SPEC) 
	$(LINK) -o dns $(MODEL_SPEC) $(LDFLAGS) $(MPILIB)

dnsa:  $(MODEL_ALPHA) 
	$(LINK) -o dnsa $(MODEL_ALPHA) $(LDFLAGS) $(MPILIB)

dnsgrid:  $(MODEL_GRID) 
	$(LINK) -o dnsgrid $(MODEL_GRID) $(LDFLAGS) $(MPILIB)

dnsi:  $(MODEL_IMPULSE) 
	$(LINK) -o dnsi $(MODEL_IMPULSE) $(LDFLAGS) $(MPILIB)

dnsslow:  $(MODEL_SPEC_SLOW) 
	$(LINK) -o dns $(MODEL_SPEC_SLOW) $(LDFLAGS) $(MPILIB)


clean: 
	rm -f *.d *.o *.mod dns  

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends

ns_alpha.o: ns.F90 ns.o
	$(COMPILE) -DALPHA_MODEL -o ns_alpha.o ns.F90

params.h:  
	echo "Run gridsetup.py to generage a params.h file!"
	exit 1

transpose.h:
	echo "Run gridsetup.py to generage a transpose.h file!"
	exit 1


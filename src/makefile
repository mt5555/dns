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
   ##F90 = /usr/local/lf9560/bin/lf95
   #F90 = lf95
   #FFLAGS = -g -O -DLINUX $(MPI) -DUSE_FFT99
   ##FFLAGS =  -O -DLINUX $(MPI) -DUSE_FFT99
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
   #LSFLIB = -L /lsf/lib64 -lbat -llsf 
   #LSFFLAGS = -I/usr/include -I/lsf/share/include -DUSE_LSFTIME

   LDFLAGS =  -l scs $(LSFLIB)  #-Lstk -lstk
   FFLAGS = -64 -Ofast -macro_expand  -DIRIX64  -DNDEBUG -DUSE_SGIFFT $(MPI)
   CPPFLAGS = 
   CC = cc -64 $(LSFFLAGS)
endif
ifeq ($(UNAME),OSF1)
   F90 = f90
   ifeq ($(HOST),mahi)
      MPI = 
      MPILIB =  
      LSFLIB = 
      LSFFLAGS = 
   else
   ifeq ($(HOST),truchas)
      MPI = -DUSE_MPI -I/usr/opt/MPI190/include
      MPILIB = -L/usr/opt/MPI190/lib -lfmpi -lmpi 
   else	
      MPI = -DUSE_MPI -I/usr/opt/MPI/include
      MPILIB = -L/usr/opt/MPI/lib -lfmpi -lmpi -lelan
      LSFLIB = -L/lsf/lib -lbat -llsf 
      LSFFLAGS = -I/usr/include -I/lsf/include -DUSE_LSFTIME
   endif
   endif


   LDFLAGS =  $(MPILIB) -lcxml $(LSFLIB)  #-Lstk -lstk 

   #FFLAGS = -recursive -fast -DOSF1 $(MPI) -DUSE_FFT99
   #FFLAGS = -recursive -fast -DOSF1  $(MPI) -DUSE_STKFFT  
   FFLAGS = -fast  -DOSF1 $(MPI) -DNDEBUG -DUSE_CPQFFT  
   CC = cc $(LSFFLAGS)
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) 



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o  mpi.o fft_interface.o lsftime.o
#modules/routines which use OBJS1
OBJS2 = transpose.o structf.o sforcing.o solve.o ghost.o bc.o
#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o 
OPS    = ghostops.o fftops.o test.o integrals.o
FILEIO = fileio.o fileioc.o
OBJS3 = $(INITS) $(OPS) $(FILEIO)

MODEL_GRID       = ns_grid.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_GHOST      = ns_ghost.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_PSIVOR     = ns_psivor.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_IMPULSE    = ns_impulse.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_SPEC_SLOW  = ns_decomp3.o dns.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_SPEC       = ns.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_ALPHA      = ns_alpha.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3)
MODEL_SHALLOW    = shallow.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3)


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

dnsghost:  $(MODEL_GHOST) 
	$(LINK) -o dnsghost $(MODEL_GHOST) $(LDFLAGS) $(MPILIB)

dnsvor:  $(MODEL_PSIVOR) 
	$(LINK) -o dnsvor $(MODEL_PSIVOR) $(LDFLAGS) $(MPILIB)

dnsi:  $(MODEL_IMPULSE) 
	$(LINK) -o dnsi $(MODEL_IMPULSE) $(LDFLAGS) $(MPILIB)

dnss:  $(MODEL_SHALLOW) 
	$(LINK) -o dnss $(MODEL_SHALLOW) $(LDFLAGS) $(MPILIB)

dnsslow:  $(MODEL_SPEC_SLOW) 
	$(LINK) -o dns $(MODEL_SPEC_SLOW) $(LDFLAGS) $(MPILIB)

splitspec:  splitspec.o fileioc.o
	$(LINK) -o splitspec splitspec.o fileioc.o

analysis:  analysis.o $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis analysis.o $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)



clean: 
	rm -f temp* *.d *.o *.mod dns dnss dnsa dnsgrid dnsghost dnsi dnsvor

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends

ns_alpha.o: ns.F90 ns.o
	$(COMPILE) -DALPHA_MODEL -o ns_alpha.o ns.F90

params.h:  
	echo "Run gridsetup.py to generage a params.h file!"
	exit 1



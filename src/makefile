# Makefile for dns
#
#

UNAME = $(shell uname)
ARCH = $(shell uname -m)

ifeq ($(UNAME),Linux)
   F90 = lf95
   LDFLAGS = --staticlink
   MPI = -DUSE_MPI -I/contrib/lam/include
   MPILIB =  -L/contrib/lam/lib -lpmpi -llamf77mpi -lmpi -llam 
   FFLAGS = -O -DLINUX $(MPI)
endif
ifeq ($(UNAME),SunOS)
   F90 = f95
   MPILIB = -L/usr/local/src/mpich/lib/solaris/ch_p4  -lmpi -lsocket -lnsl
   MPI = -DUSE_MPI -I/usr/local/src/mpich/include
   LDFLAGS = 
   # default -fast is -xvector=yes, but this introduces bugs in test.F90!
   FFLAGS = -fast -xvector=no  -DSUNOS $(MPI)
endif
ifeq ($(UNAME),IRIX64)
   F90 = f90
   MPILIB = -lmpi 
   LDFLAGS = 
   FFLAGS = -g -64 -O1 -macro_expand -DUSE_MPI  -DIRIX64
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) $(LDFLAGS)



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o fft_interface.o mpi.o
#modules/routines which use OBJS1
OBJS2 = transpose.o
#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o 
OPS    = fftops.o test.o integrals.o
MODEL  = ns3d.o bc.o
FILEIO = fileio.o
OBJS3 = $(MODEL) $(INITS) $(OPS) $(FILEIO)



#rules
%.o: %.f90
	$(COMPILE)  $< 
%.o: %.F90
	$(COMPILE)  $< 
%.o: %.f
	$(COMPILE) -o $@ $< 

dnslink = $(OBJS1) $(OBJS2) $(OBJS3)
dns:  $(dnslink) dns.o
	$(LINK) -o dns dns.o $(dnslink) $(MPILIB)


clean: 
	rm -f *.o *.mod dns  

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends

fft_interface.o:  fft99_interface.o
	cp -f fft99_interface.o fft_interface.o
#fft_interface.o:  compaq_fft_interface.o
#	cp -f fft99_interface.o fft_interface.o

user_params.h:  
	echo "Run gridsetup.py to generage a user_params.h file!"
	exit 1








 
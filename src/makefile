# Makefile for dns
#
#

UNAME = $(shell uname)
ARCH = $(shell uname -m)

ifeq ($(UNAME),Linux)
   F90 = lf95
   LDFLAGS = --staticlink
   FFLAGS = -O
   MPIF90 = mpi.F90
endif
ifeq ($(UNAME),SunOS)
   F90 = f95
   MPILIB = -L/usr/local/src/mpich/lib/solaris/ch_p4  -lmpi -lsocket -lnsl
   MPIINC = -DMPI 
   MPIF90 = /usr/local/src/mpich/include/mpif.f90
   LDFLAGS = 
   FFLAGS = -O -DMPI
endif
ifeq ($(UNAME),IRIX64)
   F90 = f90
   MPILIB = -lmpi 
   LDFLAGS = 
   FFLAGS = -g -O1 -cpp -DMPI 
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) $(LDFLAGS)



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o fft_interface.o mpi.o
#modules/routines which use OBJS1
OBJS2 = transform.o
#modules/routines which use OBJS1,2
INITS  = init_model.o init_mpi.o init_data.o init_grid.o 
OPS    = fftops.o test.o
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


#mpi.o: $(MPIF90)
#	$(COMPILE) -o mpi.o $(MPIF90)
fft_interface.o:  fft99_interface.F90 macros.h
	$(COMPILE) -o fft_interface.o fft99_interface.F90
#fft_interface.o:  compaq_fft_interface.F90 params.o
#	$(COMPILE) compaq_fft_interface.F90




 
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
HOST = $(shell hostname)


ifeq ($(UNAME),Linux)
   MPI = 
   MPILIB = 
   ifeq ($(HOST),autrey)
      MPI = -DUSE_MPI -I/contrib/lam/include
      MPILIB =  -L/contrib/lam/lib -lpmpi -llamf77mpi -lmpi -llam 
   endif
   ifeq ($(HOST),shankara1)
      MPI = -DNDEBUG  -DMPI_HAS_REAL8 -DUSE_MPI -I/opt/scali/include
      MPILIB =  -L/opt/scali/lib -lfmpi -lmpi
   endif
   #VENDOR=LAHEY
   #VENDOR=PGI
   VENDOR=INTEL

   ifeq ($(VENDOR),LAHEY)
   #F90 = /usr/local/lf9560/bin/lf95
   F90 = lf95
   FFLAGS = -g  -DLINUX $(MPI) -DUSE_FFT99
   #FFLAGS =  -tpp -O -DLINUX $(MPI) -DUSE_FFT99
   LDFLAGS = --staticlink  #-Lstk -lstk
   CC = gcc
   endif

   ifeq ($(VENDOR),PGI)
   #F90 = pgf90
   CC=pgcc
   #FFLAGS = -Mrecursive -fast -DLINUX $(MPI)  -DUSE_FFT99 
   FFLAGS = -Mrecursive -g -DLINUX $(MPI)  -DUSE_FFT99 
   #LDFLAGS = -Lstk -lstk
   endif 


   ifeq ($(VENDOR),INTEL)
   F90 = ifc 
   CC=gcc
   FFLAGS =  -auto -w90  -Vaxlib -O1 -tpp7  -DLINUX $(MPI)  -DUSE_FFT99 
   #FFLAGS =  -w90  -Vaxlib -O1 -tpp7  -DLINUX $(MPI)  -DUSE_FFT99 
   #LDFLAGS = -Lstk -lstk
   endif

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
   LSFLIB = -L /lsf/lib64 -lbat -llsf 
   LSFFLAGS = -I/usr/include -I/lsf/share/include -DUSE_LSFTIME

   ifeq ($(HOST),guyot.acl.lanl.org)
      LSFLIB = 
      LSFFLAGS =	
   endif

   # OPENMP flag: -mp 
   LDFLAGS =  -l scs $(LSFLIB)  #-Lstk -lstk
   # When SGI adds padding to data in params.h, MPI breaks:  
   FFLAGS = -64 -Ofast -OPT:reorg_common=OFF -macro_expand  -DIRIX64  -DNDEBUG -DUSE_SGIFFT $(MPI)
   CPPFLAGS = 
   CC = cc -64 $(LSFFLAGS)
endif
ifeq ($(UNAME),OSF1)
   F90 = f90
   # defaults for asci Q: MPI paths loaded via modules
   # to enable openMP  -omp ?? 
   MPI = -DUSE_MPI 
   MPILIB = -lfmpi -lmpi -lelan
   #LSFLIB = -L/lsf/lib -lbat -llsf 
   #LSFFLAGS = -I/usr/include -I/lsf/include -DUSE_LSFTIME

   ifeq ($(HOST),milkyway.lanl.gov)
      MPI = -DUSE_MPI -I/usr/opt/MPI196/include
      MPILIB = -L/usr/opt/MPI196/lib -lfmpi -lmpi 
      #MPI = 
      #MPILIB =
      LSFLIB = 
      LSFFLAGS = 
   endif

   LDFLAGS =  $(MPILIB) -lcxml $(LSFLIB)  #-Lstk -lstk 

   #flag for OPENMP:  -omp
   #FFLAGS = -automatic -fast -DOSF1 $(MPI) -DUSE_FFT99
   #FFLAGS = -automatic -fast -DOSF1  $(MPI) -DUSE_STKFFT  
   # shallow water model cannot use -automatic
   #FFLAGS =  -fast -automatic  -DOSF1 $(MPI) -DNDEBUG -DUSE_CPQFFT  
   FFLAGS =  -fast -automatic  -DOSF1 $(MPI) -DNDEBUG -DUSE_CPQFFT  
   CC = cc $(LSFFLAGS)
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) 



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o  mpi.o fft_interface.o lsftime.o
#modules/routines which use OBJS1
OBJS2 = transpose.o  solve.o ghost.o bc.o  sforcing.o spectrum.o ellipse.o 

#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o cases2v.o cases3v.o 
OPS    = ghostops.o fftops.o test.o integrals.o timecontrol.o tracers.o 
FILEIO = fileio.o fileioc.o
OBJS3 = $(INITS) $(OPS) $(FILEIO)

#diagnostics for 3d turbulence runs:
TURB_DIAG = structf.o turb_diag.o isoave.o

#diagnostics for VXPAIR runs:
PSIVOR_DIAG = psivor_diag.o  

#diagnostics for shallow water:
SHALLOW_DIAG =  shallow_diag.o 


# n_var>=3 models:  (u,v,w) or (u,v,h)
# but they can still be run on 2D problems with nz=1
MODEL_GRID       = ns_grid.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_IBMGRID    = ns_ibm.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_GHOST      = ns_ghost.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_IMPULSE    = ns_impulse.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC_SLOW  = ns_decomp3.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC       = ns.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_ALPHA      = ns_alpha.o dns.o   $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SHALLOW    = shallow.o dns.o   $(OBJS1) $(OBJS2) $(OBJS3)   $(SHALLOW_DIAG)

# n_var>=2 models  (u,v) and (w,psi)
MODEL_PSIVOR     = ns_psivor.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3) $(PSIVOR_DIAG)



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

dnsibm:  $(MODEL_IBMGRID) 
	$(LINK) -o dnsibm $(MODEL_IBMGRID) $(LDFLAGS) $(MPILIB)

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
analysis_rms:  analysis_rms.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_rms analysis_rms.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)
analysis_isoave:  analysis_isoave.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_isoave analysis_isoave.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)



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



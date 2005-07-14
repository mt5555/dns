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
# use this if compiler needs additional options to compile fixed format code:
F77FIXED=


ifeq ($(UNAME),Linux)
   #VENDOR=LAHEY
   VENDOR=PGI
   #VENDOR=INTEL
   MPI = 
   MPILIB = 
   ifeq ($(HOST),mauve.lanl.gov)
      VENDOR=INTEL8
      MPI = -DNDEBUG -DMPI_HAS_REAL8 -DUSE_MPI 
      MPILIB =  -lmpi
   endif
	ifeq ($(HOST),cnlsrocks)
      VENDOR=PGI
      MPI = -DNDEBUG -DMPI_HAS_REAL8 -DUSE_MPI -I/usr/local/pkg/mpich-gm/mpich-gm-pgi/mpich-gm-1.2.6.14/include
      MPILIB =  -L/usr/local/pkg/mpich-gm/mpich-gm-pgi/mpich-gm-1.2.6.14/lib -lfmpich  -lmpich  -L/usr/gm/lib/ -lgm
                 #-lpthread
   endif


   ifeq ($(HOST),s854018)
      VENDOR=INTEL8
      #MPI = -DUSE_MPI -DMPI_UNDERSCORE
      #MPILIB =   -llamf77mpi -lmpi -llam -lpthread
   endif
   ifeq ($(HOST),kiwi1)
      VENDOR=INTEL8
      MPI = -DUSE_MPI  -I/usr/include/lam  -DMPI_UNDERSCORE
      MPILIB =   -L/usr/lib/lam/lib -llamf77mpi -lmpi -llam 
      #MPI = -DUSE_MPI -DUSE_MPI_IO -DMPI_HAS_REAL8 -I$(HOME)/lampi/sulaco/include
      #MPILIB =  -Xlinker -rpath $(HOME)/lampi/sulaco/lib -L$(HOME)/lampi/sulaco/lib -lmpi
   endif
   ifeq ($(HOST),dosadi)
      VENDOR=INTEL8
      MPI = -DUSE_MPI -DMPI_UNDERSCORE
      MPILIB =   -llamf77mpi -lmpi -llam 
   endif
   ifeq ($(HOST),reddish)
      VENDOR=PGI
      MPI = -I/usr/pgi/linux86-64/5.2/include -DUSE_MPI -DMPI_HAS_REAL8
      MPILIB=  -L/usr/pgi/linux86-64/5.2/lib -lfmpich -lmpich 
      #MPI = -I/usr/lib/mpi/mpi_pgi/include -DUSE_MPI -DMPI_HAS_REAL8 -Dmpi_wtime=pmpi_wtime
      #MPILIB = -L/usr/lib/mpi/mpi_pgi/lib  -lfmpi -lmpi 
   endif
   ifeq ($(HOST),brain)
      VENDOR=INTEL
      MPI =  -DNDEBUG -DMPI_HAS_REAL8 -DUSE_MPI -I/opt/scali/include
      MPILIB =  -L/opt/scali/lib -lfmpi -lmpi
   endif
   ifeq ($(HOST),shankara12)
      VENDOR=INTEL
      MPI = -DUSE_MPI -I/users/contrib/lam/include
      MPILIB =   -L/users/contrib/lam/lib -llamf77mpi -lmpi -llam -lpthread
   endif
   ifeq ($(HOST),pinkish.lanl.gov)
      #MPI =  -DNDEBUG -DMPI_HAS_REAL8 -DUSE_MPI -I/usr/mpich-gm.icc/include
      #MPILIB =-L/usr/mpich-gm.icc/lib  -lfmpich  -lmpich  -L/usr/ib -lgm -lpthread
      MPI =  -DUSE_MPI -DMPI_HAS_REAL8 -DUSE_LAMPI -I/home/taylorm/lampi/pinkish/include
      MPILIB = -L /home/taylorm/lampi/pinkish/lib -lmpi
   endif


   ifeq ($(VENDOR),LAHEY)
   #F90 = /usr/local/lf9562/bin/lf95
   F90 = lf95
   FFLAGS = -g  -DLINUX $(MPI) -DUSE_FFT99
   #FFLAGS =  -tpp -O -DLINUX $(MPI) -DUSE_FFT99
   LDFLAGS = --staticlink  #-Lstk -lstk
   CC = gcc
   endif

   ifeq ($(VENDOR),PGI)
   F90 = pgf90
   CC=pgcc
   #FFLAGS = -Bstatic -Mrecursive -fast -DLINUX $(MPI)  -DUSE_FFT99 
   FFLAGS = -Mrecursive -g -DLINUX $(MPI)  -DUSE_FFT99 
   #LDFLAGS = -Lstk -lstk
   endif 


   ifeq ($(VENDOR),INTEL)
   F90 = ifc
   CC=gcc
   FFLAGS =  -auto -w90  -Vaxlib  -O1 -tpp7  -DLINUX $(MPI)  -DUSE_FFT99 
   endif

   ifeq ($(VENDOR),INTEL8)
   F90 = ifort
   CC=gcc
   #FFLAGS =  -auto -w90  -Vaxlib  -O1 -tpp7  -DLINUX $(MPI)  -DUSE_FFT99 
   FFLAGS =  -auto -w90  -Vaxlib  -O1   -DLINUX $(MPI)  -DUSE_FFT99 
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
ifeq ($(UNAME),AIX)
   F90 = mpxlf95
   F77FIXED = -qfixed
   FFLAGS = -WF,"-DAIX,-Dabort=punt,-DUSE_MPI,-DUSE_FFT99,-DNDEBUG"  -qsuffix=cpp=F90  -g -O3 -qmaxmem=-1 -qstrict 
   CC = cc -DAIX
endif
ifeq ($(UNAME),OSF1)
   F90 = f90

   #defaults for QA
   UDM = 
   UDMLIB = 
#      UDM = -I$(HOME)/udm/udm_Nov_19/include -DHAVE_UDM 
#      UDMLIB = -L$(HOME)/udm/udm_Nov_19/lib \
#         -L/usr/projects/ups/v-02-03-02_gamma_io/ALPHA64_mpi/lib \
#         -ludm -lhdf5 -lmpio -lmpi -lelan -lm \
#        -lrt -free -lcxxstd -lcxx -lexc -lc
   LDFLAGS =  -lcxml $(LSFLIB) $(UDMLIB) 

   # check for MPI module loaded:
   MOD_MPIR12 = $(shell echo $$MPI_LD_FLAGS | grep MPI)
   MOD_lampi = $(shell echo $$MPI_LD_FLAGS | grep lampi)
   ifneq ($(MOD_MPIR12),)
      # R12
      MPI =  -DUSE_MPI -DUSE_MPI_IO $(MPI_COMPILE_FLAGS)
      MPILIB =  $(MPI_LD_FLAGS) -lfmpi -lmpio -lmpi  -lelan
   else
      ifneq ($(MOD_lampi),)
         # lampi module
         MPI =  -DUSE_MPI -DUSE_MPI_IO $(MPI_COMPILE_FLAGS)
         MPILIB =  $(MPI_LD_FLAGS) -lmpi -lelan
      else
         # lampi
         MPI =  -DUSE_MPI -DUSE_MPI_IO -I/users/taylorm/lampi/qsc-stable/include
         MPILIB =  -L/users/taylorm/lampi/qsc-stable/lib -lmpi  -lelan
      endif	
   endif


   ifeq ($(HOST),milkyway.lanl.gov)
      MPI =  -DUSE_MPI_IO  -DUSE_MPI -I/usr/opt/MPI196/include
      MPILIB = -L/usr/opt/MPI196/lib -lfmpi -lmpi 
      UDM =  
      UDMLIB = 
      LDFLAGS =  -lcxml $(LSFLIB) $(UDMLIB) #-Lstk -lstk 
   endif

   # note: shallow water model cannot use -automatic
   FFLAGS =  -fast -automatic  -DOSF1 $(MPI) -DNDEBUG -DUSE_CPQFFT  $(UDM)
   CC = cc $(LSFFLAGS)
endif

COMPILE=$(F90) -c $(FFLAGS)
COMPILE77=$(F90) $(F77FIXED) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) 



#independent modules/routines (dont use any other modules):
OBJS1 = params.o util.o  fft99.o  mpi.o fft_interface.o lsftime.o
#modules/routines which use OBJS1
OBJS2 = transpose.o  solve.o ghost.o bc.o  sforcing.o spectrum.o ellipse.o 

#modules/routines which use OBJS1,2
INITS  = init_mpi.o init_data.o init_grid.o cases2v.o cases3v.o 
OPS    = ghostops.o fftops.o test.o integrals.o timecontrol.o tracers.o particles.o 
FILEIO = fileio.o fileioc.o udm_interface.o
OBJS3 = $(INITS) $(OPS) $(FILEIO)

#diagnostics for 3d turbulence runs:
TURB_DIAG = pdf.o turb_diag.o isoave.o subcubes.o

#diagnostics for compressible DNS runs:
CDNS_DIAG = pdf.o cdns_diag.F90

#diagnostics for VXPAIR runs:
PSIVOR_DIAG = psivor_diag.o  

#diagnostics for shallow water:
SHALLOW_DIAG =  shallow_diag.o 


# n_var>=3 models:  (u,v,w) or (u,v,h)
# but they can still be run on 2D problems with nz=1
MODEL_GRID       = ns_grid.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_IBMGRID    = ns_ibm.o  ibm.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_GHOST      = ns_ghost.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_CDNS       = cns.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(CDNS_DIAG)
MODEL_IMPULSE    = ns_impulse.o  dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC_SLOW  = ns_decomp3.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC       = ns.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC_SCALE = ns_aspect.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_ALPHA      = ns_alpha.o dns.o   $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_LERAY      = ns_leray.o dns.o   $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SHALLOW    = shallow.o dns.o   $(OBJS1) $(OBJS2) $(OBJS3)   $(SHALLOW_DIAG)

# n_var>=2 models  (u,v) and (w,psi)
MODEL_PSIVOR     = ns_psivor.o  dns.o $(OBJS1) $(OBJS2) $(OBJS3) $(PSIVOR_DIAG)



#rules
%.o: %.f90
	$(COMPILE)  $< 
%.o: %.F90
	$(COMPILE)  $< 
%.o: %.f
	$(COMPILE77) -o $@ $< 
%.o: %.c
	$(CC) -c -o $@ $< 

dns:  $(MODEL_SPEC) 
	$(LINK) -o dns $(MODEL_SPEC) $(LDFLAGS) $(MPILIB)
dnsz: $(MODEL_SPEC_SCALE) 
	$(LINK) -o dnsz $(MODEL_SPEC_SCALE) $(LDFLAGS) $(MPILIB)

dnsa:  $(MODEL_ALPHA) 
	$(LINK) -o dnsa $(MODEL_ALPHA) $(LDFLAGS) $(MPILIB)
dnsl:  $(MODEL_LERAY) 
	$(LINK) -o dnsl $(MODEL_LERAY) $(LDFLAGS) $(MPILIB)

dnsibm:  $(MODEL_IBMGRID) 
	$(LINK) -o dnsibm $(MODEL_IBMGRID) $(LDFLAGS) $(MPILIB)

dnsgrid:  $(MODEL_GRID) 
	$(LINK) -o dnsgrid $(MODEL_GRID) $(LDFLAGS) $(MPILIB)

dnsghost:  $(MODEL_GHOST) 
	$(LINK) -o dnsghost $(MODEL_GHOST) $(LDFLAGS) $(MPILIB)

cdns:  $(MODEL_CDNS) 
	$(LINK) -o cdns $(MODEL_CDNS) $(LDFLAGS) $(MPILIB)

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

conv:  conv.o $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
	$(LINK) -o conv conv.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

convert:  convert.o $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
	$(LINK) -o convert convert.o $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG) $(LDFLAGS) $(MPILIB)

analysis:  analysis.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis analysis.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis2:  analysis2.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis2 analysis2.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_rms:  analysis_rms.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_rms analysis_rms.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_isoave:  analysis_isoave.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_isoave analysis_isoave.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_sk: analysis_sk.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_sk analysis_sk.o  $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)


clean: 
	rm -f temp* *.d *.o *.mod dns dnss dnsa dnsgrid dnsghost dnsi dnsvor \
                            cdns dnsibm conv analysis_rms     \
                            analysis_isoave analysis analysis_sk splitspec \
			convert	

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends

ns_alpha.o: ns.F90 ns.o
	$(COMPILE) -DALPHA_MODEL -o ns_alpha.o ns.F90
ns_leray.o: ns.F90 ns.o
	$(COMPILE) -DALPHA_MODEL -DLERAY_MODEL -o ns_leray.o ns.F90

params.h:  
	echo "Run gridsetup.py to generage a params.h file!"
	exit 1



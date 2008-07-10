# Makefile for dns
#
# FFT OPTIONS:
#     USE_FFT99  (default)
#     USE_FFTW
#     USE_FFTWGURU
#     USE_STKFFT
#     USE_SGIFFT
#
#

UNAME = $(shell uname)
ARCH = $(shell uname -m)
HOST = $(shell hostname)
# use this if compiler needs additional options to compile fixed format code:
F77FIXED=
MPI = 
MPILIB = 

USE_ALPHA = -DALPHA_MODEL
USE_LERAY = -DALPHA_MODEL -DLERAY_MODEL

ifeq ($(UNAME),Linux)
   VENDOR=LAHEY
   #VENDOR=PGI
   #VENDOR=PGIMPI
   #VENDOR=INTEL
   #VENDOR=INTELMPI
   #VENDOR=GNU
   #VENDOR=GNUMPI

   FFT = -DUSE_FFT99
   FFTLIB =

   ifeq ($(ARCH),ppc64)
      VENDOR=IBM
      F90 = mpixlf90
      CC=mpixlc -DF_NO_UNDERSCORE
      FFLAGS = -WF,"-DNDEBUG,-DUSE_MPI,-DMPI_HAS_REAL8,-DUSE_FFT99,-Dflush=flush_,-DF_NO_UNDERSCORE"  -qsuffix=cpp=F90   -O2 -qmaxmem=-1 -qstrict -I.
      LDFLAGS = 
      F77FIXED = -qfixed 
   endif

   ifeq ($(ARCH),x86_64)
      # PGI or Intel:
      ifeq ($(PGI),)
         VENDOR=INTELMPI
         MPI = -DNDEBUG -DUSE_MPI  -DMPI_HAS_REAL8 $(MPI_COMPILE_FLAGS)
         MPILIB =  $(MPI_LD_FLAGS) -lmpi 	
      else
         # PGI is defined:
         VENDOR=PGIMPI
         # needed if >2GB per process:   -DMALLOC_Q
         MPI=  -DNDEBUG -DUSE_MPI -DMPI_HAS_REAL8 -DUSE_MPI_IO 
         #FFT = -I/home/mataylo/fftw/include -DUSE_FFTWGURU
         #FFTLIB = -L/home/mataylo/fftw/lib  -lfftw3
      endif 
   
   endif

   ifeq ($(ARCH),ia64)
      VENDOR=INTEL
      MPI = -DNDEBUG -DUSE_MPI  -DMPI_HAS_REAL8
      MPILIB = -lmpi
   endif

    ifeq ($(HOST),dosadi)
      VENDOR=GNUMPI
      MPI = -DUSE_MPI -DMPI_HAS_REAL8 -DUSE_MPI_IO
      #FFT = -DUSE_FFTWGURU
      #FFTLIB =  -lfftw3
      #FFT = -DUSE_STK
      #FFTLIB = -Lstk -lstk
   endif

   ifeq ($(HOST),serrano1)
      VENDOR=INTELMPI
      MPI = -DNDEBUG -DUSE_MPI -DMPI_HAS_REAL8
      MPILIB = 
   endif

   ifeq ($(HOST),mauve.lanl.gov)
      VENDOR=INTEL
      MPI = -DNDEBUG -DMPI_HAS_REAL8 -DUSE_MPI 
      MPILIB =  -lmpi
   endif

   ifeq ($(HOST),nostromo)
      #VENDOR=GNU
      VENDOR=GNUMPI
      MPI = -DNDEBUG -DUSE_MPI -DMPI_HAS_REAL8
   endif
   ifeq ($(HOST),mora)
      VENDOR=INTEL
   endif

   ifeq ($(HOST),sulaco)
      VENDOR=INTEL
#      VENDOR=PGI
## Note: LAM 7.1 added MPI_REAL8, but it is broken for collectives
## so we need 
#      MPI =  -DUSE_MPI -DUSE_MPI_IO -DMPI_UNDERSCORE -DMPI_HAS_REAL8 -DMPI_REAL8_LAM
#      MPILIB =   -llammpio -llamf77mpi -lmpi -llam -lutil -lpthread
   endif

   ifeq ($(HOST),pfe1.lanl.gov)
      VENDOR=INTEL

      MPI = -DNDEBUG -DUSE_MPI -DMPI_HAS_REAL8 $(MPI_COMPILE_FLAGS)
      MPILIB = $(MPI_LD_FLAGS) -lfmpich -lmpich -lgm -lpthread
   endif

   ifeq ($(HOST),reddish)
      VENDOR=PGI
      MPI = -I/usr/pgi/linux86-64/5.2/include -DUSE_MPI -DMPI_HAS_REAL8
      MPILIB=  -L/usr/pgi/linux86-64/5.2/lib -lfmpich -lmpich
      #MPI = -I/usr/lib/mpi/mpi_pgi/include -DUSE_MPI -DMPI_HAS_REAL8 -Dmpi_wti\
me=pmpi_wtime
      #MPILIB = -L/usr/lib/mpi/mpi_pgi/lib  -lfmpi -lmpi
   endif


   ifeq ($(VENDOR),LAHEY)
   #F90 = /usr/local/lf9562/bin/lf95
   F90 = lf95
   FFLAGS = -g  $(MPI) $(FFT)
   #FFLAGS =  -tpp -O $(MPI) $(FFT)
   LDFLAGS = --staticlink  $(FFTLIB)
   CC = gcc
   endif

   ifeq ($(VENDOR),PGI)
   F90 = pgf90
   CC=pgcc
   #FFLAGS = -Bstatic -Mrecursive -fast  $(MPI)  $(FFT) 
   FFLAGS = -Mrecursive -g  $(MPI)  $(FFT) 
   LDFLAGS = $(FFTLIB)
   endif 

   ifeq ($(VENDOR),PGIMPI)
   # this usually means Cray XT3 system
   #F90 = mpif90
   F90 = ftn
   CC = mpicc
   FFLAGS =  -Mrecursive  -fast  $(MPI)  $(FFT) 
   LDFLAGS = $(FFTLIB)
   endif 



   ifeq ($(VENDOR),INTEL)
   F90 = ifc
   CC=gcc
   FFLAGS =  -quiet -auto -w90  -Vaxlib  -O1    $(MPI)  $(FFT) 
   #FFLAGS =  -w90  -Vaxlib -O1 -tpp7   $(MPI)  $(FFT) 
   LDFLAGS = $(FFTLIB)
   endif

   ifeq ($(VENDOR),INTELMPI)
   F90 = mpif90
   CC=gcc
   FFLAGS =  -auto -w90  -Vaxlib  -O1  $(MPI) $(FFT) 
   LDFLAGS = $(FFTLIB)
   endif

   ifeq ($(VENDOR),GNU)
   F90 = gfortran
   CC=gcc
   FFLAGS =  -DGFORTRAN -O2  -ffixed-line-length-132 $(MPI) -I/usr/include $(FFT)  
   LDFLAGS = $(FFTLIB)
   endif

   ifeq ($(VENDOR),GNUMPI)
   F90 = mpif90
   CC=gcc
   FFLAGS =  -DGFORTRAN -O2  -ffixed-line-length-132 $(MPI) -I/usr/include $(FFT)  
   LDFLAGS = $(FFTLIB)
   endif



endif


ifeq ($(UNAME),Darwin)
#   MPI =  -DUSE_MPI  -DMPI_REAL8_LAM    -I/sw/include/lammpi
    MPI = -DUSE_MPI -DMPI_HAS_REAL8 -DUSE_MPI_IO
#   MPILIB = -L/sw/lib/  -llammpio -llamf77mpi -lmpi -llam 
   F90 = mpif90
   FFLAGS =  -O2  -DGFORTRAN -ffixed-line-length-132 $(MPI) 
#   CC = cc -DG77_UNDERSCORE

#   CC = cc -DF_NO_UNDERSCORE
#   F90 = xlf95
#   F77FIXED = -qfixed
#   FFLAGS = -WF,"-Dflush=flush_,-DNDEBUG"  -qsuffix=cpp=F90   -O2 -qmaxmem=-1 -qstrict 
#   LDFLAGS = -lSystemStubs
endif


ifeq ($(UNAME),SunOS)
   F90 = f95
   #MPILIB = -L/usr/local/src/mpich/lib/solaris/ch_p4  -lmpi -lsocket -lnsl
   #MPI = -DUSE_MPI -I/usr/local/src/mpich/include
   # default -fast is -xvector=yes, but this introduces bugs in test.F90!
   FFLAGS = -fast -xvector=no -xarch=v9   $(MPI) -DUSE_STKFFT
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
   CC = cc -64 $(LSFFLAGS)
endif


ifeq ($(UNAME),AIX)
   F90 = mpxlf95
   F77FIXED = -qfixed
   FFLAGS = -WF,"-Dflush=flush_,-DUSE_MPI,-DNDEBUG,-DMPI_HAS_REAL8"  -qsuffix=cpp=F90  -g -O3 -qmaxmem=-1 -qstrict 
   CC = cc -DF_NO_UNDERSCORE

   USE_ALPHA = -WF,"-DALPHA_MODEL"
   USE_LERAY = -WF,"-DALPHA_MODEL -DLERAY_MODEL"
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

   # R12
   MPI =  -DUSE_MPI -DUSE_MPI_IO $(MPI_COMPILE_FLAGS) -DMPI_HAS_REAL8
   MPILIB =  $(MPI_LD_FLAGS) -lfmpi -lmpio -lmpi  -lelan

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

#diagnostics for 3d aspect ratio runs:
ASPECT_DIAG = aspect_diag.o isoave.o subcubes.o

#diagnostics for 3d aspect ratio runs:
SLOW_DIAG = pdf.o slow_diag.o isoave.o subcubes.o

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
MODEL_SPEC_SLOW  = ns_d) $(FILEIO)
ecomp3.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC       = ns.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC_PENCIL = ns_xpencil.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC_PENCIL2 = ns_xpencil2.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
MODEL_SPEC_SCALE = ns.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(ASPECT_DIAG)
MODEL_SPEC_SL    = ns.o dns.o  $(OBJS1) $(OBJS2) $(OBJS3) $(SLOW_DIAG)
MODEL_SLAV       = ns_slaving.o dns.o $(OBJS1) $(OBJS2) $(OBJS3) $(TURB_DIAG)
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
dnsb: $(MODEL_SPEC_SCALE) 
	$(LINK) -o dnsb $(MODEL_SPEC_SCALE) $(LDFLAGS) $(MPILIB)
dnsslav:  $(MODEL_SLAV) 
	    $(LINK) -o dnsslav $(MODEL_SLAV) $(LDFLAGS) $(MPILIB)
dnssl: $(MODEL_SPEC_SL) 
	$(LINK) -o dnssl $(MODEL_SPEC_SL) $(LDFLAGS) $(MPILIB)
dnsp: $(MODEL_SPEC_PENCIL) 
	$(LINK) -o dnsp $(MODEL_SPEC_PENCIL) $(LDFLAGS) $(MPILIB)
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

analysis_slow:  analysis_slow.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_slow analysis_slow.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis2:  analysis2.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis2 analysis2.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_rms:  analysis_rms.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_rms analysis_rms.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_isoave:  analysis_isoave.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_isoave analysis_isoave.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_sk: analysis_sk.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_sk analysis_sk.o  $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)

analysis_spec2d: analysis_spec2d.o $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3)
	$(LINK) -o analysis_spec2d analysis_spec2d.o  $(TURB_DIAG) $(OBJS1) $(OBJS2) $(OBJS3) $(LDFLAGS) $(MPILIB)


clean: 
	rm -f temp* *.d *.o *.mod dns dnsb dnsslav dnsp dnss dnsa \
                            dnsgrid dnsghost dnsi dnsvor \
                            cdns dnsibm conv analysis_rms     \
                            analysis_isoave analysis analysis_sk \
                            analysis_spec2d splitspec \
			    convert	

dep:
	./mkdep.py *.F90 *.f  > .depends

-include .depends


ns_alpha.o: ns.F90 ns.o
	$(COMPILE) $(USE_ALPHA) -o ns_alpha.o ns.F90
ns_leray.o: ns.F90 ns.o
	$(COMPILE) $(USE_LERAY) -o ns_leray.o ns.F90

params.h:  
	echo "Run gridsetup.py to generage a params.h file!"
	exit 1


# fft_inteface, when using FFTW-GURU, depends on transpose.o, 
# but make dep system cannot detect this
fft_interface.o: macros.h params.o transpose.o makefile 

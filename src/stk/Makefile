####
#### $Id: Makefile,v 1.5 2001-08-20 18:33:27 mt Exp $
####

default:
	@echo "Usage: $(MAKE) sun | sgi | hp | ibm | alpha | linux | linuxpgc | macosx | generic | clean"

SRC_COMPLEX=cfft.c cfftAnalysis.c cfftSynthesis.c \
cfftAnalysisR2.c cfftAnalysisR3.c cfftAnalysisR4.c cfftAnalysisR5.c \
cfftAnalysisR8.c cfftAnalysisRg.c cfftSynthesisR2.c cfftSynthesisR3.c \
cfftSynthesisR4.c cfftSynthesisR5.c cfftSynthesisR8.c cfftSynthesisRg.c 

SRC_REAL=rfft.c rfftAnalysis.c rfftSynthesis.c \
rfftAnalysisR2.c rfftAnalysisR3.c rfftAnalysisR4.c rfftAnalysisR5.c \
rfftAnalysisR8.c rfftAnalysisRg.c rfftSynthesisR2.c rfftSynthesisR3.c \
rfftSynthesisR4.c rfftSynthesisR5.c rfftSynthesisR8.c rfftSynthesisRg.c 

SRC_FORTRAN=fortran.c

INC=cfft.h rfft.h real.h constants.h

OBJ_COMPLEX=$(SRC_COMPLEX:.c=.o)
OBJ_REAL=$(SRC_REAL:.c=.o)
OBJ_FORTRAN=$(SRC_FORTRAN:.c=.o)

$(OBJ_COMPLEX): cfft.h real.h constants.h
$(OBJ_REAL): rfft.h real.h constants.h
$(OBJ_FORTRAN): rfft.h cfft.h real.h

OBJ=$(OBJ_COMPLEX) $(OBJ_REAL) $(OBJ_FORTRAN)

TARGET=libstk.a t_fft test

sun:
#	@$(MAKE) CC=cc CFLAGS="-mt -xO5 -fsimple=2 -native -dalign -xarch=v9 -xrestrict" AR=/usr/ccs/bin/ar RANLIB=true LIBS="-L. -lstk -lm" $(TARGET)
	@$(MAKE) CC=cc CFLAGS="-xO5 -fsimple=2 -native -dalign -xarch=v9 -xrestrict" AR=/usr/ccs/bin/ar RANLIB=true LIBS="-L. -lstk -lm" $(TARGET)

sgi:
	@$(MAKE) CC=cc CFLAGS="-O2 -64 -OPT:alias=RESTRICT" AR=ar RANLIB=true LIBS="-L. -lstk -lm" $(TARGET)

hp:
	@$(MAKE) CC=c89 CFLAGS="-D_HPUX_SOURCE +e +DA2.0W +O3" AR=ar RANLIB=true LIBS="-L. -lstk -lm -lpthread" $(TARGET)

ibm:
	@$(MAKE) CC=xlc_r CFLAGS="-O3 -qstrict -qarch=pwr3 -qalias=allp -qlanglvl=ansi" AR="ar" RANLIB=true LIBS="-L. -lstk -lm" $(TARGET)

alpha:
	@$(MAKE) CC=cc CFLAGS="-O5 -arch ev6 -ansi_alias -std1" AR=ar RANLIB=true LIBS="-L. -lstk -lm -lpthread" $(TARGET)

linux:
	@$(MAKE) CC=cc CFLAGS="-O3" AR=ar RANLIB=true LIBS="-L. -lstk -lm -lpthread" $(TARGET)

linuxpgc:
	@$(MAKE) CC=pgcc CFLAGS="-fast -Msafeptr" AR=ar RANLIB=true LIBS="-L. -lstk -lm -lpthread" $(TARGET)

macosx:
	@$(MAKE) CC=cc CFLAGS="-O4" AR=ar RANLIB=ranlib LIBS="-L. -lstk -lm" $(TARGET)

generic:
	@$(MAKE) CC=cc CFLAGS=-O AR=ar RANLIB=ranlib LIBS="-L. -lstk -lm" $(TARGET) 

.c.o:
	$(CC) $(CFLAGS) -c  $<

$(OBJ): $(INC)

libstk.a: $(OBJ)
	$(AR) -rv libstk.a $(OBJ)
	$(RANLIB) libstk.a

t_fft: libstk.a t_fft.c
	$(CC) $(CFLAGS) t_fft.c -o $@ $(LIBS)

test: t_fft
	@echo "running test procedure..."
	@./t_fft

clean:
	rm -f  *.o libstk.a t_fft

tarfile:
	tar -cvf stk.tar $(SRC_COMPLEX) $(SRC_REAL) $(SRC_FORTRAN) $(INC) t_fft.c Makefile README LICENSE


# Makefile for dns
#
#

UNAME = $(shell uname)
ARCH = $(shell uname -m)

ifeq ($(UNAME),Linux)
   F90 = lf95
   LDFLAGS = --staticlink
endif

COMPILE=$(F90) -c $(FFLAGS)
LINK=$(F90) $(FFLAGS) $(LDFLAGS)

OBJS = util.o poisson.o

#independent modules:
BASICMODS = params.mod 

#modules which depend on the above, but not on each other:
ADVMODS = fft_interface.mod

#rules
%.mod: %.F90
	$(COMPILE) $< 
%.o: %.F90
	$(COMPILE) -c $< 
%.o: %.f
	$(COMPILE) -c -o $@ $< 


dnsargs = $(BASICMODS) $(ADVMODS) 
dns: $(dnsargs)
	$(LINK) $(LDFLAGS) $(dnsargs)


clean: 



fft_interface.mod:  fft99_interface.F90 params.mod
	$(COMPILE) fft99_interface.F90
#fft_interface.mod:  compaq_fft_interface.F90 params.mod
#	$(COMPILE) compaq_fft_interface.F90




 
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

OBJS = util.o poisson.o init_model.o ns3d.o

#independent modules (dont depend on any other modules):
MODS1 = params.mod
#modules which depend on MODS1, but not on each other
MODS2 = fft_interface.mod


#rules
%.mod: %.F90
	$(COMPILE) $< 
%.o: %.F90
	$(COMPILE) -c $< 
%.o: %.f
	$(COMPILE) -c -o $@ $< 


dnslink = $(OBJS)
dns: $(OBJS) $(MODS1) $(MODS2) 
	$(LINK) $(OBJS)


clean: 

dns.o:         dns.F90 params.mod
derivatives.o: derivataves.F90 params.mod
init_model.o:  init_model.F90 params.mod fft_interface.mod
ns3d.o:        ns3d.F90 params.mod
poisson.o:     poisson.F90 params.mod fft_interface.mod




fft_interface.mod:  fft99_interface.F90 params.mod
	$(COMPILE) fft99_interface.F90
#fft_interface.mod:  compaq_fft_interface.F90 params.mod
#	$(COMPILE) compaq_fft_interface.F90

ns3d.o:  ns3d.F90 params.mod



 
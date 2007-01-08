#!/usr/bin/env python
from math import *
#E=1.9; eps=3.8; kmax=1930*2*pi
#E=.5; eps=.0668; kmax=5792.
#E=1.79; eps=.506; kmax=482.
#E=1.97; eps=.591; kmax=482.
E=3.53; eps=1.302; kmax=482

Cr=1.0
C=.5

l=E**(3./2)/eps
param=l*kmax/Cr
Rl=2.582*param**(2./3.)
ntime=param/(.5*C)
nu = eps**(1/3.) * (Cr/kmax)**(4./3.)
print "nu=",nu
print "R_lambda=",Rl
print "t_e/dt = ",ntime


#!/usr/bin/env python
from string import *
import os, commands, getopt, sys

def nfactor(n,factor):
   nfacs=0
   nleft=n
   if (0==(n % factor)):
      [nleft,nfacs]=nfactor(n/factor,factor)
      nfacs=nfacs+1
   return [nleft,nfacs]

def fullfactor(n):
   facs=range(3)
   [nleft,facs[0]]=nfactor(n,2)
   [nleft,facs[1]]=nfactor(nleft,3)
   [nleft,facs[2]]=nfactor(nleft,5)

   if (nleft<>1):
      print "dimension: ",n[i],"must only have factors of 2,3 and 5"
      sys.exit(1)

   if n<>1 and (facs[0]+facs[1]+facs[2])<2:
      print "dimension: ",n,"must have at least 2 factors"
      sys.exit(1)
	
   return facs









if (len(sys.argv) <> 7):
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz'
   sys.exit(1) 



n=range(3)
ncpu=range(3)
nslab=range(3)
lmn=range(3)
for i in range(3):
   ncpu[i] = atoi(sys.argv[1+i])
for i in range(3):
   n[i] = atoi(sys.argv[4+i])


for i in range(3):
   if (0 <> n[i] % ncpu[i]):
      print "n[",i,"]",n[i],"not divisable by ncpu[",i,"]=",ncpu[i]
      sys.exit(1)
   nslab[i]=n[i]/ncpu[i]		

for i in range(3):
   j=(i+2) % 3
   if (0 <> nslab[j] % ncpu[i]):
      print "ncpu[",i,"]=",ncpu[i],"does not divide nslab[",j,"]:",nslab[j]
      sys.exit(1)
   lmn[i]=ncpu[i]/nslab[j]



# factor n into 2's, 3's and 5's.
facs=[range(3),range(3),range(3)]
for i in range(3):
   facs[i]=fullfactor(n[i])


# find a decompostion which looks like this:
# nx= ncpu_x*nslabx
# ny= ncpu_y*nslaby
# nz= ncpu_y*nslabz
#
# ncpu_x*l=nslabz
# ncpu_y*m=nslabx
# ncpu_z*n=nslaby
#


print "! number of prognostic variables:" 
print "integer,parameter :: n_var=3"

print "! parallel decomposition: "
print "integer,parameter :: ncpu_x=",ncpu[0]
print "integer,parameter :: ncpu_y=",ncpu[1]
print "integer,parameter :: ncpu_z=",ncpu[2]

print "! dimensions of grid and data:"
print "integer,parameter :: nx=",nslab[0]+2
print "integer,parameter :: ny=",nslab[1]+2
print "integer,parameter :: nz=",nslab[2]+2
print "integer,parameter :: nx1=1,nx2=",nslab[0]
print "integer,parameter :: ny1=1,ny2=",nslab[1]
print "integer,parameter :: nz1=1,nz2=",nslab[2]




    





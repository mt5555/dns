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

   if n<>1 and n<=4:
      print "dimension: ",n,"must be > 4"
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

use_x_z=1
use_x_y=1
for i in range(3):
   j=(i+2) % 3

   if (0 <> nslab[j] % ncpu[i]):
      print "ncpu[",i,"]=",ncpu[i],"does not divide nslab[",j,"]:",nslab[j]
      if i==0:
         print "cant use TRANSPOSE_X_SPLIT_Z"
         use_x_z=0
      else:
         sys.exit(1)

   # also try j=1
   if i==0:  
       j=1
       if (0 <> nslab[j] % ncpu[i]):
          print "ncpu[",i,"]=",ncpu[i],"does not divide nslab[",j,"]:",nslab[j]
          print "cant use TRANSPOSE_X_SPLIT_Y"
          use_x_y=0

       if use_x_y:   # use this value when possible
          use_x_z=0  
          j=1
          print 'Using TRANSPOSE_X_SPLIT_Y.  Updating transpose.h'
          

       if use_x_z:
          j=(i+2) % 3  # go back to original value
          print 'Using TRANSPOSE_X_SPLIT_Z.  Updating transpose.h'
          update_define(use_x_z)

   lmn[i]=ncpu[i]/nslab[j]



# factor n into 2's, 3's and 5's.
facs=[range(3),range(3),range(3)]
for i in range(3):
   facs[i]=fullfactor(n[i])

fout=open("transpose.h",'w')
fout.write("! Specify the 2D x-coordinate parallel decomposition\n")
fout.write("! (y and z 2D decompositions are fixed\n")
fout.write("\n")
fout.write("! default is TRANSPOSE_X_SPLIT_Y, which can be used with nslabz>=1 \n")
fout.write("! TRANSPOSE_X_SPLIT_Z gives the original parallel decomposition\n")
fout.write("! which cant be used with nslabz=1\n")
fout.write("\n")


if (use_x_z):
   fout.write("#define TRANSPOSE_X_SPLIT_Z\n")
   fout.write("#undef TRANSPOSE_X_SPLIT_Y\n")
else:
   fout.write("#undef TRANSPOSE_X_SPLIT_Z\n")
   fout.write("#define TRANSPOSE_X_SPLIT_Y\n")
fout.close()

print "Grid per cpu: ",nslab[0],nslab[1],nslab[2]," Updating user_params.h"


fout=open("user_params.h",'w')
fout.write("! number of prognostic variables:\n" )
fout.write("integer,parameter :: n_var=3\n")

fout.write("! parallel decomposition:\n ")
fout.write("integer,parameter :: ncpu_x="+str(ncpu[0])+"\n" )
fout.write("integer,parameter :: ncpu_y="+str(ncpu[1])+"\n")
fout.write("integer,parameter :: ncpu_z="+str(ncpu[2])+"\n")

fout.write("! dimensions of grid and data:\n")
fout.write("integer,parameter :: nx="+str(nslab[0]+2)+"\n")
fout.write("integer,parameter :: ny="+str(nslab[1]+2)+"\n")
zpad=2
if (nslab[2]==1):
   zpad=0
fout.write("integer,parameter :: nz="+str(nslab[2]+zpad)+"\n")
fout.write("integer,parameter :: nx1=1,nx2="+str(nslab[0])+"\n")
fout.write("integer,parameter :: ny1=1,ny2="+str(nslab[1])+"\n")
fout.write("integer,parameter :: nz1=1,nz2="+str(nslab[2])+"\n")
fout.close()


    





#!/usr/bin/env python
from string import *
import os, commands, getopt, sys
####################################################################################
#
#  gridsetup python script.  main program.
#
#
####################################################################################
padset=0
pad1=[0,0,0]
pad2=[2,2,2]

if (len(sys.argv)== 13):
   padset=1
   for i in range(3):
      pad1[i]=atoi(sys.argv[10+i])
      pad2[i]=atoi(sys.argv[7+i])
elif len(sys.argv)==10:
   padset=1
   for i in range(3):
      pad2[i]=atoi(sys.argv[7+i])
elif len(sys.argv)== 7:
   i=0
   #nothing
else:
   #arg       0      1      2      3     4  5  6   7     8   9      10      11   12         
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz'
   print 'OR:'
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz  padx pady padz'
   print 'OR:'
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz  padx pady padz offsetx offsety offsetz'
   sys.exit(1) 






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





n=range(3)
ncpu=range(3)
nslab=range(3)
lmn=range(3)
for i in range(3):
   ncpu[i] = atoi(sys.argv[1+i])
for i in range(3):
   n[i] = atoi(sys.argv[4+i])

# override default of 2 in z direction for 2D problems
if not padset:
   if n[2]==1 :
      pad2[2]=0
      


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
          print 'Using TRANSPOSE_X_SPLIT_Y. '

       if use_x_z:
          use_x_y=0
          j=(i+2) % 3  # go back to original value
          print 'Using TRANSPOSE_X_SPLIT_Z. '


          


   lmn[i]=ncpu[i]/nslab[j]



# factor n into 2's, 3's and 5's.
facs=[range(3),range(3),range(3)]
for i in range(3):
   facs[i]=fullfactor(n[i])


print sys.argv[1:-1]
print "Grid per cpu: ",nslab[0],nslab[1],nslab[2]


fout=open("/tmp/transpose.h",'w')
fout.write("! Specify the 2D x-coordinate parallel decomposition\n")
fout.write("! (y and z 2D decompositions are fixed\n")
fout.write("\n")
fout.write("! default is TRANSPOSE_X_SPLIT_Y, which can be used with nslabz>=1 \n")
fout.write("! TRANSPOSE_X_SPLIT_Z gives the original parallel decomposition\n")
fout.write("! which cant be used with nslabz=1\n")
fout.write("\n")

cmd = "diff /tmp/transpose.h transpose.h"
(status,out) = commands.getstatusoutput(cmd) 

if (status != 0):
   print "Creating new transpose.h file"
   cmd = "mv -f /tmp/transpose.h transpose.h"
   status = os.system(cmd)
else:
   print "transpose.h file unchanged"




if (use_x_z):
   fout.write("#define TRANSPOSE_X_SPLIT_Z\n")
   fout.write("#undef TRANSPOSE_X_SPLIT_Y\n")
else:
   fout.write("#undef TRANSPOSE_X_SPLIT_Z\n")
   fout.write("#define TRANSPOSE_X_SPLIT_Y\n")
fout.close()


fout=open("/tmp/params.h",'w')
fout.write("! number of prognostic variables:\n" )
fout.write("integer,parameter :: n_var=3\n")

fout.write("! parallel decomposition:\n ")
fout.write("integer,parameter :: ncpu_x="+str(ncpu[0])+"\n" )
fout.write("integer,parameter :: ncpu_y="+str(ncpu[1])+"\n")
fout.write("integer,parameter :: ncpu_z="+str(ncpu[2])+"\n")

tot=[0,0,0]
for i in range(3):
   tot[i]=nslab[i]+pad1[i]+pad2[i]
   pad2[i]=pad1[i]+nslab[i]
   pad1[i]=1+pad1[i]


fout.write("! dimensions of grid and data:\n")
fout.write("integer,parameter :: nx="+str(tot[0])+"\n")
fout.write("integer,parameter :: ny="+str(tot[1])+"\n")
fout.write("integer,parameter :: nz="+str(tot[2])+"\n")
fout.write("integer,parameter :: nx1="+str(pad1[0])+",nx2="+str(pad2[0])+"\n")
fout.write("integer,parameter :: ny1="+str(pad1[1])+",ny2="+str(pad2[1])+"\n")
fout.write("integer,parameter :: nz1="+str(pad1[2])+",nz2="+str(pad2[2])+"\n")
fout.close()


cmd = "diff /tmp/params.h params.h"
#status = os.system(cmd)  #display output on screen
(status,out) = commands.getstatusoutput(cmd) 


if (status != 0):
   print "Creating new params.h file"
   cmd = "mv -f /tmp/params.h params.h"
   status = os.system(cmd)
else:
   print "params.h file unchanged"


    





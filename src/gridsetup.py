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
n_var=3


if (len(sys.argv)>= 14):
   n_var=atoi(sys.argv[13])

if (len(sys.argv)>= 13):
   for i in range(3):
      pad1[i]=atoi(sys.argv[10+i])

if len(sys.argv)>=10:
   padset=1
   for i in range(3):
      pad2[i]=atoi(sys.argv[7+i])

if len(sys.argv)>= 7:
   i=0
   #nothing
   
if len(sys.argv)<7:
   #arg       0      1      2      3     4  5  6   7     8   9      10      11   12         
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz'
   print 'OR:'
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz  padx pady padz'
   print 'OR:'
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz  padx pady padz offsetx offsety offsetz'
   print 'OR:'
   print 'gridsetup ncpu_x ncpu_y ncpu_z nx ny nz  padx pady padz offsetx offsety offsetz  n_var'
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
      print "WARNING: dimension: ",n,"must only have factors of 2,3 and 5"
#      sys.exit(1)

#   if n<>1 and n<=4:
#      print "dimension: ",n,"must be > 4"
#      sys.exit(1)
	
   return facs





n=range(3)
ncpu=range(3)
nslab=range(3)
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
      print "ERROR: n[",i,"]",n[i],"not divisable by ncpu[",i,"]=",ncpu[i]
      sys.exit(1)
   nslab[i]=n[i]/ncpu[i]		

for i in range(3):
   if (i==0):
      j=1
   if (i==1):
      j=0
   if (i==2):
      j=1

   if (0 <> nslab[j] % ncpu[i]):
      print "WARNING: ncpu[",i,"]=",ncpu[i],"does not divide nslab[",j,"]:",nslab[j]




# factor n into 2's, 3's and 5's.
facs=[range(3),range(3),range(3)]
for i in range(3):
   facs[i]=fullfactor(n[i])


print sys.argv[1:]
print "Total cpus:   ",ncpu[0]*ncpu[1]*ncpu[2],
print "Grid per cpu: ",nslab[0],nslab[1],nslab[2]




fout=open("/tmp/params.h",'w')
fout.write("! number of prognostic variables:\n" )
fout.write("integer,parameter :: n_var="+str(n_var)+"\n")

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
fout.write("integer,private,parameter :: nxd="+str(tot[0])+"\n")
fout.write("integer,private,parameter :: nyd="+str(tot[1])+"\n")
fout.write("integer,private,parameter :: nzd="+str(tot[2])+"\n")
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
   cmd = "rm -f /tmp/params.h "
   status = os.system(cmd)
   print "params.h file unchanged"


    





#!/usr/bin/env python
from string import *
import os, commands, getopt, sys
####################################################################################
#
#  ./lsfcheck.py   filename N
#
# reads a number M from the first line of filename,
# replaces it with M+1
#
# output:     FILE_ERROR        error accessing file
#             SMALLER           M < N
#             LARGER            M > N
#             EQUAL             M==N
# 
#
####################################################################################

if (len(sys.argv)!= 3):
   print 'usage:   ./lsfcheck.py filename N'
   print 'output:  FILE_ERROR, M_LESS_N, M_GREATER_N, M_EQUAL_N'
   sys.exit(1)
else:
   file=sys.argv[1]
   target=atoi(sys.argv[2])


try:
   fid=open(file,'r')
   line=fid.readline()
   
   fvalue=atoi(line)
   fvalue=fvalue+1
   fid.close()
   fid=open(file,'w')
   fid.write(str(fvalue))

   if (fvalue<target):
      print 'SMALLER'
   elif (fvalue>target):
      print 'LARGER'
   elif (fvalue==target):
      print 'EQUAL'
   else:
      print 'FILE_ERROR'
   

except IOError,e:
   print 'FILE_ERROR'
   sys.exit(1)





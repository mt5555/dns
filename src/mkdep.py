#!/usr/bin/env python
from string import *
import os, commands, getopt, sys



#######################################################################
#
# take a name.F90 and return name.o
#
#######################################################################
def getdoto(file):
    tmp=split(file,".")
    extension=tmp[-1]
    fileo=file[:-(len(extension)+1)] + ".o"
    return fileo


#######################################################################
#
# parse file, looking for use module  or #include "name.h"
#
#######################################################################
def parse(file):
    fid=open(file)
    line=fid.readline()
    deps=[]

    while line:
        line=rstrip(line)
        tmp=split(line)
        if (len(tmp)>=2):
            keyword=tmp[0]
            use= -1 <> find(keyword,"use")
            inc= -1 <> find(keyword,"include")
            if (use or inc):
                name=tmp[1]
                name=strip(name)
                if name[0]=='"' :
                    name=name[1:]
                if name[-1]=='"' :
                    name=name[:-1]
                # dont include system directories:
                if (name[0:4] <>"/usr") and (name[0:4] <> "/opt") and (name[0:4] <> "mpif") :
                    if use:
                        name = name + ".o"
                    if name not in deps:
                        deps.append(name)
                        
        line=fid.readline()
    return deps


if (len(sys.argv) <= 1):
   print 'Run this script via "make dep"'
   exit 

for file in sys.argv[1:]:
    fileo=getdoto(file)
    deps=parse(file)
    l=len(deps)
    print fileo + ":", 
    if (l>0):
        for i in range(0,len(deps)):
            print deps[i],
    print "makefile",
    print ""

        
    

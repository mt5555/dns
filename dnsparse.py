#!/usr/bin/env python
from string import *
import os, commands, getopt, sys, exceptions

class search_failed(exceptions.Exception):
    def __init__(self,args=None):
        self.args=args
class eof(exceptions.Exception):
    def __init__(self,args=None):
        self.args=args




#
# look for key1.  if we find key2 instead, raise error
#
def lookfor1(fid,key1,key2="",allow_eof=0):
    line=fid.readline()
    while line:
        pos=find(line,key1)
        if (-1 <> pos ):
            sline = line[pos+len(key1):-1]
            return sline
#            if (len(sline)==0):
#                return "shortline"
#            else: 
#               return sline[0]
        if (len(key2)>0):
           pos=find(line,key2)
           if (-1 <> pos ):
               print "error looking for: "+key1
               raise search_failed,"run not complete, found: "+key2
        line=fid.readline()

    if (allow_eof==1):
        raise eof,"EOF"
    raise search_failed,"Search failed  string='"+key1+"'"
    return 0

    


tsolve={}
try:
    while 1:
        startstr="Running parallel.  NCPUS="
        str = lookfor1(sys.stdin,startstr,"",1)
        str=split(str)
        np=atoi(str[0])
        nx=atoi(str[2])
        ny=atoi(str[4])
        nz=np/nx/ny  # sometimes nz is *** in the output file

        str=lookfor1(sys.stdin,"Global grid:",startstr)
        str=split(str)
        N=atoi(str[0])
        
        key=(N,np,nx,ny,nz)                        
        # now look for time
        str=lookfor1(sys.stdin,"dns_solve:",startstr)
        str=split(str)
        time=atof(str[5])  # time in min

        if tsolve.has_key(key):
            x=[]
            x.append(tsolve[key])
            x.append(time)
            tsolve[key]=x
        else:
            tsolve[key]=[]
            tsolve[key].append(time)

        
        
except search_failed,e:
    print "search failed"
    print e
    sys.exit(1)
    
except eof,e:
    print '%% time(m) per timestep'
    x=tsolve.keys()
    x.sort()
    res_last=0
    for k in x:
        if (k[0]!=res_last):
            print '%% res = %i' % k[0]
            res_last=k[0]
            count=0
        count=count+1
        best=min(tsolve[k])   # time for 1 timesetp, seconds
        print "ncpu%i(%i)=%i; ncpux%i(%i)=%i; ncpuz%i(%i)=%i;  time%i(%i)=%e  "   %          (k[0],count,k[1],k[0],count,k[2],k[0],count,k[4],k[0],count,best)


    sys.exit(0)


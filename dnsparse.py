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

        # create empty list if needed, then append new data
        tsolve.setdefault(key,[]).append(time)

        
        
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
        print "ncpu%i(%i)=%i;  time%i(%i)=%e; %% (%2i,%i,%4i)  res/nz=%i "   %          (k[0],count,k[1],k[0],count,best,k[2],k[3],k[4],res_last/k[4])

    tbest={}
    tbest_nx={}
    tbest_nz={}
    # key=(N,np,nx,ny,nz)                        
    for k in x:
        key=(k[0],k[1])
        if (tbest.has_key(key)):
            tbest[key] = tbest[key]+ tsolve[k]
            tbest_nx[key] = tbest_nx[key]+ [ k[2] ]
            tbest_nz[key] = tbest_nz[key]+ [ k[4] ]
        else:
            tbest[key]=tsolve[k]
            tbest_nx[key] = [ k[2] ]
            tbest_nz[key] = [ k[4] ]


    x=tbest.keys();
    x.sort();
    res_last=0
    for k in x:
        if (k[0]!=res_last):
            print '%% res = %i' % k[0]
            print 'nbest%i=[]; tbest%i=[]; eddypd%i=[];' % (k[0],k[0],k[0])
            res_last=k[0]
            count=0
        count=count+1
        best=min(tbest[k])

        i = tbest[k].index(best)
        print "%% res=%i NCPU=%i best was: (%i,%i,%i)  res/nx=%i res/nz=%i" % (res_last,k[1],tbest_nx[k][i],1,tbest_nz[k][i],res_last/tbest_nx[k][i],res_last/tbest_nz[k][i])

        eperd = 3333*res_last/512.  # number of timesteps
        eperd = eperd*best     # time in min for 1 eddy turnover
        eperd = 24*60./eperd   # eddy turnover per day
        print "nbest%i(%i)=%5i; tbest%i(%i)=%e; eddypd%i(%i)=%e;" % (k[0],count,k[1],k[0],count,best,k[0],count,eperd)


    sys.exit(0)


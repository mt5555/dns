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
tfft={}
tr_x={}
tr_y={}
tr_z={}
tr_tot={}
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
        time=atof(str[5])  # time in min per timestep
        time_tot=atof(str[1])  # total time in min
#        nstep=int(.1+time_tot/time)  # number of timesteps timed
        nstep=3
        
        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"transpose_to_z ",startstr)
        str=split(str)
        time_to_z=atof(str[0])/nstep/12
        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"transpose_from_z ",startstr)
        str=split(str)
        time_from_z=atof(str[0])/nstep/20
        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"transpose_to_x ",startstr)
        str=split(str)
        time_to_x=atof(str[0])/nstep/24
        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"transpose_from_x ",startstr)
        str=split(str)
        time_from_x=atof(str[0])/nstep/12
        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"transpose_to_y ",startstr)
        str=split(str)
        time_to_y=atof(str[0])/nstep/32
        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"transpose_from_y ",startstr)
        str=split(str)
        time_from_y=atof(str[0])/nstep/36

        # now look for tranpose time.  note typo in output
        str=lookfor1(sys.stdin,"traspose total ",startstr)
        str=split(str)
        time_tr=atof(str[0])/nstep



        # now look for FFT time.  note typo in output
        str=lookfor1(sys.stdin,"FFT ",startstr)
        str=split(str)
        time_fft=atof(str[0])
        str=lookfor1(sys.stdin,"iFFT ",startstr)
        str=split(str)
        time_fft += atof(str[0])
        time_fft = time_fft/nstep

        # create empty list if needed, then append new data
        tsolve.setdefault(key,[]).append(time)
        tfft.setdefault(key,[]).append(time_fft)
        tr_x.setdefault(key,[]).append((time_to_x + time_from_x) /2)
        tr_y.setdefault(key,[]).append((time_to_y + time_from_y) /2)
        tr_z.setdefault(key,[]).append((time_to_z + time_from_z) /2)
        tr_tot.setdefault(key,[]).append((time_tr))
        
        
        
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
    tbest_fft={}
    tbest_nx={}
    tbest_nz={}
    tbest_tr_x={}
    tbest_tr_y={}
    tbest_tr_z={}
    tbest_tr_tot={}
    # key=(N,np,nx,ny,nz)                        
    for k in x:
        key=(k[0],k[1])
        if (tbest.has_key(key)):
            tbest[key] = tbest[key]+ tsolve[k]
            tbest_nx[key] = tbest_nx[key]+ [ k[2] ]
            tbest_nz[key] = tbest_nz[key]+ [ k[4] ]
            tbest_tr_x[key] = tbest_tr_x[key] + tr_x[k]
            tbest_tr_y[key] = tbest_tr_y[key] + tr_y[k]
            tbest_tr_z[key] = tbest_tr_z[key] + tr_z[k]
            tbest_tr_tot[key] = tbest_tr_tot[key] + tr_tot[k]
            tbest_fft[key] = tbest_fft[key] + tfft[k]
        else:
            tbest[key]=tsolve[k]
            tbest_nx[key] = [ k[2] ]
            tbest_nz[key] = [ k[4] ]
            tbest_tr_x[key] = tr_x[k]
            tbest_tr_y[key] = tr_y[k]
            tbest_tr_z[key] = tr_z[k]
            tbest_tr_tot[key] = tr_tot[k]
            tbest_fft[key] = tfft[k]


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
        print "tr_x%i(%i)=%e; tr_y%i(%i)=%e; tr_z%i(%i)=%e;" %         (k[0],count,tbest_tr_x[k][i],k[0],count,tbest_tr_y[k][i],k[0],count,tbest_tr_z[k][i] )
        print "tr_tot%i(%i)=%e; " %         (k[0],count,tbest_tr_tot[k][i])
        print "fft%i(%i)=%e; " %            (k[0],count,tbest_fft[k][i])

        # bi-direciton bandwidth: each transpose is sending 

    sys.exit(0)


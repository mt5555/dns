#!/usr/bin/env python
from string import *
from cStringIO import StringIO
import os, commands, getopt, sys


filenameu = "temp.u"
filenamev = "temp.v"
filenamew = "temp.w"

fu = open(filenameu,mode='rb')
fv = open(filenamev,mode='rb')
fw = open(filenamew,mode='rb')
fi = open('interleave.uvw',mode='wb')



i=0
ntot = 0
bsize = 8*1024*1024   # read 64MB at a time.  need 384MB to run
while 1:
    i=i+1
    print 'ntot = %i' % ntot
    datau = fu.read(8*bsize)
    datav = fv.read(8*bsize)
    dataw = fw.read(8*bsize)
    nread = len(datau)
    if ( nread % 8 != 0):
        print 'Error reading file - did not contain multiple of 8 bytes'
        sys.exit(0)
    nread = nread/8
    ntot += nread
    if (nread == 0 ):
        print 'hit EOF. stopping'
        print 'processed N=%i real*8 numbers from each file' % ntot
        fi.close() 
        sys.exit(0)

    outbuf=StringIO()
    for j in range(nread):
        j1 = j*8
        j2 = j1 + 8
        outbuf.write(datau[j1:j2])
        outbuf.write(datav[j1:j2])
        outbuf.write(dataw[j1:j2])
    
    fi.write(outbuf.getvalue())
    






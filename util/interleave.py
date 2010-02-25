#!/usr/bin/env python
from string import *
import os, commands, getopt, sys


filenameu = "temp.u"
filenamev = "temp.v"
filenamew = "temp.w"

fu = open(filenameu,mode='rb')
fv = open(filenamev,mode='rb')
fw = open(filenamew,mode='rb')
fi = open('interleave.uvw',mode='wb')

i=0
while 1:
    i=i+1
    if (i % 100000 == 0 ):
        print 'processing N=%i' % i
    data = fu.read(8)
    if (data == '' ):
        print 'hit EOF. stopping'
        print 'processed N=%i numbers from each file' %i
        sys.exit(0)
    fi.write(data)
    data = fv.read(8)
    fi.write(data)
    data = fw.read(8)
    fi.write(data)
    






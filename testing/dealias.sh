#!/bin/tcsh -f
#
#
# Check aliasing errors:
#
# 1. on an NxNxN grid:
#    generate a U, apply truncation
#    output spectral coefficients
#    output UN = omega cross u (dealised)
#
# 2. on an MxMxM grid (M=2N)
#    read in U, UN
#    compute UN2 (on MxMxM grid)
#    compare UN2 with UN.  should agree in all modes in N-truncation
#
#
#


set NCPU = 1

if ( $NCPU == 1 ) then
   set MPIRUN = " "
else
   set MPIRUN = "mpirun -np $NCPU "
endif

#./gridsetup.py 1 1 $NCPU 32 32 32
#make -j2 convert ; mv convert convert.32

./gridsetup.py 1 1 $NCPU 96 96 96
make -j2 convert ; mv convert convert.64


#set method = "fft-dealias"
set method = "fft-sphere"
#set method = "fft-phase"

echo $method 
rm -f /tmp/temp0000.0000* /tmp/dealias*inp /tmp/temp1.out /tmp/temp2.out
sed s/METHOD/"$method"/  ../testing/dealias.inp > /tmp/dealias.inp


$MPIRUN ./convert.32 -so -cout nlout  -i /tmp/dealias.inp -d /tmp  | tee /tmp/temp1.out
$MPIRUN ./convert.64 -si -cout nlin  -i /tmp/dealias.inp -d /tmp   | tee /tmp/temp2.out

grep "number of retained modes:" /tmp/temp1.out
grep "number of modes" /tmp/temp2.out
grep "number of non-zero" /tmp/temp2.out
grep "max error over " /tmp/temp2.out








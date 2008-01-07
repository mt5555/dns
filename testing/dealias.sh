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
   ./gridsetup.py 1 1 $NCPU 32 32 32
   make dns
else
   set MPIRUN = "mpirun -np $NCPU "
   ./gridsetup.py 1 1 $NCPU 32 32 32
   make -j 2 dns
endif







#
# fft-dealias >=:  -.17e-5  -.56e-7  -.19e-8  -.66e-10 -.82e-11 -.15e-10
#             >    -.10e-4  -.34e-6  -.11e-7  -.41e-9  -.21e-10 -.28e-10
# fft-sphere    :  -.19e-2  -.58e-4  -.17e-5  -.47e-7  -.16e-8  -.94e-10

set i=0
loop:

if ( $i == 0 ) set method = "fft-dealias"
if ( $i == 1 ) set method = "fft-sphere"
if ( $i == 2 ) set method = "fft-phase"

set dt = 0.005
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
#$MPIRUN./dns -i /tmp/dealias.inp -d /tmp | tee  /tmp/temp.out
$MPIRUN./dns -i /tmp/dealias.inp -d /tmp >   /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out

set dt = .0028
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
$MPIRUN./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out

set dt = .00157
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
$MPIRUN./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out

set dt = .000883
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
$MPIRUN./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out


set dt = .000497
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
$MPIRUN./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out



echo
echo 



@ i += 1
if ( $i < 3 ) goto loop







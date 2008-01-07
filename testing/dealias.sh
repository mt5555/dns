#!/bin/tcsh -f
./gridsetup.py 1 1 2 32 32 32
make -j 2 dns

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
#mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp | tee  /tmp/temp.out
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >   /tmp/temp.out
grep "ke:" /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out

set dt = .0025
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out

set dt = .00125
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out

set dt = .000625
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out


set dt = .000001
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out



echo
echo 



@ i += 1
if ( $i < 3 ) goto loop







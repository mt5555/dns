#!/bin/tcsh -f
./gridsetup.py 1 1 2 32 32 32
make -j 2 dns



set i=0
loop:

if ( $i == 1 ) set method = "fft-dealias"
if ( $i == 0 ) set method = "fft-sphere"
if ( $i == 2 ) set method = "fft-phase"

set dt = .005
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
#mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp | tee  /tmp/temp.out
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >   /tmp/temp.out
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


set dt = .0003125
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out


set dt = .00015626
echo $method :  dt = $dt
rm -f /tmp/temp.out /tmp/dealias*inp
sed s/DT/$dt/ ../testing/dealias.inp | \
sed s/METHOD/"$method"/  > /tmp/dealias.inp
mpirun -np 2 ./dns -i /tmp/dealias.inp -d /tmp >  /tmp/temp.out
grep "entire run:" /tmp/temp.out
grep "per timestep, min/max" /tmp/temp.out



@ i += 1
if ( $i < 3 ) goto loop







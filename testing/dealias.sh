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


set NCPU1 = 1
set NCPU2 = 1
set NCPU = 0
@ NCPU = $NCPU1 * $NCPU2

if ( $NCPU == 1 ) then
   set MPIRUN = " "
else
   set MPIRUN = "mpirun -np $NCPU "
endif


# on M grid, modes 2N-M will alias
# contaiminated = M-(2N-M) = 2(M-N)
# we want 2(M-N) >= N      2M >= 3N  M>= 3N/2
set M = 72   # allows for N<=48
set N = 40

if ! (-x convert.big)  then
  ./gridsetup.py 1 1 1 $M $M $M
  make -j2 convert ; mv convert convert.big
endif

./gridsetup.py $NCPU1 1 $NCPU2 $N $N $N
make -j2 convert 



#set method = "fft-sphere"
#set method = "fft-23sphere"
#set method = "fft-dealias"
set method = "fft-phase"

echo $method 
rm -f /tmp/temp0000.0000* /tmp/dealias*inp /tmp/temp1.out /tmp/temp2.out
sed s/METHOD/"$method"/  ../testing/dealias.inp > /tmp/dealias.inp


$MPIRUN ./convert -smax 999 -so -cout nlout  -i /tmp/dealias.inp -d /tmp  | tee /tmp/temp1.out
./convert.big -si -cout nlin  -i /tmp/dealias.inp -d /tmp   | tee /tmp/temp2.out
echo 
echo =======================================================
echo "Initial condition code N=$N"
grep "number of retained modes:" /tmp/temp1.out
echo "Errors compared to N=$M calculation:"
grep "number of modes" /tmp/temp2.out
grep "number of non-zero" /tmp/temp2.out
grep "max error over " /tmp/temp2.out

exit

RESULTS: (for N=32, M=72. )
                   # modes      
                   tested     error
             
fft-sphere         14363      2.5e-3     aliasing error in waves 2-15
fft-sphere(new)    12893      2.0e-3     aliasing in waves 4-15

fft-23dealias       5041      1e-18
fft-23dealias(new)  4945      3.0e-17 

fft-dealais         9261      5.6e-17

fft-phase(new)     12893      6.7e-17


check for zero error:  ff
N=48  fft-phase, fft-23dealias, fft-dealias  error:  ~ 1e-17
N=40  fft-phase, fft-23dealias, fft-dealias  error:  ~ 1e-17
N=24  fft-phase, fft-23dealias, fft-dealias  error:  ~ 1e-16

                      

notes: 
fft-sphere: 
  -smax 8   no aliasing error  8*2 = 16, resolved on our grid
  -smax 9   aliasing error     9*2 = 18, aliasing to 14,15
  -smax 10                    10*2 = 20, aliasing to 12-15
  -smax 999                   15*2 = 30, aliasing to 2-15






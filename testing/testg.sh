#/bin/csh -f

cd ../src
set refin=../testing/ghost2d.inp
set refout=../testing/ghost2d.out
set tmp=/tmp/temp.out

if ($#argv == 0 ) then
   echo "./test.sh [1,p,makeref]"
   echo " 1 = run 1 simple 2D test case"
   echo " p = run some parallel cases"
   echo " makeref  = generate new reference output, 2D"
   exit
endif

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 128 128 1 
   make dnsgrid; rm -f $refout 
   ./dnsgrid -d /tmp < $refin > $refout
  cat $refout

endif

if ($1 == 1) then

./gridsetup.py 1 1 1 128 128 1  2 2 0 2 2 0
make dnsgrid >& /dev/null ;  rm -f $tmp ; ./dnsgrid -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 1 1 128 128 1  2 2 0 2 2 0
make dnsghost >& /dev/null ;  rm -f $tmp ; ./dnsghost -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout


endif

if ($1 == p) then


./gridsetup.py 2 1 1 128 128 1  2 2 0 2 2 0
make dnsghost >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dnsghost -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 2 1 128 128 1  2 2 0 2 2 0
make dnsghost >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dnsghost -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 2 2 1 128 128 1  2 2 0 2 2 0
make dnsghost >& /dev/null ;  rm -f $tmp ; mpirun -np 4 ./dnsghost -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout


endif


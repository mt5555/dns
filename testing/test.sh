#/bin/csh -f

cd ../src
set refin=../testing/reference2d.inp
set refout=../testing/reference2d.out
set tmp=/tmp/temp.out

if ($#argv == 0 ) then
   echo "./test.sh [1,2,3,3p]"
   echo " 1 = run 1 simple 2D test case"
   echo " 2 = run lots of 2D test cases (different dimensions)"
   echo " 3 = run several 3D test cases (different dimensions)"
   echo " 3p = run several 3D test cases in parallel"
   exit
endif

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 128 128 1 
   make ; rm -f $refout 
   ./dns < $refin > $refout
  cat $refout

endif

if ($1 == 3p) then

  ./gridsetup.py 1 1 2 128 128 8  
  make >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout

  ./gridsetup.py 1 1 4 128 128 8  3 2 1  5 4 3
  make >& /dev/null ;  rm -f $tmp ; mpirun -np 4 ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout



endif

if ($1 == 3) then

  ./gridsetup.py 1 1 1 128 128 6  0 0 0  2 2 2
  make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout

  ./gridsetup.py 1 1 1 128 128 6  3 2 1  5 4 3
  make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout

  ./gridsetup.py 1 1 1 128 128 6  2 2 2  0 0 0 
  make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout

endif

if ($1 == 1) then

./gridsetup.py 1 1 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout

endif

if ($1 == 2) then

./gridsetup.py 1 1 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout



./gridsetup.py 1 1 1 128 128 1  0 0 0  2 2 0
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout


./gridsetup.py 1 1 1 128 128 1  2 2 2  0 0 0 
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout



./gridsetup.py 1 1 1 128 128 1  2 2 1  0 0 0 
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout


./gridsetup.py 1 1 1 128 128 1  2 3 4  0 0 0 
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout


./gridsetup.py 1 1 1 128 128 1  4 2 3  0 0 0 
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout



./gridsetup.py 1 1 1 128 128 1  0 0 0  2 3 1
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout


./gridsetup.py 1 1 1 128 128 1  0 0 0  3 2 1
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout



./gridsetup.py 1 1 1 128 128 1  0 0 0  1 2 3
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout



./gridsetup.py 1 1 1 128 128 1  2 3 4  1 2 3
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 1 1 128 128 1  4 2 3  2 3 1
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout




endif


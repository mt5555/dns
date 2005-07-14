#/bin/csh -f

cd ../src
set refin=../testing/reference2d.inp
set refout=../testing/reference2d.out
set tmp=/tmp/temp.out
set rundir=../testing/2d

if ($#argv == 0 ) then
   echo "./test.sh [1,2,3,3p]"
   echo " "
   echo "Run a 2D test problem, with a 2 or 3 dimensional grid "
   echo " "
   echo " 1 = run 1 simple 2D test case"
   echo " 1r = run 1 simple 2D restart test case"
   echo " 1p = run 2D restart test case in parallel (not yet coded)"
   echo " "
   echo " 2 = run lots of 2D test cases (different dimensions)"
   echo " 3 = run several 3D test cases (different dimensions)"
   echo " 3p = run several 3D test cases in parallel (2 and 4 cpus)"
   echo " makeref  = generate new reference output, 2D"
   exit
endif

#set MPIRUN = "mpirun -x LD_LIBRARY_PATH -wd $PWD "
set MPIRUN = "mpirun"

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 128 128 1 
   make ; rm -f $refout 
   ./dns < $refin > $refout
  cat $refout
  mv reference2d0000.0000.u $rundir/restart.u
  mv reference2d0000.0000.v $rundir/restart.v

endif

if ($1 == 3p) then


  ./gridsetup.py 1 1 2 128 128 8  
  echo "compiling..."
  make >& /dev/null   
  echo "done compiling. status=" $status
  rm -f $tmp ; $MPIRUN -np 2 ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout

  ./gridsetup.py 1 1 4 128 128 8  3 2 1  5 4 3
  make >& /dev/null ;  rm -f $tmp ; $MPIRUN -np 4 ./dns < $refin > $tmp 
  ../testing/check.sh $tmp $refout

  ./gridsetup.py 2 2 1 128 128 8  
  make >& /dev/null ;  
  rm -f $tmp ; $MPIRUN -np 4 ./dns < $refin > $tmp 
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

if ($1 == 1r) then

./gridsetup.py 1 1 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; 
echo "running restart dns case"
./dns -r -d $rundir < $refin > $tmp 
../testing/check.sh $tmp $refout

endif

if ($1 == 1p) then

./gridsetup.py 1 1 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; 
echo "running restart dns 1x1x1"
$MPIRUN -np 1 ./dns -r -d $rundir -i  $refin > $tmp 
../testing/check.sh $tmp $refout


./gridsetup.py 8 1 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; 
echo "running restart dns 8x1x1"
$MPIRUN -np 8 ./dns -r -d $rundir -i  $refin > $tmp 
../testing/check.sh $tmp $refout

exit

./gridsetup.py 1 8 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; 
echo "running restart dns 1x8x1"
$MPIRUN -np 8 ./dns -r -d $rundir -i  $refin > $tmp 
../testing/check.sh $tmp $refout

endif


if ($1 == 1) then

./gridsetup.py 1 1 1 128 128 1  
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 1 1 128 128 1  
make dnsgrid >& /dev/null ;  rm -f $tmp ; ./dnsgrid < $refin > $tmp 
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


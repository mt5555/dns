#/bin/csh -f

cd ../src
set refrestart=../testing/restart3d.inp
set refin=../testing/reference3d.inp
set refout=../testing/reference3d.out
set tmp=/tmp/temp.out

if ($#argv == 0 ) then
   echo "./test3d.sh [1,2,3,3p]"
   echo " 1 = run 2 (with and without restart) simple 3D test case"
   echo " g = run 1 simple 3D test case using dnsgrid"
   echo " 2 = run lots of 3D test cases (different dimensions)"
   echo " p = run several 3D test cases in parallel (2 and 4 cpus)"
   echo " makeref  = generate new reference output, 3D"
   exit
endif

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 32 32 32
   make ; rm -f $refout 
   ./dns < $refin > $refout
  cat $refout
  cd ../testing/3d
  mv reference3d0000.0000.u restart.u
  mv reference3d0000.0000.v restart.v
  mv reference3d0000.0000.w restart.w
  

endif

if ($1 == 1) then
./gridsetup.py 1 1 1 32 32 32

echo "without restart:"
make >& /dev/null ;  rm -f $tmp ; ./dns < $refin > $tmp 
../testing/check.sh $tmp $refout

echo "with restart:"
make >& /dev/null ;  rm -f $tmp ; ./dns < $refrestart > $tmp 
../testing/check.sh $tmp $refout


endif


if ($1 == g) then
./gridsetup.py 1 1 1 32 32 32
make dnsgrid >& /dev/null ;  rm -f $tmp ; ./dnsgrid < $refrestart > $tmp 
../testing/check.sh $tmp $refout

endif


if ($1 == 2) then

./gridsetup.py 1 1 1 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; ./dns < $refrestart > $tmp 
../testing/check.sh $tmp $refout

endif

if ($1 == p) then
#set in = $refrestart
#echo USING RESTART
set in = $refin
echo NOT USING RESTART

./gridsetup.py 1 1 2 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dns < $in > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 2 1 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dns < $in > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 2 1 1 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dns < $in > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 2 1 2 32 32 32 2 3 4 4 3 2 
make >& /dev/null ;  rm -f $tmp ; mpirun -np 4 ./dns < $in > $tmp 
../testing/check.sh $tmp $refout


endif





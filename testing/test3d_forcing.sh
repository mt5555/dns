#/bin/csh -f

cd ../src
set refin=../testing/reference3d_f.inp
set refout=../testing/reference3d_f.out
set rundir=../testing/3df
set tmp=/tmp/temp.out

set MPIRUN  = "mpirun -np"
if (`uname` == OSF1) then
   set MPIRUN =  "prun -n"
endif


if ($#argv == 0 ) then
   echo "./test3d_forcing.sh [1,2,p,...]  [dns,dnsp,...]"
   echo 
   echo 'to add options:  ./test3d_forcing.sh 2 "dnsp -v3"'
   echo 
   echo " 1 = run dns without restart, simple 3D test case"
   echo " r = run dns with restart, simple 3D test case"
   echo " s = run dns with passive scalars, with and without restart, simple 3D test case"
   echo " 2 = run lots of 3D test cases (different dimensions)"
   echo " p  = run several 3D test cases in parallel (2 and 4 cpus)"
   echo " pr = run several 3D test cases in parallel, with restart"
   echo " pz = run several 3D test cases in parallel, with compressed restart"
   echo " ps = run several 3D test cases in parallel, with spec  restart"

   echo " makeref  = generate new reference output, 3D"
   exit
endif

set code = dns
if ($#argv == 2) then
set code = "$2"
endif
echo "==============================================================="
echo "Running code: " $code
echo "==============================================================="

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 32 32 32
   make $code; rm -f $refout 
   ./$code -d $rundir reference3d  -i $refin | tee $refout
   ./$code -s -d $rundir reference3ds  -i $refin | tee $refout
   ./$code -zo -d $rundir reference3dz  -i $refin | tee $refout

  cd $rundir
  mv reference3d0000.0000.u restart.u
  mv reference3d0000.0000.v restart.v
  mv reference3d0000.0000.w restart.w

  mv reference3ds0000.0000.us restart.us
  mv reference3ds0000.0000.vs restart.vs
  mv reference3ds0000.0000.ws restart.ws

  mv reference3dz0000.0000.uc restart.uc
  mv reference3dz0000.0000.vc restart.vc
  mv reference3dz0000.0000.wc restart.wc

endif

rm -f dns


if ($1 == 1) then
./gridsetup.py 1 1 1 32 32 32

echo "***********************************************************"
echo "without restart:"
make $code >& /dev/null ;  rm -f $tmp ; ./$code -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout
endif


if ($1 == r) then
./gridsetup.py 1 1 1 32 32 32

echo "***********************************************************"
echo "with restart:"
make $code >& /dev/null ;  rm -f $tmp ; ./$code -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
echo "with spectral restart:"
make $code >& /dev/null ;  rm -f $tmp ; ./$code -s -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

endif




if ($1 == 2) then

./gridsetup.py 1 1 1 32 32 32 2 2 0
make $code >& /dev/null ;  rm -f $tmp ; ./$code -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 1 1 32 32 32 2 3 4 4 3 2 
make $code >& /dev/null ;  rm -f $tmp ; ./$code -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout


endif






if ($1 == s) then

echo "***********************************************************"
echo "with 3 passive scalars"
./gridsetup.py 1 1 1 32 32 32 2 2 0 0 0 0 6
make $code >& /dev/null ;  rm -f $tmp ; ./$code  -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
echo "with 3 passive scalars and restart"
./gridsetup.py 1 1 1 32 32 32 2 2 0 0 0 0 6
make $code >& /dev/null ;  rm -f $tmp ; ./$code -r  -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

endif









if ($1 == p || $1 == pr || $1 == pz || $1 == ps) then

if ($1 == pr) then
   set opt = "-r"
   echo USING RESTART
else if ($1 == ps) then
   set opt = "-s -r" 
   echo USING SPEC RESTART   
else if ($1 == pz) then
   set opt = "-r -zi" 
   echo USING COMPRESSED RESTART   
else
   set opt = ""
   echo NOT USING RESTART
endif
if (`uname` == OSF1) then
   set opt = "$opt -b -mio "
endif
echo command line options:  $opt


echo "***********************************************************"
./gridsetup.py 1 1 2 32 32 32 2 2 0
make $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 1 4 32 32 32 2 3 4 4 3 2 
make $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./$code $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 2 1 32 32 32 2 2 0
make $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 1 32 32 32 2 2 0
make $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 2 32 32 32 2 3 4 4 3 2 
make $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./$code $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout


endif





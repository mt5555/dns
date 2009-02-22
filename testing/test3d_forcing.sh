#/bin/csh -f

cd ../src
set refin=../testing/reference3d_f.inp
set refout=../testing/reference3d_f.out
set rundir=../testing/3df
set tmp=/tmp/temp.out



if ($#argv == 0 ) then
   echo "./test3d_forcing.sh [1,2,p,...]  [dns,dnsp,...] [options]"
   echo 
   echo 
   echo " 1 = run dns without restart, simple 3D test case"
   echo " r = run dns with restart, simple 3D test case"
   echo " s = run dns with passive scalars, with and without restart, simple 3D test case"
   echo " 2 = run lots of 3D test cases (different dimensions)"
   echo " p  = run several 3D test cases in parallel (2 and 4 cpus)"
   echo " pr = run several 3D test cases in parallel, with restart"
   echo " prbig = run several 3D test cases in parallel, with restart, 64 cpus"
   echo " pz = run several 3D test cases in parallel, with compressed restart"
   echo " ps = run several 3D test cases in parallel, with spec  restart"

   echo " makeref  = generate new reference output, 3D"
   exit
endif

set code = dns
if ($#argv >= 2) then
  set code = "$2"
endif
set opt = "-d $rundir -i $refin"
if ($#argv >= 3) then
  set opt = "$3 $4 $5 $6 $opt"
endif

set MPIRUN  = "mpirun -np"
if (`uname` == OSF1) then
   set MPIRUN =  "prun -n"
endif
if (`uname` == Darwin) then
   set MPIRUN =  "mpiexec -np"
   set opt = "-b $opt"
endif
if (`arch` == x86_64) then
   set MPIRUN =  "yod -VN -sz "
endif


echo "==============================================================="
echo "Running code: " $code   
echo "MPIRUN = " $MPIRUN
echo "with options:" $opt
echo "==============================================================="

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 32 32 32
   make -j2 $code; rm -f $refout 
   ./$code $opt   reference3d   | tee $refout
   ./$code $opt -s  reference3ds   | tee $refout
   ./$code $opt -zo  reference3dz   | tee $refout

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
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout
endif


if ($1 == r) then
./gridsetup.py 1 1 1 32 32 32

echo "***********************************************************"
echo "with restart:"
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt -r  reference3d   >& $tmp 
../testing/check.sh $tmp $refout



echo "***********************************************************"
echo "with spectral restart:"
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt -s -r  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

endif




if ($1 == 2) then

./gridsetup.py 1 1 1 32 32 32 2 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt -r  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 1 1 32 32 32 2 3 4 4 3 2 
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt -r  reference3d   >& $tmp 
../testing/check.sh $tmp $refout


endif






if ($1 == s) then

echo "***********************************************************"
echo "with 3 passive scalars"
./gridsetup.py 1 1 1 32 32 32 2 2 0 0 0 0 6
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt   reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
echo "with 3 passive scalars and restart"
./gridsetup.py 1 1 1 32 32 32 2 2 0 0 0 0 6
make -j2 $code >& /dev/null ;  rm -f $tmp ; ./$code $opt -r   reference3d   >& $tmp 
../testing/check.sh $tmp $refout

endif





if ($1 == prbig) then
   set opt = "-r $opt"
   echo USING RESTART

echo "***********************************************************"
./gridsetup.py 1 1 32 32 32 32 2 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 32 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 1 32 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 32 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout


echo "***********************************************************"
./gridsetup.py 2 1 32 32 32 32 2 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 64 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 32 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 64 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout


echo "***********************************************************"
./gridsetup.py 2 1 16 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 32 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 4 1 8 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 32 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout


echo "***********************************************************"
./gridsetup.py 1 1 16 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 16 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 16 1 1 32 32 32 2 3 4 4 3 2 
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 16 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout



endif



if ($1 == p || $1 == pr || $1 == pz || $1 == ps) then

if ($1 == pr) then
   set opt = "-r $opt"
   echo USING RESTART
else if ($1 == ps) then
   set opt = "-s -r $opt" 
   echo USING SPEC RESTART   
else if ($1 == pz) then
   set opt = "-r -zi $opt" 
   echo USING COMPRESSED RESTART   
   echo REQUIRES 1 x 1 x N parallel decompostion 
else
   echo NOT USING RESTART
endif
if (`uname` == OSF1) then
   set opt = "$opt -b -mio $opt"
endif
echo command line options:  $opt


echo "***********************************************************"
./gridsetup.py 1 1 2 32 32 32 0 0 2
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 1 2 32 32 32 2 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 1 2 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 1 4 32 32 32 2 3 4 4 3 2 
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 2 1 32 32 32 2 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 1 32 32 32 2 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt   reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 1 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./$code $opt   reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 2 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./$code $opt   reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 2 32 32 32 0 2 0
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./$code $opt -nov3  reference3d   >& $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 2 32 32 32 2 3 4 4 3 2 
make -j2 $code >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./$code $opt  reference3d   >& $tmp 
../testing/check.sh $tmp $refout


endif





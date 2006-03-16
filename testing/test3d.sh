#/bin/csh -f

cd ../src
set refin=../testing/reference3d.inp
set refout=../testing/reference3d.out
set rundir=../testing/3d
set tmp=/tmp/temp.out

set EXE =  ./dns
set EXEG = ./dnsgrid
set MPIRUN  = "mpirun -np"
#   set MPIRUN =  "prun -n"
if (`hostname` == mauve.lanl.gov ) then
   set EXE = "mpirun -np 1 ./dns"
   set EXEG = "mpirun -np 1 ./dnsgrid"
endif



if ($#argv == 0 ) then
   echo "./test3d.sh [1,2,p]"
   echo " 1 = run dns and dnsgrid with and without restart, simple 3D test case"
   echo " s = run dns with passive scalars, with and without restart, simple 3D test case"
   echo " 2 = run lots of 3D test cases (different dimensions)"
   echo " p  = run several 3D test cases in parallel (2 and 4 cpus)"
   echo " pr = run several 3D test cases in parallel, with restart"
   echo " pudm = run several 3D test cases in parallel, with udm restart"
   echo " ps = run several 3D test cases in parallel, with spec  restart"

   echo " makeref  = generate new reference output, 3D"
   exit
endif

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 32 32 32
   make ; rm -f $refout 
   ./dns -d $rundir reference3d  -i $refin > $refout
   ./dns -s -d $rundir reference3ds  -i $refin > $refout
  cat $refout
  cd $rundir
  mv reference3d0000.0000.u restart.u
  mv reference3d0000.0000.v restart.v
  mv reference3d0000.0000.w restart.w
  mv reference3ds0000.0000.us restart.us
  mv reference3ds0000.0000.vs restart.vs
  mv reference3ds0000.0000.ws restart.ws

endif

rm -f dns dnsgrid


if ($1 == 1) then
./gridsetup.py 1 1 1 32 32 32

echo "***********************************************************"
echo "without restart:"
make >& /dev/null ;  rm -f $tmp ; $EXE -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
echo "dnsgrid without restart:"
./gridsetup.py 1 1 1 32 32 32
make dnsgrid >& /dev/null ;  rm -f $tmp ; $EXEG  -d $rundir reference3d -i $refin > $tmp 
../testing/check.sh $tmp $refout


echo "***********************************************************"
echo "with restart:"
make >& /dev/null ;  rm -f $tmp ; $EXE -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
echo "with spectral restart:"
make >& /dev/null ;  rm -f $tmp ; $EXE -s -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout



echo "***********************************************************"
echo "dnsgrid with restart:"
./gridsetup.py 1 1 1 32 32 32
make dnsgrid >& /dev/null ;  rm -f $tmp ; $EXEG -r -d $rundir reference3d -i $refin > $tmp 
../testing/check.sh $tmp $refout




endif




if ($1 == 2) then

./gridsetup.py 1 1 1 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; $EXE -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 1 1 32 32 32 2 3 4 4 3 2 
make >& /dev/null ;  rm -f $tmp ; $EXE -r -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout


endif






if ($1 == s) then

echo "***********************************************************"
echo "with 3 passive scalars"
./gridsetup.py 1 1 1 32 32 32 2 2 0 0 0 0 6
make >& /dev/null ;  rm -f $tmp ; $EXE  -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
echo "with 3 passive scalars and restart"
./gridsetup.py 1 1 1 32 32 32 2 2 0 0 0 0 6
make >& /dev/null ;  rm -f $tmp ; $EXE -r  -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

endif









if ($1 == p || $1 == pr || $1 == pudm || $1 == ps) then

if ($1 == pr) then
   set opt = "-r"
   echo USING RESTART
else if ($1 == ps) then
   set opt = "-s -r" 
   echo USING SPEC RESTART   
else if ($1 == pudm) then
   set opt = "-r -ui" 
   echo USING UDM RESTART   
else
   set opt = ""
   echo NOT USING RESTART
endif
if (`uname` == OSF1) then
   set opt = "$opt -b -mio "
endif

set opt = "$opt -mio "
echo command line options:  $opt


echo "***********************************************************"
./gridsetup.py 1 1 2 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./dns $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 1 2 1 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./dns $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 1 32 32 32 2 2 0
make >& /dev/null ;  rm -f $tmp ; $MPIRUN 2 ./dns $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout

echo "***********************************************************"
./gridsetup.py 2 1 2 32 32 32 2 3 4 4 3 2 
make >& /dev/null ;  rm -f $tmp ; $MPIRUN 4 ./dns $opt -d $rundir reference3d  -i $refin > $tmp 
../testing/check.sh $tmp $refout


endif





#/bin/csh -f

cd ../src
set refin=../testing/vxpair.inp
set refout=../testing/vxpair.out
set tmp=/tmp/temp.out

if ($#argv == 0 ) then
   echo "./test.sh [1,p,makeref]"
   echo " 1 = run 1 simple 2D test case"
   echo " p = run some parallel cases"
   echo " makeref  = generate new reference output, 2D"
   exit
endif

if ($1 == makeref) then

   ./gridsetup.py 1 1 1 65 65 1 2 2 0 2 2 0
   make dnsvor; rm -f $refout 
   ./dnsvor -d /tmp < $refin > $refout
  cat $refout

endif

if ($1 == 1) then

   ./gridsetup.py 1 1 1 65 65 1 2 2 0 2 2 0
make dnsvor >& /dev/null ;  rm -f $tmp ; ./dnsvor -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout


endif

if ($1 == p) then


./gridsetup.py 2 1 1 65 65 1  2 2 0 2 2 0
make dnsvor >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dnsvor -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 1 2 1 65 65 1  2 2 0 2 2 0
make dnsvor >& /dev/null ;  rm -f $tmp ; mpirun -np 2 ./dnsvor -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout

./gridsetup.py 2 2 1 65 65 1  2 2 0 2 2 0
make dnsvor >& /dev/null ;  rm -f $tmp ; mpirun -np 4 ./dnsvor -d /tmp < $refin > $tmp 
../testing/check.sh $tmp $refout


endif


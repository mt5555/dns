#/bin/csh -f

cd ../src
set refin=benchmark1.inp

#if ($#argv == 0 ) then
#   echo "./bencht.sh ?"
#   exit
#endif

./gridsetup.py 1 1 1 64 64 64
#make >& /dev/null ; 
make dns
cd ../benchmark
../src/dns < $refin 



#/bin/csh -f

cd ../src
set refin=../benchmark/benchmark1.inp
set refout=../benchmark/benchmark.out

#if ($#argv == 0 ) then
#   echo "./bencht.sh ?"
#   exit
#endif

./gridsetup.py 1 1 1 64 64 64
make >& /dev/null ; 
cd ../benchmark
../src/dns < $refin 



#/bin/csh -f

cd ../src
set refin=benchmark1.inp

if ($#argv == 0 ) then
   echo "./bench.sh [96,512]"
   exit
endif

if ($1 == 96) then
   ./gridsetup.py 1 1 1 96 96 96
   make dns
   cd ../benchmark
   ../src/dns < benchmark96.inp
endif

if ($1 == 512) then
   ./gridsetup.py 1 1 1 512 512 2
   make dns
   cd ../benchmark
   ../src/dns < benchmark512.inp
endif



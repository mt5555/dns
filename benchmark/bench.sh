#/bin/csh -f

cd ../src
set refin=benchmark1.inp

if ($#argv == 0 ) then
   echo "./bencht.sh [64,96,256]"
   exit
endif

if ($1 == 64) then
   ./gridsetup.py 1 1 1 64 64 64
   make dns
   cd ../benchmark
   ../src/dns < $refin 
endif

if ($1 == 96) then
   ./gridsetup.py 1 1 1 96 96 96
   make dns
   cd ../benchmark
   ../src/dns < benchmark256.inp
endif

if ($1 == 256) then
   ./gridsetup.py 1 1 1 256 8 256
   make dns
   cd ../benchmark
   ../src/dns < benchmark256.inp

endif

if ($1 == 1024) then
   ./gridsetup.py 1 1 1 1024 1024 1
   make dns
   cd ../benchmark
   ../src/dns < benchmark1024.inp

endif


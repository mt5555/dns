#/bin/csh -f

cd ../src
set refin=benchmark1.inp

if ($#argv != 3) then
   echo " run a NxNxN simulation on 'ncpus' for ndelt timesteps "
   echo "./bench.sh N ncpus  ndelt"
   echo "(to run without mpi, use N=0)"
   exit
endif


set ncpus = $2
set command = ../src/dns
if ($ncpus > 0) then
#   set command = "mpirun.lam -np $ncpus ../src/dns"
   set command = "mpirun -np $ncpus ../src/dns"
   if (`hostname` == brain) then
        set command = "mpirun -np $ncpus -npn 2 ../src/dns"
#       set command = "mpirun -np $ncpus -npn 1 ../src/dns"
   endif
   echo $command
else
   set ncpus = 1
endif


./gridsetup.py 1 1 $ncpus $1 $1 $1  2 2 0
make dns
cd ../benchmark
sed s/NDELT/$3/ step.inp.sed > benchmark.inp
$command -i benchmark.inp


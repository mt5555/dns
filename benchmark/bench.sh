#/bin/csh -f

#cd ../src
set refin=benchmark1.inp
set benchdir=~/dns/benchmark
set exe=`pwd`/dns
echo $exe

if ($#argv != 3) then
   echo " run a NxNxN simulation on 'ncpus' for ndelt timesteps "
   echo "./bench.sh N ncpus  ndelt"
   echo "(to run without mpi, use N=0)"
   exit
endif


set ncpus = $2
set command = $exe
if ($ncpus > 0) then
   #set command = "mpirun.lam -np $ncpus $exe"
   set command = "mpirun -np $ncpus $exe"
   if (`uname` == OSF1) then
      #set command = "prun -n $ncpus $exe"
      set command = "/users/taylorm/liblampi/bin/mpirun -np $ncpus $exe"
      setenv MPI_VERSION lampi
   endif
   if (`hostname` == brain) then
        set command = "mpirun -np $ncpus -npn 2 $exe"
        #set command = "mpirun -np $ncpus -npn 1 $exe"
   endif
   if (`hostname` == node001) then
        set command = "/home/gmpi.pgi/bin/mpirun -np $ncpus  $exe"
   endif
   if (`hostname` == pinkish.lanl.gov) then
        set pinkopts = "-o 16 --gm --nper 2"
        set command = "mpirun $pinkopts -np $ncpus   $exe"
   endif
   echo $command
else
   set ncpus = 1
endif


./gridsetup.py 1 1 $ncpus $1 $1 $1  2 2 0
make dns
cd $benchdir
sed s/NDELT/$3/ step.inp.sed > benchmark.inp

if (`hostname` == pinkish.lanl.gov) then
   set pinkhost = `mpirun $pinkopts -np $ncpus hostname`
   set pinknum = `echo $pinkhost | sed 's/n//' `
   echo $pinkhost $pinknum
   bpsh $pinknum mkdir /tmp/taylorm
   bpcp benchmark.inp {$pinkhost}:/tmp/taylorm
   $command -i /tmp/taylorm/benchmark.inp -d /tmp/taylorm
else
   $command -i benchmark.inp
endif










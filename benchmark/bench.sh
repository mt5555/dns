#/bin/csh -f

#cd ../src
set refin=benchmark1.inp
set benchdir=~/dns/benchmark
if (! -e $benchdir) then
   set benchdir=~/lanl/dns/benchmark
endif
if (! -e $benchdir) then
   echo  $benchdir does not exist
   exit 1
endif
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
      set command = "prun -n $ncpus $exe"
      module load MPI_64bit_Thread_Safe_R12
      #set command = "lampirun -np $ncpus $exe"
      #setenv MPI_VERSION lampi

      module list
   endif
   if (`hostname` == brain) then
        set command = "mpirun -np $ncpus -npn 2 $exe"
        #set command = "mpirun -np $ncpus -npn 1 $exe"
   endif
   if (`hostname` == sulaco) then
        set command = "lampirun -np $ncpus $exe"
   endif
   if (`hostname` == node001) then
        set command = "/home/gmpi.pgi/bin/mpirun -np $ncpus  $exe"
   endif
   if (`hostname` == pinkish.lanl.gov) then
        setenv MPI_VERSION lampi
        if ($MPI_VERSION == lampi) then
            set pinknum = 0
            set ncpus2 = 0
            @ ncpus2 = ($ncpus + 1) / 2
            set echo
	    set nodes = `pushlib -t $pinknum+$ncpus2`
	    pushlib $nodes
            set command = "lampirun -np $ncpus -H $nodes   $exe"
        else
            set pinkopts = "-o 16 --gm --nper 2"
            set command = "mpirun $pinkopts -np $ncpus   $exe"
            set pinkhost = `mpirun $pinkopts -np $ncpus hostname`
            set pinknum = `echo $pinkhost | sed 's/n//' `
        endif
   endif
   echo $command
else
   set ncpus = 1
endif


./gridsetup.py 1 1 $ncpus $1 $1 $1  2 2 0
make dns
cd $benchdir
sed s/NDELT/$3/ step.inp.sed > benchmark.inp

echo $command
if (`hostname` == pinkish.lanl.gov) then
   set nodesc = `echo $nodes | tr "," " " `
   foreach nl ($nodesc)
      bpcp benchmark.inp {$nl}:/tmp/taylorm
   end
   $command -d /tmp/taylorm -i /tmp/taylorm/benchmark.inp 
else
   $command -i benchmark.inp
endif










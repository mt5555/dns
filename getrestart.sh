#! /bin/csh -f
#
# get some restart files from HPSS or specified directory
# using the specified basename, find the file of the form basename????.????.u
# with the newest time stamp (????.????)
#
#  $1 = basename of restart file. 
#  $2 = HPSS   get file off of HPSS
#     = path   get file from specified path
#
# status=0  success
# status=1  failed
#
set name = $1
set fpath = $2


if ($fpath == HPSS) then

   #search HPSS for newest restart file
   set resnamew = `psi ls dns/{$name}\*.w | sort | tail -1`
   if ($resnamew =="") then
      echo "Error finding restart file.  Exit"
      exit 1
   else
      set nametime = `basename $resnamew .w`
      set resnameu = `psi ls  dns/{$nametime}\*.u | sort | tail -1`
      set resnamev = `psi ls  dns/{$nametime}\*.v | sort | tail -1`
      echo "Using restart files: " 
      echo $resnameu
      echo $resnamev
      echo $resnamew
   endif
   set resnameu2 = `basename $resnameu`
   set resnamev2 = `basename $resnamev`
   set resnamew2 = `basename $resnamew`

   # check to see if files are left over from last run: 
   if !(-e $resnameu2) then
      psi get $resnameu
   endif
   if !(-e $resnamev2) then
      psi get $resnamev
   endif
   if !(-e $resnamew2) then
      psi get $resnamew
   endif

   \rm -f restart.*
   \ln -s $resnameu2  restart.u
   \ln -s $resnamev2  restart.v
   \ln -s $resnamew2  restart.w
   if !(-e restart.u) then
      echo "No restart.w file"
      exit 1
   endif

else

   #search $fpath for newest restart file
   set resnamew = `\ls {$fpath}/{$name}\*.w | sort | tail -1`
   if ($resnamew =="") then
      echo "Error finding restart file.  Exit"
      exit 1
   else
      set nametime = `basename $resnamew .w`
      set resnameu = `ls  {$fpath}/{$nametime}\*.u | sort | tail -1`
      set resnamev = `ls  {$fpath}/{$nametime}\*.v | sort | tail -1`
      echo "Using restart files: " 
      echo $resnameu
      echo $resnamev
      echo $resnamew
   endif
   set resnameu2 = `basename $resnameu`
   set resnamev2 = `basename $resnamev`
   set resnamew2 = `basename $resnamew`

   # check to see if files are left over from last run: 
   if !(-e $resnameu2) then
      \cp $resnameu .
   endif
   if !(-e $resnamev2) then
      \cp $resnamev .
   endif
   if !(-e $resnamew2) then
      \cp $resnamew .
   endif

   \rm -f restart.*
   \ln -s $resnameu2  restart.u
   \ln -s $resnamev2  restart.v
   \ln -s $resnamew2  restart.w
   if !(-e restart.u) then
      echo "No restart.w file"
      exit 1
   endif

endif

exit 0












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
   set resnamew = `psi ls dns/{$name}/{$name}\*.h5 | sort | tail -1`
   if ($resnamew =="") then
      echo "Error finding restart file.  Exit"
      exit 1
   else
      set nametime = `basename $resnamew .w`
      echo "Using restart file: " 
      echo $resnamew
   endif
   set resnamew2 = `basename $resnamew`

   # check to see if files are left over from last run: 
   if !(-e $resnamew2) then
      psi get $resnamew
   endif

   \rm -f restart.*
   \ln -s $resnamew2  restart.h5
   if !(-e restart.h5) then
      echo "No restart.h5 file"
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
      echo "Using restart file: " 
      echo $resnamew
   endif

   \rm -f restart.*
   \ln -s $resnamew  restart.h5
   if !(-e restart.h5) then
      echo "No restart.h5 file"
      exit 1
   endif

endif

exit 0












#! /bin/csh -f
#
#save data to HPSS
#  $1 = run name
#  $2 = u,v,w,h5  or diag    if present, save only the .u, .v or .w files
#                            (call 3 times, from 3 jobs, to save 3x faster) 
#

set name = $1
set ext = all
if ($#argv >= 2 ) then
   set ext = $2
endif


psi <<EOF
cd dns
mkdir $name
EOF

if ( ( $ext == all ) || ( $ext == diag ) ) then
   echo 'saving diagnostics'
   psi save -d dns/$name $name*.scalars-turb
   psi save -d dns/$name $name*.isostr $name*.isow2s2 $name*.iso1
   psi save -d dns/$name  $name*.sf
   psi save -d dns/$name $name*.s2v2
   psi save -d dns/$name $name*.jpdf
   psi store -d dns/$name  $name*.scalars
   psi store -d dns/$name  $name*.spec
   psi store -d dns/$name $name*.spect
   if ( $ext == all ) then
      echo 'saving uvw'
      psi save -d dns/$name $name*.h5 $name*.u  $name*.v   $name*.w
   endif
else
   echo 'saving .' $ext ' files'
   #psi save -d dns/$name {$name}0003.6781.$ext $name*.$ext
   psi save -d dns/$name $name*.$ext
endif




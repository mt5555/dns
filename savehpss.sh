#! /bin/csh -f
#
#save data to HPSS
#  
#

set name = $1

psi <<EOF
cd dns
mkdir $name
cd $name
save $name*.scalars-turb
save $name*.isostr $name*.isow2s2 $name*.iso1
save $name*.sf
save $name*.s2v2
save $name*.jpdf
store $name*.scalars
store $name*.spec
store $name*.spect
save $name*.h5 $name*.u $name*.v $name*.w
EOF




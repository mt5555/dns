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
save $name*.spec
save $name*.spect
save $name*.isostr
save $name*.sf
save $name*.s2v2
save $name*.jpdf

store $name*.scalars

save $name*.h5
save $name*.u
save $name*.v
save $name*.w
EOF


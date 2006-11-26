#!/bin/tcsh -f 

if ($1 == 0 ) then
   exit 0
endif

set SRC = $2
set code = $3
set mesh = "$4"
set EXE = $5

cd $SRC
rm -f $code
rm -f $EXE
./gridsetup.py $mesh 2 2 0 

make $code 

\cp -f $code $EXE



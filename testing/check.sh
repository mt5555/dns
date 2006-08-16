#/bin/csh -f


grep -a -E 'min/max|(u,v,w)|ke|Ea|w_xx' $1 > /tmp/new.out
grep -a -E 'min/max|(u,v,w)|ke|Ea|w_xx' $2 > /tmp/ref.out
diff -b -B -d /tmp/new.out /tmp/ref.out

if ($status == 0) then
   echo "output is identical"
endif




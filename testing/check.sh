#/bin/csh -f


grep -E 'min/max|(u,v,w)|ke|Ea|w_xx' $1 > /tmp/new.out
grep -E 'min/max|(u,v,w)|ke|Ea|w_xx' $2 > /tmp/ref.out
diff /tmp/new.out /tmp/ref.out

if ($status == 0) then
   echo "output is identical"
endif




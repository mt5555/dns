#/bin/csh -f


grep -E '(u,v,w)|ke' $1 > /tmp/new.out
grep -E '(u,v,w)|ke' $2 > /tmp/ref.out
diff /tmp/new.out /tmp/ref.out

if ($status == 0) then
   echo "output is identical"
endif




%
%########################################################################
%#  read in crossing times, plot
%########################################################################

%name='/ccs/scratch/taylorm/vxpair/vx4096b.cross';
name='/ccs/taylorm/dns/src/vxpair/vx6144d0000.0000.cross';

fid=endianopen(name,'r');

label=[];
time=[];
while (1) 
  [nc,count]=fread(fid,1,'float64');
  if (count~=1) break; end; 
  if (nc>0)
     data=fread(fid,nc,'float64');  % particle lable
     label=[label;data];

     data=fread(fid,nc,'float64');      % crossing time
     time=[time;data];
  end       
end
fclose(fid)

% cross(1,:) = tracer number  (1 through 10)
% cross(2,:) = crossing time.  

maxl=max(label);
minl=min(label);

for i=minl:maxl
  index=find(label==i);
  ti=times(index);
  l=length(ti);
  ti=ti(2:l)-ti(1:l-1);
  ti
  if (length(ti)>0) 
     plot(ti,'o-'); hold on; 
  end       
end
hold off;

  







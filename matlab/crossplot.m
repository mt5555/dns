%
%########################################################################
%#  read in crossing times, plot
%########################################################################

%name='/ccs/scratch/taylorm/vxpair/vx4096b.cross';
%name='/ccs/taylorm/dns/src/vxpair/vx6144d0000.0000.cross';
name='/ccs/taylorm/dns/src/vxpair/vx4096d0000.0000.cross';
fid=endianopen(name,'r');

label=[];
time=[];
ycoord=[];
while (1) 
  [nc,count]=fread(fid,1,'float64');
  if (count~=1) break; end; 
  if (nc>0)
     data=fread(fid,nc,'float64');  % particle lable
     label=[label;data];
     data=fread(fid,nc,'float64');      % crossing time
     time=[time;data];
     data=fread(fid,nc,'float64');      % crossing time
     ycoord=[ycoord;data];
  end       
end
fclose(fid)


% cross(1,:) = tracer number  (1 through 10)
% cross(2,:) = crossing time.  

maxl=max(label);
minl=min(label);

ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  

subplot(1,1,1)
minl=500;
for i=minl:maxl
  index=find(label==i);
  ti=time(index);
  y=ycoord(index);

  ti0=[0;ti];
  l=length(ti0);
  if (length(ti0)>1) 
     delt=ti0(2:l)-ti0(1:l-1);
     tc=ti0(2:l);
     plot(tc,delt,[ccol(i-minl+1),'o-']); hold on; 
  end       
end
hold off;

  







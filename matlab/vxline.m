%
%########################################################################
%#  plotting the .vxline data
%########################################################################
clear;

%ts=input('time=? ');



name='/home/taylorm/ccs/dns/src/vxpair/vx4096c';
times=[0:.1:10];


ccol=[ 'b','g','r','c','m','y', 'b','g','r','c','m','y' ];  

emode=[];
timev=[];
k=0;

figure(1); clf;
subplot(2,1,1)
hold on;
subplot(2,1,2)
hold on;

for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.vxline'];
  fid=endianopen(fname,'r');

  if (fid>=0) 
     disp(fname)

     [nline,count]=fread(fid,1,'float64');
     if (count~=1) break; end;
     time=fread(fid,1,'float64');
     xc=fread(fid,1,'float64');
     yc=fread(fid,1,'float64');
     yline=fread(fid,nline,'float64');
     wline=fread(fid,nline,'float64');
     k=k+1;
     timev(k)=time;
     subplot(2,1,1)
     for i=2:9
        if (wline(i)>=0) 
           plot(time,wline(i),['.',ccol(i)]);
        end
     end
     ylabel('vorticity')

     subplot(2,1,2)
     plot(time,xc,'.')
     plot(time,yc,'.')
     xlabel('time')
     ylabel('x, y coordinates')
  else
     disp(sprintf('error openting file %s',fname))
  end

end

hold off;
      

return



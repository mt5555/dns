%
%########################################################################
%#  plotting ellipses
%########################################################################

%ts=input('time=? ');

%name='../src/temp';
%name='../src/temp0000.0000.ellipse';
%name='../src/vxpair/vx4096b0009.0000.ellipse';
%name='/home/taylorm/vxpair/vx2048a0050.0000.ellipse';
%name='/data/vxpair/vx2048c0000.0000.ellipse';

name='/ccs/taylorm/dns/src/vxpair/vx6144c';
times=[0:.1:19];


for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.ellipse'];
  disp(fname)
  % cant use endianopen() because first number is time
  fid=fopen(fname,'r','l');
  if (fid>=0) 

     [time,count]=fread(fid,1,'float64');
     if (count~=1) break; end;
     [nell,count]=fread(fid,1,'float64');
     np=fread(fid,1,'float64');

     figure(1); clf;   hold on

     for i=1:nell
        x=fread(fid,np,'float64');
        y=fread(fid,np,'float64');
        x(np+1)=x(1);
        y(np+1)=y(1);
        plot(x,y)
       %plot(x,y,'r.')
     end     
     axis([1 3 0 1.5]);
     axis equal
     title(sprintf('time=%f',time))
     hold off;
     %'pause' ;     pause
     fclose(fid);
  end
end

return



%
%########################################################################
%#  plotting ellipses
%########################################################################

%ts=input('time=? ');

%name='../src/temp';
%name='../src/temp0000.0000.ellipse';
name='../src/vxpair/vx2048a0010.0000.ellipse';
fid=fopen(name,'r');

  
while(1)

  [time,count]=fread(fid,1,'float64')
  if (count~=1) break; end;
  
  nell=fread(fid,1,'float64')
  np=fread(fid,1,'float64')
  figure(1)
  clf
  hold on
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
  %    'pause' ;     pause
end
fclose(fid);
return



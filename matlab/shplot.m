
range=.0:.2:75;
range=.0:.2:75;
name='../src/temp'


mkpr=0;            % make ps and jpeg files
mkcontour=1;       % use pcolor or contour


s=findstr(name,'/');
s=s(length(s));
shortname=name(s+1:length(name));

for i=range
  ts=i;
  ts = sprintf('%9.5f',10000+ts);
  ts=ts(2:10);

  tsf=[name,ts,'.vor']
  [x,y,z,q,time]=getfield(tsf);
  nx=length(x);
  ny=length(y);
  nz=length(z);
  q = reshape(q,nx,ny,nz);
  qmax=max(max(max(abs(q))));
  disp(sprintf('max vor=                %f ',qmax));
  vor = squeeze(q(:,:,1))';

  
  tsf=[name,ts,'.h']
  [x,y,z,h,time]=getfield(tsf);
  
  h = reshape(h,nx,ny,nz);
  h = squeeze(h(:,:,1))';
  

  %
  %2D field.  options set above:
  %mkcontour=0,1    use pcolor, or contour plot
  %
  figure(1)

  time=time*2*pi*14;
  subplot(2,1,1)
  ts=sprintf('vorticity  time=%.2f  max=%e',time,qmax)
  pcolor(x,y,vor)
  shading interp
  title(ts);
  axis square

  subplot(2,1,2)
  pcolor(x,y,h)
  shading interp
  axis square
  title('height')
  
  
  disp('pause...')
  pause
end


  
if (mkpr) 
  if (mplot>0)
    orient tall
  end
  pname=[name,'.vor.ps'];
  disp('creating ps...')
  print('-depsc',pname);
  pname=[name,'.vor.jpg'];
  disp('creating jpg...')
  print('-djpeg','-r 96',pname);
end


fidu=fopen('test-0-0-0-0000.0000.data');

%
%########################################################################
%#  plotting output file
%########################################################################

%ts=input('time=? ');

range=50:10:150.00;
%range=29:1.0:1000.0;
%range=[.00];
%name='/ccs/taylorm/dns/src/vxpair/vx6144b_';
%name='../src/temp';
name='/data/vxpair/vx2048a';

usefig=1;
mkpr=1;            % make ps and jpeg files
mkcontour=1;       % use pcolor or contour


s=findstr(name,'/');
s=s(length(s));
shortname=name(s+1:length(name));

for i=range
  ts=i;
  ts = sprintf('%9.5f',10000+ts);
  ts=ts(2:10);

  fname=[name,ts,'.vor']
  fidvor=fopen(fname,'r');
  time=fread(fidvor,1,'float64')
  data=fread(fidvor,3,'float64');
  nx=data(1);
  ny=data(2);
  nz=data(3);
  
  x=fread(fidvor,nx,'float64');
  y=fread(fidvor,ny,'float64');
  z=fread(fidvor,nz,'float64');
  
  q = fread(fidvor,nx*ny*nz,'float64');
  tmp = fread(fidvor,1,'float64');
  tmp=size(tmp);
  if (tmp(1)~=0) 
    disp('Error reading input file...')
  end
  fclose(fidvor);
  q = reshape(q,nx,ny,nz);
  qmax=max(max(max(q)));
  vor = squeeze(q(:,:,1));
  disp(sprintf('max vor=                %f ',qmax));

    
  
  fname=[name,ts,'.psi']
  fidvor=fopen(fname,'r');
  time=fread(fidvor,1,'float64')
  data=fread(fidvor,3,'float64');
  nx=data(1);
  ny=data(2);
  nz=data(3);
  
  x=fread(fidvor,nx,'float64');
  y=fread(fidvor,ny,'float64');
  z=fread(fidvor,nz,'float64');
  
  q = fread(fidvor,nx*ny*nz,'float64');
  tmp = fread(fidvor,1,'float64');
  tmp=size(tmp);
  if (tmp(1)~=0) 
    disp('Error reading input file...')
  end
  fclose(fidvor);
  q = reshape(q,nx,ny,nz);
  psi = squeeze(q(:,:,1));

  

  stitle=sprintf('%s    time=%.2f  max=%f',shortname,time,qmax)
  figure(usefig)
  
    %
    %  2D field.  options set above:
    %  mkcontour=0,1    use pcolor, or contour plot
    %
    clf
    subplot(2,1,1)

    if (mkcontour==0)
      pcolor(x,y,vor')
      shading interp
    else
      v = -12:1:3;
      v=2.^v;
      %v = 4.5*[3/4, 1/2, .4 , 1/4, 1/8]
      contour(x,y,vor',v)
      hold on
      contour(x,y,vor',[.005 .005],'r')
      hold off
    end
    
    title(stitle);
    axis equal
    axis([0,max(x),0,max(y)]);


    subplot(2,1,2)

    if (mkcontour==0)
      pcolor(x,y,psi')
      shading interp
    else
      v=20;                             % use 20 contours
      contour(x,y,psi',v)
    end
    axis equal
    axis([0,max(x),0,max(y)]);

    
    

  
    if (mkpr) 
      orient tall
      pname=[name,'.vor.ps'];
      disp('creating ps...')
      print('-depsc',pname);
      pname=[name,'.vor.jpg'];
      disp('creating jpg...')
      print('-djpeg','-r 96',pname);
    end

    'pause'
    pause
end
return




%
%########################################################################
%#  plotting output file
%########################################################################


%ts=input('time=? ');

range=1:.5:4;
%range=50:5:150.00;
%range=29:1.0:1000.0;
%range=[125:10:200];
%name='../src/vxpair/vx6144b_';
%name='/ccs/taylorm/dns/src/vxpair/vx4096d';
%name='/ccs/scratch/taylorm/vxpair/vx6144e';
%name='/scratch2/taylorm/vx12288b/vx12288b';
%name='/ccs/taylorm/dns/src/vxpair/vx4500a';
%name='/ccs/taylorm/dns/src/vxpair/vx4500a';
%name='/ccs/scratch/taylorm/kras/vx2560a/vx2560a';
%name='/ccs/scratch/taylorm/kras/vx2560c/vx2560c';
%name='/ccs/scratch/taylorm/kras/vx5120a/vx5120a';
name='/home/scratch/kras/vx7680a/vx7680a'
%name='/home/scratch/kras/vx1280a/vx1280a';


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

  % get the data size:
  [x,y,z]=getfield(fname);
  nx=length(x);
  ny=length(y);

  % whole plot, but subsample data:
  % subsample=1;
  % skip=4;   x1=1; x2=nx; y1=1; y2=ny;

  % zoom in:  
  subsample=1;   skip = 1;
  x1=floor(nx/4); x2=floor(3*nx/4);
  y1=floor(ny/4); y2=floor(3*ny/4);

  rangex=x1:skip:x2;
  rangey=y1:skip:y2; 


  [x,y,z,vor,time]=getfield(fname);
  qmax=max(max(max(vor)));
  vor = squeeze(vor(:,:,1));
  if (subsample==1) vor=vor(rangex,rangey); end;
  disp(sprintf('max vor=                %f ',qmax));
    

  fname=[name,ts,'.psi']
  [x,y,z,psi,time]=getfield(fname);
  psi = squeeze(psi(:,:,1));
  if (subsample==1) psi=psi(rangex,rangey); end;


  if (subsample==1) 
    x=x(rangex);
    y=y(rangey);
  end

  

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
      v = -7:2:5;
      v=2.^v;
      %v = 4.5*[3/4, 1/2, .4 , 1/4, 1/8]
      %v=  [3.3938    2.2626    1.1313    0.5656];
      contour(x,y,vor',v,'b')
      hold on
      contour(x,y,vor',-v,'g')
      contour(x,y,vor',[.005 .005],'r')
      hold off
    end
    
    title(stitle);
    axis equal
    %axis([0,max(x),0,max(y)]);


    subplot(2,1,2)

    if (mkcontour==0)
      pcolor(x,y,psi')
      shading interp
    else
      v=20;                             % use 20 contours
      contour(x,y,psi',v)
    end
    axis equal
    %axis([0,max(x),0,max(y)]);

    
    

  
    if (mkpr) 
      orient tall
      pname=[name,'.vor.ps'];
      disp('creating ps...')
      print('-depsc','-r 300',pname);
      pname=[name,'.vor.jpg'];
      disp('creating jpg...')
      print('-djpeg','-r 200',pname);
    end

    'pause'
    pause
end
return




%
%########################################################################
%#  plotting output file
%########################################################################


%ts=input('time=? ');

range=1:1:4;
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
name='/ccs/scratch/taylorm/kras/vx5120a/vx5120a';
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
  [x,y,z,vor,time]=getfield(fname);
  qmax=max(max(max(vor)));

  vor = squeeze(vor(:,:,1));
  disp(sprintf('max vor=                %f ',qmax));
    

  fname=[name,ts,'.psi']
  [x,y,z,psi,time]=getfield(fname);
  psi = squeeze(psi(:,:,1));

  subsample=4;
  if (subsample>1) 
    nx=length(x);
    ny=length(y);
    vor=vor(1:subsample:nx,1:subsample:ny);
    psi=psi(1:subsample:nx,1:subsample:ny);
    x=x(1:subsample:nx);
    y=y(1:subsample:ny);
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
      v = -7:1:5;
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



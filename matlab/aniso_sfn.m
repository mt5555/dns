function[]=aniso_sfn(Dl,Dt,Dlnorm,Dtnorm,ndelta,ndir,dir_max,r_val,nx,delx_over_eta,xx)
% this function is called from isoave.m

cdir = ['k','b','g','r','m']; %colors, one for each m-comp.

%cdir=[ 'k',',','k' ];  % x,y,zz
%cdir=[cdir, 'g','g','g','g','g','g'];  % face diagonals
%cdir=[cdir, 'r','r','r','r'];      % body diagonals
%cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,2,0) directions
%cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,1,2) directions
%cdir=[cdir, 'b','b','b','b','b','b','b','b','b','b','b','b'];      % 12 (1,2,2) directions
%cdir=[cdir, 'y','y','y','y','y','y','y','y','y','y','y','y'];      % 12 (1,3,0) directions
%cdir=[cdir, 'y','y','y','y','y','y','y','y','y','y','y','y'];      % 12 (1,1,3) directions

%perform angle-average

if (ndir==3)
  w=ones([1,3])/3;
else
  equalw=0;
  if (equalw) 
    % put this in to use equally weighted:
    w=ones([1,ndir])/ndir;
  else
    % get the weights:
    wname=sprintf('../src/voronoi/isoave.weights%i',ndir);
    w=textread(wname,'%f');
    % take every other weight
    w=2*w(1:2:length(w));
  end
end
if (abs(1-sum(w))>1e-7) 
  disp('error: weights do not sum to 1')
  return;
end

msize=8;   % marker size
  yltave=0*xx;       %avg over m after angle average
yttave=0*xx;         %avg over m after angle average
for j=1:5              % sphere-harm comps m=-2...2
  ylave=0*xx;
ytave=0*xx;
  for dir=1:ndir         % directions
  x = r_val(:,dir);      % units of box length
  x_plot=x*nx*delx_over_eta;  % units of r/eta
  xx_plot = xx*nx*delx_over_eta;

%matlab calculated
yl = Dlnorm(:,dir,j);
yt = Dtnorm(:,dir,j);

ylave = ylave + w(dir)*spline(x,yl,xx);
ytave = ytave + w(dir)*spline(x,yt,xx);
end
yltave = yltave + ylave;
yttave= yttave + ytave;

  figure(10)
  loglog(xx_plot,abs(ylave),['.',cdir(j)],'MarkerSize',msize);   hold on;
set(gca,'fontsize',16)
title('Slt, matlab')
  legend('-2','2','-1','1','0')
figure(11)
  loglog(xx_plot,abs(ytave),['.',cdir(j)],'MarkerSize',msize);   hold on;
set(gca,'fontsize',16)
title('Stt, matlab')
legend('-2','2','-1','1','0')



%  figure(12)
%    loglog(x_plot,abs(yl),['.',cdir(j)],'MarkerSize',msize);   hold on;
%title('Slt, fortran')

%figure(13)
%  loglog(x_plot,yt,['.',cdir(j)],'MarkerSize',msize);   hold on;
%title('Stt, fortran')

end

yltave = yltave/5;
yttave= yttave/5;

figure(10)
 loglog(xx_plot,abs(yltave),['.-','c'],'MarkerSize',msize);  
 
figure(11)
loglog(xx_plot,abs(yttave),['.-','c'],'MarkerSize',msize)

figure(11)
     xs = r_val(:,dir_max)
     yl = Dl(:,dir_max,j)
     yt = Dt(:,dir,j)
     loglog(x_plot, abs(yl));hold on;

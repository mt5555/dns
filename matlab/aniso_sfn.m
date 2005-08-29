function[]=aniso_sfn(Dl,Dt,Dlt1_wt,Dlt2_wt,ndelta,ndir,dir_max,r_val,nx,delx_over_eta,xx)
% this function is called from isoave.m. It computes the mixed structure function in the j=2 speherical harmonic basis. It also computes the mixed structure function in the special direction idir_max in which the shear is maximum.

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

ylt1ave=0*xx;       %avg over m after angle average
ylt2ave=0*xx;         %avg over m after angle average

for j=1:5              % sphere-harm comps m=-2...2

ylt1_d=0*xx;
ylt2_d=0*xx;
for dir=1:ndir         % directions
     x = r_val(:,dir);      % units of box length
     x_plot=x*nx*delx_over_eta;  % units of r/eta
     xx_plot = xx*nx*delx_over_eta;

     ylt1 = Dlt1_wt(:,dir,j);
     ylt2 = Dlt2_wt(:,dir,j);

     ylt1_d = ylt1_d + w(dir)*spline(x,ylt1,xx);
     ylt2_d = ylt2_d + w(dir)*spline(x,ylt2,xx);
end

ylt1ave = ylt1ave + ylt1_d;
ylt2ave= ylt2ave + ylt2_d;

figure(10)
loglog(xx_plot,abs(ylt1_d),['.',cdir(j)],'MarkerSize',msize);   hold on;
set(gca,'fontsize',16)
title('Slt1, matlab')
legend('-2','2','-1','1','0')

figure(11)
loglog(xx_plot,abs(ylt2_d),['.',cdir(j)],'MarkerSize',msize);   hold on;
set(gca,'fontsize',16)
title('Slt2, matlab')
legend('-2','2','-1','1','0')

end

ylt1ave = ylt1ave/5;
ylt2ave= ylt2ave/5;

figure(10)
loglog(xx_plot,abs(ylt1ave),['.-','c'],'MarkerSize',msize);  
legend('-2','2','-1','1','0','avg')
figure(11)
loglog(xx_plot,abs(ylt2ave),['.-','c'],'MarkerSize',msize)
legend('-2','2','-1','1','0','avg')



figure(15)
   for dir = 1:ndir
     xs = r_val(:,dir);
     x_plot=xs*nx*delx_over_eta;
     ylt1 = Dl(:,dir,1);
     ylt2 = Dl(:,dir,2);
     if (dir==dir_max)
     dir
     loglog(x_plot, abs(ylt1),['.-','r'],'MarkerSize',msize);hold on;
else
    loglog(x_plot, abs(ylt1),['.-','b'],'MarkerSize',msize);hold on;
    end
end


figure(16)
     xs = r_val(:,dir_max);
     x_plot=xs*nx*delx_over_eta;
     ylt1 = Dl(:,dir_max,1);
     ylt2 = Dl(:,dir_max,2);
     loglog(x_plot, abs(ylt1),['.-','r'],'MarkerSize',msize);hold on;
    loglog(x_plot, abs(ylt2),['.-','b'],'MarkerSize',msize);hold on;
end


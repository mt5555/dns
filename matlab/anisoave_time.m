%
% Calculation of time and angle-averaged (spherical harmoic weighted) mixed structure function
%S.Kurien
%

mu=0;
ke=0;
nx=1;
delx_over_eta=1;
eta = 1/(nx*delx_over_eta);
ext='.isostr4';
ndir = 73;

%subcubes
nxdecomp=2;
nydecomp=2;
nxdecomp=2;

%name='/auto/nest/u/skurien/dns/src/skhel512a'
name='/nh/u/skurien/projects/shear/skhel512a'
nx=512; delx_over_eta=2.5; epsilon=2.72;teddy = 1.05; %Rl=250

cdir = ['k','b','g','r','m','y','c','k']; %colors, one for each m-comp.
times=[7];

ndir_use = 0;

teddy=1;

spharm=1;


xx=(1:.5:(nx./2.5)) / nx;        %interpolation points
xx_plot = xx*nx*delx_over_eta;    %units of r/neta  
lex = round(length(xx_plot)/4);
lenx = length(xx_plot);

ndir_use=73;
%
% use only 49 directions:
if (ndir_use>0) 
ndir=ndir_use; 
end




if (ndir==3)
  w=ones([1,3])/3;
else
  equalw=0;
  if (equalw) 
    % put this in to use equally weighted:
    w=ones([1,ndir])/ndir;
    disp(sprintf('Using equall weights for spherical integration'))
  else
    % get the weights:
    wname=sprintf('../src/voronoi/isoave.weights%i',ndir);
    disp(sprintf('Reading Voronio weights from file:  %s',wname))
    w=textread(wname,'%f')
     length(w)
    % take every other weight
   w=2*w(1:2:length(w));
  end
end
if (abs(1-sum(w))>1e-7) 
  disp('error: weights do not sum to 1')
  return;
end

%  close all;



%times_plot=[];

ylt1_time = zeros(5,length(xx));  %time average for each m after angle-average
ylt2_time = zeros(5,length(xx));  %time average for each m after angle-average

ylt1_s = 0*x;  %time average of Y_00 after angle-average
ylt2_s = 0*x;  %time average of Y_00 after angle-average

for t=times

  for subz = 0:(nxdecomp-1)
  for suby = 0:(nydecomp-1)
  for subx = 0:(nzdecomp-1)

    if (nxdecomp*nydecomp*nzdecomp>1) 
      ext_sub = sprintf('_%d%d%d',subx,suby,subz)
    else
     ext_sub = [];
    end

   sprintf('%f',t);
   tstr=sprintf('%10.4f',t+10000);
   fname=[name,tstr(2:10)];
   disp([fname,ext]);
   [nx,ndelta,ndir,r_val,ke,epsilon,mu,...
    D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt,...
    SP_lll,SN_lll,SP1_ltt,SP2_ltt,SN1_ltt,SN2_ltt,H_ltt,H_tt,D_lltt,Dl,Dt,...
    h_epsilon] ...
   = readisostr([fname,ext,ext_sub] );

[ux,dir_max] = read_ux([fname,'.ux',ext_sub]);

%
% xx is given in units of the box length
% but r_val is given in units of delx.  convert to box length units
%

r_val=r_val/nx;
xx_plot = xx*nx*delx_over_eta;        % units of r/eta

lambda=sqrt(10*ke*mu/epsilon);       % single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu;


disp(sprintf('KE:      %f  2pi units: %f',ke,ke*4*pi*pi));
disp(sprintf('epsilon: %f  2pi units: %f',epsilon,epsilon*4*pi*pi));
disp(sprintf('h_epsilon: %f  2pi units: %f',h_epsilon,h_epsilon*2*pi));
disp(sprintf('mu:      %f  2pi units: %f',mu,mu*4*pi*pi));
disp(sprintf('eta:     %f  2pi units: %f',eta,eta*2*pi));
disp(sprintf('lambda:  %f  2pi units: %f',lambda,lambda*2*pi));
disp(sprintf('delx/eta %f',delx_over_eta));
disp(sprintf('R_l:     %f',R_lambda));
disp(sprintf('ndir:     %f',ndir));
disp(' ')


[Dlt1_wt, Dlt2_wt]=sphere_harm_weight(Dl(:,:,1),Dl(:,:,2),spharm);


msize=8;   % marker size


for j=1:(2*spharm+1)              % sphere-harm comps -2j-1 ... 2j+1

  ylt1_s=0*xx;  
  ylt2_s=0*xx;
    
  for dir=1:ndir         % directions
     x = r_val(:,dir);      % units of box length
     x_plot=x*nx*delx_over_eta;  % units of r/eta
     xx_plot = xx*nx*delx_over_eta;


%sfn weighted by mixed spherical harmonic Y_1m and angle-averaged
     ylt1 = Dlt1_wt(:,dir,j);
     ylt2 = Dlt2_wt(:,dir,j);

     ylt1_s = ylt1_s + w(dir)*spline(x,ylt1,xx);
     ylt2_s = ylt2_s + w(dir)*spline(x,ylt2,xx);
  

%raw sfn in special direction, un-weighted by spherical harmonic 
  if (dir==dir_max)

     ylt1_d = Dl(:,dir,1);
     ylt2_d = Dl(:,dir,2);
     ind = subx + 2*suby + 4*subz + 1  %index of subcube
       figure(1)
       loglog(x, abs(ylt1_d),['+',cdir(ind)],'MarkerSize',msize);hold on;
       loglog(x, abs(ylt2_d),['--',cdir(ind)],'MarkerSize',msize);hold on;
  end
       
  end

% angle-averaged sfn weighted by mixed spherical harmonic Y_1m, time accumulate 
ylt1_time(j,:) = ylt1_time(j,:) + ylt1_s;
ylt2_time(j,:) = ylt2_time(j,:) + ylt2_s;
 
end


end
end
end

end  %matches for t=times

figure(1)
%legend('S_{LT1},1', 'S_{LT2},1','S_{LT1},2', 'S_{LT2},2','S_{LT1},3', 'S_{LT2},3','S_{LT1},4', 'S_{LT2},4', 'S_{LT1},5', 'S_{LT2},5','S_{LT1},6', 'S_{LT2},6','S_{LT1},7', 'S_{LT2},7','S_{LT1},8', 'S_{LT2},8');
 title(sprintf('Time = %d; subcubes mixed sfn in special directions',t));

%time average
ylt1_time(:,:) = ylt1_time(:,:)/length(times);
ylt2_time(:,:) = ylt2_time(:,:)/length(times);


for j=1:(2*spharm+1)

j
cdir(j)

figure(22)
loglog(xx_plot,abs(ylt1_time(j,:)),['.-',cdir(j)],'MarkerSize',msize);hold on;
legend('m=-1','1','0')
title('Time average S_{LT1} projected on j=1')

figure(23)
loglog(xx_plot,abs(ylt2_time(j,:)),['.-',cdir(j)],'MarkerSize',msize);hold on;
legend('m=-1','1','0')
title('Time average S_{LT2} projected on j=1')

figure(24)
loglog(xx_plot,abs(ylt1_time(j,:)+ ylt2_time(j,:))/2,['.-',cdir(j)],'MarkerSize',msize);hold on; 
legend('m=-1','1','0')
title('Average of S_{LT1} and S_{LT2} projected on j=1')

end

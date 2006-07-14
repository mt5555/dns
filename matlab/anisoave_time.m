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
nxdecomp=1;
nydecomp=1;
nzdecomp=1;

%name='/auto/nest/u/skurien/dns/src/skhel512a'
%name='/nh/u/skurien/projects/shear/skhel512a'
%name='/nh/u/skurien/projects/shear/livescu/512-3b/512-3b'
%name='/nh/u/skurien/projects/shear/livescu/S=1.275/fort.298';
name='/nh/u/skurien/projects/shear/livescu/S=10.2/fort.266';
nx=512; delx_over_eta=1.855; epsilon=5.398;teddy = 1.05; %Rl=250

cdir = ['k','b','g','r','m','c','k--','b--','r--','m--','c--','k.-','b.-','r.-']; %colors, one for each m-comp.
col = 0;
%times=[5:0.2:7.2];

times=[0];

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
     w=textread(wname,'%f');
     length(w);
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


%ylt1_s = 0*x;  %time average of Y_00 after angle-average
%ylt2_s = 0*x;  %time average of Y_00 after angle-average

ylt1_ds_time = 0*xx;
  ylt2_ds_time = 0*xx;
  ylt_ds_time = 0*xx; 


for t=times
col = col+1

  for subz = 0:(nxdecomp-1)
  for suby = 0:(nydecomp-1)
  for subx = 0:(nzdecomp-1)

    index = subx + 2*suby + 4*subz + 1;  % index of subcube
    
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

figure(100)
plot(index,dir_max,'ro'); hold on;
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
  
  ylt1_ds = 0*xx;
  ylt2_ds = 0*xx;
  ylt_ds = 0*xx;  
  
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

%plot str_fns in 73 directions 
ylt1_d = Dl(:,dir,1);
ylt2_d = Dl(:,dir,2);

ylt1_ds = ylt1_ds + w(dir)*spline(x,ylt1_d,xx);
ylt2_ds = ylt2_ds + w(dir)*spline(x,ylt2_d,xx);
ylt_ds = ylt_ds + w(dir)*spline(x,(ylt1_d+ylt2_d),xx);


figure(99)
loglog(x_plot, abs(ylt1_d), cdir(col),'MarkerSize',msize);hold on;
title('|D_{LT2}| in 73 directions')
figure(98)
loglog(x_plot, abs(ylt2_d), cdir(col),'MarkerSize',msize);hold on;
title('|D_{LT3}| in 73 directions')
figure(97)
loglog(x_plot, abs(ylt1_d+ylt2_d), cdir(col),'MarkerSize',msize);hold on;
title('|D_{LT2} + D_{LT3}| in 73 directions')

  if (dir==dir_max)

%  if (dir==2)
     ylt1_d = Dl(:,dir,1);
     ylt2_d = Dl(:,dir,2);

       figure(1)
%       subplot(2,1,1)
       loglog(x_plot, abs(ylt2_d./ylt1_d),cdir(col),'MarkerSize',msize);hold on;
       title('Ratio D_{LT3}/D_{LT2} computed in max shear direction')
       figure(2)
       loglog(x_plot,abs(ylt1_d),cdir(col),'Markersize',msize);hold on;
       title('|D_{LT2}| computed in max shear direction')
       figure(3)
       loglog(x_plot,abs(ylt2_d),cdir(col),'Markersize',msize);hold on;
        title('|D_{LT3}| computed in max shear direction')
     figure(4)
       loglog(x_plot,abs(ylt1_d+ylt2_d),cdir(index),'Markersize',msize);hold on;
       title('|D_{LT}|= |D_{LT2}+D_{LT3}| computed in max shear direction')
       
       
%       loglog(x, (ylt2_d),['--',cdir(index)],'MarkerSize',msize);hold on;
%       title(sprintf('++ +ve S_{LT1}; -- +veS_{LT2}; In special directions'))
%       subplot(2,1,2)
%       loglog(x, -(ylt1_d),['+',cdir(index)],'MarkerSize',msize);hold on;
%       loglog(x, -(ylt2_d),['--',cdir(index)],'MarkerSize',msize);hold on;
%       title(sprintf('++ -ve S_{LT1}; -- -ve S_{LT2}; In special directions'))

     end %matches if (dir==..)
       
     end %matches for dir = 1:ndir

ylt1_ds_time = ylt1_ds_time + ylt1_ds;
ylt2_ds_time = ylt2_ds_time + ylt2_ds;
ylt_ds_time = ylt_ds_time + ylt_ds;

  
% angle-averaged sfn weighted by mixed spherical harmonic Y_1m, time accumulate 
ylt1_time(j,:) = ylt1_time(j,:) + ylt1_s;
ylt2_time(j,:) = ylt2_time(j,:) + ylt2_s;
 
if (0)
  
figure(20+index)
subplot(2,1,1)
loglog(xx_plot,(ylt1_time(j,:)),['.-',cdir(j)],'MarkerSize',msize);hold on;
legend('m=-1','1','0')
title(sprintf('Time average +ve S_{LT1} projected on j=1, subcube= %d',index))
%figure(40+ind)
subplot(2,1,2)
loglog(xx_plot,-(ylt1_time(j,:)),['.-',cdir(j)],'MarkerSize',msize);hold on;
legend('m=-1','1','0')
title(sprintf('Time average -ve S_{LT1} projected on j=1, subcube= %d',index))


figure(40+index)
subplot(2,1,1)
loglog(xx_plot,(ylt2_time(j,:)),['.-',cdir(j)],'MarkerSize',msize);hold on;
legend('m=-1','1','0')
title(sprintf('Time average +ve S_{LT2} projected on j=1, subcube= %d',index))
%figure(43+ind)
subplot(2,1,2)
loglog(xx_plot,-(ylt2_time(j,:)),['.-',cdir(j)],'MarkerSize',msize);hold on;
legend('m=-1','1','0')
title(sprintf('Time average -ve S_{LT2} projected on j=1, subcube= %d',index))


figure(50+index)
subplot(2,1,1)
loglog(xx_plot,(ylt1_time(j,:)+ ylt2_time(j,:))/2,['.-',cdir(j)],'MarkerSize',msize);hold on; 
legend('m=-1','1','0')
title(sprintf('Average +veS_{LT1} and S_{LT2} projected on j=1, scube = %d', index))
subplot(2,1,2)
loglog(xx_plot,-(ylt1_time(j,:)+ ylt2_time(j,:))/2,['.-',cdir(j)],'MarkerSize',msize);hold on; 
legend('m=-1','1','0')
title(sprintf('Average -ve S_{LT1} and S_{LT2} projected on j=1, scube = %d', index))

end

end %matches for j==


end  %matches for subx
end  %matches for suby
end  %matches for subz



end  %matches for t=times


%time average OF LAST SUBCUBE
ylt1_ds_time = ylt1_ds_time/length(times);
ylt2_ds_time = ylt2_ds_time/length(times);
ylt_ds_time = ylt_ds_time/length(times);

figure(95)
loglog(xx_plot,abs(ylt1_ds_time),'r-');hold on;
loglog(xx_plot,abs(ylt2_ds_time),'b-');hold on;
loglog(xx_plot,abs(ylt_ds_time),'m-');hold on;
legend('|D_{LT2}|','|D_{LT3}|','|D_{LT}|')
title('angle averaged')

ylt1_time(:,:) = ylt1_time(:,:)/length(times);
ylt2_time(:,:) = ylt2_time(:,:)/length(times);


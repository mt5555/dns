%
%
% 2/15ths law: time and angle averaged 
%
%
mu=0;
ke=0;
nx=1;
delx_over_eta=1;
eta = 1/(nx*delx_over_eta);

%name = '/home2/skurien/helicity_data/isostr_1/check256_hq_';
%pname = 'check256_hq_';
%ext = '.new.isostr';
%times=[0:1:30];
%nx = 256;

%name = '/home2/skurien/helicity_data/helical_forced/hel256_hpi2/hel256_hpi2_';
%pname = 'hel256\_hpi2\_';
%ext='.new.isostr';
%times=[4.2:0.2:15.8];
%nx = 256;

name = '/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/skhel512a';
pname = 'skhel512a';
ext='.new.isostr';
times=[2.0:0.2:3.8];
nx = 512;

[avg_eps, avg_heps, avg_delx_over_eta] = ensemble_avg_params(name,ext,times)

delx_over_eta=avg_delx_over_eta; epsilon=avg_eps; h_epsilon=avg_heps;  
teddy=1.05;



ndir_use=0;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')

% this type of averging is expensive:
time_and_angle_ave=1;

k=0;


xx=(1:.5:(nx./2.5)) / nx;
xx_plot=(1:.5:(nx./2.5)) *delx_over_eta;   % units of r/eta

y215_ave=zeros([length(xx),73]);
if (time_and_angle_ave==1) 
   y215_iso_ave=zeros([length(xx),1]);
end

times_plot=[];
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10)];
  disp([fname,ext]);
ppname = [pname,tstr(2:10),ext];
  fid=fopen([fname,ext]);
  if (fid<0) ;
    disp('error opening file, skipping...');
  else
    fclose(fid);
    times_plot=[times_plot,t];
    k=k+1;

    if (time_and_angle_ave==1) 
      klaws=3;                          % compute 2/15 laws
      plot_posneg=0;
      check_isotropy=0;
      
      [y45,y415,y43,eps,h_eps,y215]=compisoave(fname,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,0);
      
      
      mx215_iso_localeps(k)=max(y215);
      mx215_iso(k)=max(y215)*h_eps/h_epsilon;
      mn215_iso_localeps(k)=y215(1);
      mn215_iso(k) = y215(1)*h_eps/h_epsilon;


      y215_iso_ave=y215_iso_ave+y215';
      
    end    


    [nx,ndelta,ndir,r_val,ke,eps_l,mu,tmp,tmp,tmp,tmp,tmp,tmp,tmp,...
    tmp,tmp,tmp,tmp,tmp,H_ltt,H_tt,tmp,tmp,tmp,h_eps_l] ...
        = readisostr( [fname,ext] );
    
    eta_l = (mu^3 / eps_l)^.25;
    delx_over_eta_l=(1/nx)/eta_l;

    
    for dir=1:73;
      x=r_val(:,dir)/nx;                % box length
      y=-H_ltt(:,dir)./(x.^2*abs(h_eps_l)/2);
      
      y215 = spline(x,y,xx);
      
      mx215_localeps(k,dir)=max(y215);
      mx215(k,dir)=max(y215)*h_eps_l/h_epsilon;
      
      mn215_localeps(k,dir)=y215(1);
      mn215(k,dir) = y215(1)*h_eps_l/h_epsilon;

      y215_ave(:,dir)=y215_ave(:,dir)+y215';
    end
  end

end

times=times_plot;
y215_ave=y215_ave/length(times);
y215_iso_ave=y215_iso_ave/length(times);

if (0)

% averging starting at t=0:
save k215data_t0 teddy times mx215_localeps mx215_iso_localeps ...
      y215_ave y215_iso_ave xx_plot

% averging starting at t=1:
save k215data_t1 teddy times mx215_localeps mx215_iso_localeps ...
      y215_ave y215_iso_ave xx_plot

% load data from t=0 (averaged wrong!) for isoave paper:
load k215data_t0 

% load data from t=1 (for time averaged) for isoave paper:
load k215data_t1 
end

scale = 1;  % scale = 2/15 if normalizing by prefactor as well

figure(5); clf; hold on; 
for i=1:1:1 %directions

plot(times/teddy,mx215_localeps(1:length(times),i)/scale,'g-','LineWidth',1.0);hold on;
plot(times/teddy,mn215_localeps(1:length(times),i)/scale, 'r-','LineWidth',1.0);hold on;
end
plot(times/teddy,mx215_iso_localeps/scale,'b-','LineWidth',1.0);hold on;
plot(times/teddy,mn215_iso_localeps/scale,'k-','LineWidth',1.0);hold on;
ax=axis;
%axis( [ax(1),ax(2),0,1] );
plot(times,(2/15)*times./times,'m');
title('Green - dir max; Red - dir first pt; Blue - avg max; Black - avg first pt; Magenta - 2/15')   
hold off;

xlabel('time','FontSize',16)
     ylabel(ppname)
     grid on;
print -dpsc k215time.ps


figure(6); clf
scale = 2/15; % scale = 2/15 to scale out
for i=[1:1:73]
  %semilogx(xx_plot,y215_ave(:,i),'k:','LineWidth',2.0); hold on
  semilogx(xx_plot,(y215_ave(:,i))/scale,'r-','LineWidth',1.0); hold on
end
%semilogx(xx_plot,y215_iso_ave,'k','LineWidth',2.0); hold on
semilogx(xx_plot,y215_iso_ave/scale,'b','LineWidth',2.0); hold on
max(y215_iso_ave)
%axis([1 1000 -0.05 0.15])
x=1:1000; plot(x,(2/15)*x./x/scale,'k');
hold off;
title('H_{ltt} /h r^2   (2/15 law) Ensemble averaged');
ylabel(ppname);
xlabel('r/\eta','FontSize',16);

print -dpsc k215mean.ps



figure(7); clf
offset=y215_iso_ave;
stdr=0*offset;
scale = 2/15; % scale = 2/15 if need to factor out the 4/5, 1 otherwise
for i=[1:1:73]
  semilogx(xx_plot,abs(y215_ave(:,i)-offset)/offset,'k-','LineWidth',1.0); hold on
  stdr=stdr+(y215_ave(:,i)-offset).^2;
end
stdr=sqrt(stdr/15)./offset;
%axis([1 1000 -.2 .2])
%x=1:1000; plot(x,(2/15)*x./x,'k');
hold off;
title('timemean(H_{ltt}(dir)-H_{ltt}(avg))/H_{ltt}(avg) Measure of anisotropy ');
ylabel(ppname);
xlabel('r/\eta','FontSize',16);



starttime=1;
ln=find(times>=starttime);
ln=ln(1);
lnmax=length(times);

dir_use=3;
[mean(mx215(ln:lnmax,dir_use)),mean(mx215_iso(ln:lnmax))]
[std(mx215(ln:lnmax,dir_use)),std(mx215_iso(ln:lnmax)) ]








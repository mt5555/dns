%
%
% 4/5 ths law: time and angle averaged 
%
%
mu=0;
ke=0;
nx=1;
delx_over_eta=1;
eta = 1/(nx*delx_over_eta);
ext='.isostr';

%name='/scratch1/taylorm/iso12w512A0001.3847'
%nx=512; delx_over_eta=5.8615; epsilon=.2849;

%name='/scratch1/taylorm/iso12_500A0001.7723'
%nx=500; delx_over_eta=2.740; epsilon=3.5208;

%name='/scratch1/taylorm/iso12_250A0022.000'
%nx=250; delx_over_eta=.80; epsilon=3.9;


%name='/ccs/scratch/taylorm/check256_0000.8000'

%name='/ccs/scratch/taylorm/dns/iso12_5120002.7000'

%name='/ccs/scratch/taylorm/dns/iso12/iso12_512'
%nx=512; delx_over_eta=2.74; epsilon=3.89;  teddy=1.0

name='/home2/skurien/helicity_data/helical_forced/hel256_hpi2/hel256_hpi2_'
nx=256; delx_over_eta=2.97; epsilon=2.72; teddy=1.24; % R_l=186


ext='.new.isostr';
%times=[1:.1:7.0];
times = [4.2:0.2:4.8];

%name='/ccs/scratch/taylorm/dns/sc1024A/sc1024A'
%nx=2048; delx_over_eta=2.98; epsilon=3.74; teddy=1.024;
%ext='.new.isostr';
%times=[1:.1:2.0];

%name = '/home2/skurien/helicity_data/isostr_1/check256_hq_';
%ext='.new.isostr';
[avg_eps, avg_heps, avg_delx_over_eta] = ensemble_avg_params(name,ext,times)
nx=256; delx_over_eta=avg_delx_over_eta; epsilon=avg_eps; h_epsilon=avg_heps;  

%teddy=1;
%times=[0:1:30];

ndir_use=0;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')

% this type of averging is expensive:
time_and_angle_ave=1;

k=0;


xx=(1:.5:(nx./2.5)) / nx;
xx_plot=(1:.5:(nx./2.5)) *delx_over_eta;   % units of r/eta

y45_ave=zeros([length(xx),73]);
if (time_and_angle_ave==1) 
   y45_iso_ave=zeros([length(xx),1]);
end

times_plot=[];
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10)];
  disp([fname,ext]);
  fid=fopen([fname,ext]);
  if (fid<0) ;
    disp('error openining file, skipping...');
  else
    fclose(fid);
    times_plot=[times_plot,t];
    k=k+1;

    if (time_and_angle_ave==1) 
      klaws=1;                          % compute 4/5 laws
      plot_posneg=0;
      check_isotropy=0;
      
      [y45,y415,y43,eps,h_eps]=compisoave(fname,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,1);
      
      
      mx45_iso_localeps(k)=max(y45);
      mx45_iso(k)=max(y45)*eps/epsilon;
      
      y45_iso_ave=y45_iso_ave+y45';
      
    end    


    [nx,ndelta,ndir,r_val,ke,eps_l,mu,D_ll,D_lll] ...
        = readisostr( [fname,ext] );
    
    eta_l = (mu^3 / eps_l)^.25;
    delx_over_eta_l=(1/nx)/eta_l;
    
    for dir=1:73;
      x=r_val(:,dir)/nx;                % box length
      y=-D_lll(:,dir)./(x*eps_l);
      
      y45 = spline(x,y,xx);
      
      mx45_localeps(k,dir)=max(y45);
      mx45(k,dir)=max(y45)*eps/epsilon;
      
      y45_ave(:,dir)=y45_ave(:,dir)+y45';
    end
  end

end

times=times_plot;
y45_ave=y45_ave/length(times);
y45_iso_ave=y45_iso_ave/length(times);

if (0)

% averging starting at t=0:
save k45data_t0 teddy times mx45_localeps mx45_iso_localeps ...
      y45_ave y45_iso_ave xx_plot

% averging starting at t=1:
save k45data_t1 teddy times mx45_localeps mx45_iso_localeps ...
      y45_ave y45_iso_ave xx_plot

% load data from t=0 (averaged wrong!) for isoave paper:
load k45data_t0 

% load data from t=1 (for time averaged) for isoave paper:
load k45data_t1 
end


figure(8); clf; hold on; 
for i=2:2
%   plot(times/teddy,mx45_localeps(1:length(times),i),'k:','LineWidth',2.0);
   plot(times/teddy,mx45_localeps(1:length(times),i),'g-','LineWidth',2.0);
end
%plot(times/teddy,mx45_iso_localeps,'k-','LineWidth',2.0);
plot(times/teddy,mx45_iso_localeps,'b-','LineWidth',2.0);
ax=axis;
axis( [ax(1),ax(2),.5,1.0] );
plot(times,(4/5)*times./times,'k');
hold off;
ylabel(' < (u(x+r)-u(x))^3 > / (\epsilon r)','FontSize',16);
xlabel('time','FontSize',16)
print -dpsc k45time.ps


figure(9); clf
scale = 1; %scale = 4/5 if need to normalize out the 4/5th
for i=[1:1:73]
  %semilogx(xx_plot,y45_ave(:,i),'k:','LineWidth',1.0); hold on
  semilogx(xx_plot,y45_ave(:,i)/scale,'g-','LineWidth',1.0); hold on
end
%semilogx(xx_plot,y45_iso_ave/scale,'k','LineWidth',1.0); hold on
semilogx(xx_plot,y45_iso_ave/scale,'b','LineWidth',1.0); hold on
axis([1 1000 0 1.0])
x=1:1000; plot(x,(4/5)*x./x,'k');
hold off;
%title('D_{lll} / r\epsilon   (4/5 law) ');
ylabel('< (u(x+r)-u(x))^3 > / (\epsilon r)','FontSize',16);
xlabel('r/\eta','FontSize',16);

print -dpsc k45mean.ps



figure(10); clf
offset=y45_iso_ave;
stdr=0*offset;
scale = 4/5; % scale = 4/5 if need to factor out the 4/5 
for i=[1:1:73]
  %semilogx(xx_plot,(y45_ave(:,i)-offset),'k:','LineWidth',1.0); hold on
  semilogx(xx_plot,(y45_ave(:,i)-offset)/scale,'m-','LineWidth',1.0); hold on
  stdr=stdr+(y45_ave(:,i)-offset).^2;
end
stdr=sqrt(stdr/15)./offset;
%semilogx(xx_plot,y45_iso_ave-offset,'k','LineWidth',1.0); hold on
semilogx(xx_plot,(y45_iso_ave-offset)/scale,'b','LineWidth',1.0); hold on
title('Measure of anisotropy - 4/5 law');
%ylabel('< (u(x+r)-u(x))^3 > / (\epsilon r)','FontSize',16);
xlabel('r/\eta','FontSize',16);






starttime=1;
ln=find(times>=starttime);
ln=ln(1);
lnmax=length(times);

dir_use=2;
[mean(mx45(ln:lnmax,dir_use)),mean(mx45_iso(ln:lnmax))]
[std(mx45(ln:lnmax,dir_use)),std(mx45_iso(ln:lnmax)) ]


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

name='/ccs/scratch/taylorm/dns/iso12_512'
nx=512; delx_over_eta=2.74; epsilon=3.89;  teddy=1.0
%ext='.isostr001';



ndir_use=0;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')

% this type of averging is expensive:
time_and_angle_ave=1;

k=0
times=[1:.1:3.7];



xx=(1:.5:(nx./2.5)) / nx;
xx_plot=(1:.5:(nx./2.5)) *delx_over_eta;   % units of r/eta

y45_ave=zeros([length(xx),15]);
if (time_and_angle_ave==1) 
   y45_iso_ave=zeros([length(xx),1]);
end


for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10)];
  disp([fname,ext])
  k=k+1;

  if (time_and_angle_ave==1) 
    klaws=1;                            % compute 4/5 laws
    plot_posneg=0;
    check_isotropy=0;
    
    [y45,y415,y43,eps]=compisoave(fname,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,0);


    mx45_iso_localeps(k)=max(y45);
    mx45_iso(k)=max(y45)*eps/epsilon;

    y45_iso_ave=y45_iso_ave+y45';

  end    
  
  [nx,ndelta,ndir,r_val,ke,eps_l,mu,D_ll,D_lll] ...
      = readisostr( [fname,ext] );
  
  eta_l = (mu^3 / eps_l)^.25;
  delx_over_eta_l=(1/nx)/eta_l;
  
    for dir=1:15;
      x=r_val(:,dir)/nx;         % box length
      y=-D_lll(:,dir)./(x*eps_l);
      
      y45 = spline(x,y,xx);
      
      mx45_localeps(k,dir)=max(y45);
      mx45(k,dir)=max(y45)*eps/epsilon;
      
      y45_ave(:,dir)=y45_ave(:,dir)+y45';
    end
end

y45_ave=y45_ave/length(times);
y45_iso_ave=y45_iso_ave/length(times);



figure(4); clf; hold on; 
for i=2:2
   plot(times/teddy,mx45_localeps(:,i));
end
plot(times/teddy,mx45_iso_localeps,'g')
ax=axis
axis( [ax(1),ax(2),.5,1.0] );
hold off;

figure(5); clf
for i=1:3
  semilogx(xx_plot,y45_ave(:,i)); hold on
end
  semilogx(xx_plot,y45_iso_ave,'k'); hold on
axis([1 1000 0 1.0])
hold off;


starttime=1;
ln=find(times>=starttime);
ln=ln(1);
lnmax=length(times);

dir_use=2;
[mean(mx45(ln:lnmax,dir_use)),mean(mx45_iso(ln:lnmax))]
[std(mx45(ln:lnmax,dir_use)),std(mx45_iso(ln:lnmax)) ]








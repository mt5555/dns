%   
%   
% pv-velocity 2/3ths law: time and angle averaged 
%   
%   
mu=0;
ke=0;
nx=1;
delx_over_eta=1;
ext='.bisostr';
eta = 1/(nx*delx_over_eta);

%name = '~/projects/pv/data_analysis/lowforc/qg';
%name = '~/projects/pv/data_analysis/lowforc/low4/qg64/qg';
%times=[200:2:492];
%nx=64;
%name = '~/data/qg';
%name = '~/projects/pv/data_analysis/lowforc/low4/qg256/qg';
%pname = 'qg';

%times=[30:2:102];

%name = '~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/qg64';
%pname='qg64';
%times=[1700:5:1795];
%nx=64;

name = '~/projects/pv/data_analysis/lowforc/low4/qg256/bous2000/'
%pname='qg256hyper'
%times=[0:.1:3.3];
%nx=256;

name = '~/projects/INCITE_runs/Intrepid/bous_NSvisc/'
%name = '~/INCITE_runs/Intrepid/bous_NSvisc/'
%pname='n256_Ro1Fr0.01_'
%times=[1:.1:1.9];
%nx = 256;
pname='n512_Ro1Fr0.01_'
times=[.1:.1:4.8];
nx=512;


%check this subroutine, for now set averages to 1
%[avg_eps, avg_heps, avg_delx_over_eta] = ensemble_avg_params(name,ext,times)



delx_over_eta=1; epsilon=1; h_epsilon=1;Q_eps_l=1;  
teddy=1;

ndir_use=73;
%ndir_use=49;  disp('USING ONLY 49 DIRECTIONS')
%

if (ndir_use>0) ndir=ndir_use; end;

% $$$ if (ndir==3)
% $$$   w=ones([1,3])/3;
% $$$ else
% $$$   equalw=0;
% $$$   if (equalw) 
% $$$     % put this in to use equally weighted:
% $$$     w=ones([1,ndir])/ndir;
% $$$     disp(sprintf('Using equal weights for spherical integration'))
% $$$   else
% $$$     % get the weights:
% $$$     wname=sprintf('../src/voronoi/isoave.weights%i',ndir);
% $$$     disp(sprintf('Reading Voronio weights from file:  %s',wname))
% $$$     w=textread(wname,'%f');
% $$$     % take every other weight
% $$$     w=2*w(1:2:length(w));
% $$$   end
% $$$ end
% $$$ if (abs(1-sum(w))>1e-7) 
% $$$   disp('error: weights do not sum to 1')
% $$$   return;
% $$$ end



% this type of averging is expensive:
time_and_angle_ave=1;

k=0;


xx=(1:.5:(nx./2.5)) / nx;
xx_plot=(1:.5:(nx./2.5)) *delx_over_eta;   % units of r/eta

y23_ave=zeros([length(xx),73]);
if (time_and_angle_ave==1) 
   y23_iso_ave=zeros([length(xx),1]);
end

times_plot=[];
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,pname,tstr(2:10)];
  disp([fname,ext])
ppname = [pname,tstr(2:10),ext]
  fid=fopen([fname,ext])
  if (fid<0) ;
    disp('error opening file, skipping...');
  else
    fclose(fid);
    times_plot=[times_plot,t];
    k=k+1;

    if (time_and_angle_ave==1) 
      klaws=5;                          % compute 2/3 laws
      plot_posneg=0;
      check_isotropy=0;
      [y45,y415,y43,epsl,h_eps,y215,y23]=compisoave(fname,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,0);
      
      y23_iso_ave=y23_iso_ave+y23';  %accumulate the angle averaged y23 from each frame
    figure(10); plot(k*.1,max(y23),'*');hold on;
    end    



    [nx,ndelta,ndir,r_val,ke,eps_l,mu,...
     D_ll,D_lll,D1_tt,D2_tt,D1_ltt,D2_ltt,...
     SP_lll,SN_lll,SP1_ltt,SP2_ltt,SN1_ltt,SN2_ltt,H_ltt,H_tt,D_lltt,Dl,Dt,...
     h_eps_l,Q_eps_l] ...
        = readisostr( [fname,ext] );
    
    eta_l = (mu^3 / eps_l)^.25;
    delx_over_eta_l=(1/nx)/eta_l;

    Q_eps_l   
    for dir=1:ndir
      x=r_val(:,dir)/nx;               % box length
      y=Dl(:,dir,1)./(Q_eps_l*x);               % for forced data

%      y23 = w(dir)*spline(x,y,xx);
      y23 = spline(x,y,xx);
      y23_ave(:,dir)=y23_ave(:,dir)+y23';  % accumulate the y23 for each direction, from each frame.


    end
  end

end

times=times_plot;
ll=length(times);
y23_ave=y23_ave/ll;
y23_iso_ave=y23_iso_ave/ll;

scale = 1;  % scale = 2/3 if normalizing by prefactor as well


figure(6);
scale = 1; % scale = 2/3 to scale out
for i=[1:ndir]
semilogx(xx_plot,(y23_ave(:,i))/scale,'k:','LineWidth',1.0); hold on
%pause
end
semilogx(xx_plot,y23_iso_ave/scale,'k','LineWidth',2.0); hold on
%axis([1 1000 -0.05 0.15])
x=1:1000; plot(x,(2/3)*x./x/scale,'k');
hold on;
title('  (2/3 law) Ensemble averaged');
ylabel(ppname);
xlabel('r/\eta','FontSize',16);


figure(8);
scale = 1; % scale = 2/3 to scale out
loglog(xx_plot,abs(y23_iso_ave.*xx_plot')/scale,'k','LineWidth',2.0); hold on
%axis([1 1000 -0.05 0.15])
hold on;
title('  (2/3 law) Ensemble averaged');
ylabel(ppname);
xlabel('r/\eta','FontSize',16);



%print -dpsc k215mean.ps


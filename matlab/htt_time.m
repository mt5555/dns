%
%
% Second-order helical structure function;  time and angle averaged 
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

name = '/home2/skurien/helicity_data/helical_forced/hel256_hpi2/hel256_hpi2_';
pname = 'hel256\_hpi2\_';
ext='.new.isostr';
times=[4.2:0.2:15.8];
nx = 256;

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

htt_ave=zeros([length(xx),73]);
if (time_and_angle_ave==1) 
   ytt_iso_ave=zeros([length(xx),1]);
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
      
      [y45,y415,y43,eps,h_eps,y215,ytt]=compisoave(fname,ext,xx,ndir_use,klaws,plot_posneg,check_isotropy,0);
      
      
      ytt_iso_ave=ytt_iso_ave+ytt';
      
    end    


    [nx,ndelta,ndir,r_val,ke,eps_l,mu,tmp,tmp,tmp,tmp,tmp,tmp,tmp,...
    tmp,tmp,tmp,tmp,tmp,H_ltt,H_tt,tmp,tmp,tmp,h_eps_l] ...
        = readisostr( [fname,ext] );
    
    eta_l = (mu^3 / eps_l)^.25;
    delx_over_eta_l=(1/nx)/eta_l;

    
    for dir=1:73;
      x=r_val(:,dir)/nx;                % box length
      y=H_tt(:,dir);
      
      htt = spline(x,y,xx);
      
      htt_ave(:,dir)=htt_ave(:,dir)+htt';
    end
  end

end

times=times_plot;
htt_ave=htt_ave/length(times);
ytt_iso_ave=ytt_iso_ave/length(times);


figure(6); clf
for i=[1:1:73]
  semilogx(xx_plot,abs(htt_ave(:,i))./(xx_plot'.^(2/3)),'r-','LineWidth',1.0); hold on
end
semilogx(xx_plot,abs(ytt_iso_ave)./(xx_plot'.^(2/3)),'b','LineWidth',2.0); hold on
max(y215_iso_ave)
%axis([1 1000 -0.05 0.15])
title('H_{tt} /r^{4/3}  ensemble averaged');
ylabel(ppname);
xlabel('r/\eta','FontSize',16);
hold off;
print -dpsc httmean.ps














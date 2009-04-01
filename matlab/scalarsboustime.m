%########################################################################
%#  read in Boussenesque scalars   *.scalars-bous
%########################################################################
%
%
clear all;


% defaults:
nx=0;
fcor=0;
f_k=0;

fid2=-1;


%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg64/sto_high_4/hyper_nu/bous100/qg64hyper_all.scalars-bous','r');
%fid=endianopen(['~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg64/' ...
%                'sto_high_4/hyper_nu/bous500/qg64hyper_all.scalars-bous'],'r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg64/sto_high_4/hyper_nu/bous1000/qg64hyper_all.scalars-bous','r');

%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg64/iso12w/qg64_iso12w_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/bous100/qg64_sto16_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/bous200/qg64_200all.scalars-bous','r');

%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg128/sphere_visc/qg128sph_all.scalars-bous','r');
%fid=endianopen(['~/research.old/projects/pv/data_analysis/' ...
%                'lowforc/low4/qg/qg128/slab_visc/qg128slab_all.scalars-bous'],'r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg128/sphere_dealias/qg128_fftsph_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg128/no_hyper/qg128_nhyper_all.scalars-bous','r');


%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous500/qg256hyper_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous1000/qg256hyper_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous2000/qg256hyper_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg256/fcor2000_bous1000/qg256hyper_all.scalars-bous','r');

%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous1000/qg512_b1000_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous1000/new_data/qg512_b1000_newdata_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous1000/hyper4/qg512hyper4_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous2000/correct_hyper/qg512_b2000_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous2000/qg512hyper_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/fcor2000_bous20/n512_f2000b20_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/fcor2000_bous200/n512_f2000b200_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/fcor20_bous2000/n512_f20b2000_all.scalars-bous','r');
%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg512/fcor200_bous2000/n512_f200b2000_all.scalars-bous','r');

fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg640/qg640_b3000_all.scalars-bous','r');


%fid = endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/n256/n256_f2000n5_all.scalars-bous','r');
%fid = endianopen(['~/research.old/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/' ...
%                  'n256/n256_f1000n5_all.scalars-bous'],'r');
%fid = endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/Ro1Fr0/n256/n256_f5n2000_all.scalars-bous','r');
%fid = endianopen(['~/research.old/projects/pv/data_analysis/lowforc/low4/Ro1Fr0/' ...
%                  'n256/n256_f5n1000_all.scalars-bous'],'r');

%fid = endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/n256/n256high_f2000n5_all.scalars-bous','r');
%fid = endianopen(['~/research.old/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/' ...
%                  'n256/n256high_f1000n5_all.scalars-bous'],'r');

fid = endianopen('~/INCITE_runs/SW02_tests/bous128_Ro21Fr0.21_all.scalars-bous','r')

fid=endianopen('~/projects/bous640runs/qg640_b3000_all.scalars-bous','r');

fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/mauler/f27b136/bous200Lz0.2_all.scalars-bous','r')

nscalars=0;
ints=[];
time=[];
while (1) 
  [ni,count]=fread(fid,1,'float64');
  if (count~=1) break;   end
  nints=ni;
  data=fread(fid,1,'float64');
  time=[time,data];
  data=fread(fid,ni,'float64');
  % might need to take the transpose of 'data' here:
  ints=[ints,data];
  nscalars = nscalars+1;
end;  


disp(sprintf('number of integrals=%i  number of times read=%i',nints,nscalars))
fclose(fid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the scalars computed every time step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

ke=ints(1,:);
pe=ints(2,:);
ke_diss = ints(3,:);
pe_diss = ints(4,:);
pv=ints(5,:);
potens=ints(7,:);
potens_diss=ints(8,:);
potens_qg=ints(10,:);
potens_ro0fr1 = ints(12,:);
potens_ro1fr0 = ints(14,:);


tote = ke + pe;

time_2=[];


figure(1)
clf
hold on
plot(time,ke,'r-')
plot(time,pe,'g-')
plot(time,tote,'k-')
title('KE: red  PE:  green  Etot: black');
hold off

figure(2)
clf
hold on
plot(time,pv,'g-')
title('PV');
hold off

figure(3)
subplot(2,1,2);
%clf
plot(time,potens,'b','Linewidth',2);hold on;
plot(time,potens_qg,'ro','Markersize',6);
%plot(time,potens_ro0fr1,'r.-');
%plot(time,potens_ro1fr0,'k-');
%title('Potential enstrophy: (blue: Total); (red o QG); (red.-: Ro->0,Fr1); (black: Ro1, Fr->0)');
legend('$Q$','$Q_{qg}$','$Q_{q\sim f \partial_z \rho}$', '$Q_{q\sim N \omega_3}$');

figure(4)
clf
hold on
dt = diff(time);
dqdt = diff(potens)./dt;
plot(time(1:length(time)-1),dqdt,'b');
plot(time(1:length(time)), -potens_diss,'k');
plot(time(1:length(time)-1), (dqdt + potens_diss(1:length(time)-1)),'m');
mean_forc = mean((dqdt + potens_diss(1:length(time)-1)))
legend('\Delta Q/\Delta t', 'potens\_diss','difference');
grid on;
title('potential enstrophy dissipation rate dQ/dt check')

print -depsc bous-scalars.ps


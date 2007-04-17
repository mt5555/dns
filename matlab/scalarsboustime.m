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


%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg64/sto_high_4/hyper_nu/bous100/qg64hyper_all.scalars-bous','r');
%fid=endianopen(['~/projects/pv/data_analysis/lowforc/low4/qg/qg64/' ...
%                'sto_high_4/hyper_nu/bous500/qg64hyper_all.scalars-bous'],'r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg64/sto_high_4/hyper_nu/bous1000/qg64hyper_all.scalars-bous','r');

%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/iso12w/qg64_iso12w_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/bous100/qg64_sto16_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/bous200/qg64_200all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous500/qg256hyper_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous1000/qg256hyper_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg256/bous2000/qg256hyper_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg256/fcor2000_bous1000/qg256hyper_all.scalars-bous','r');


fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg/qg512/bous2000/qg512hyper_all.scalars-bous','r');

%fid = endianopen('~/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/n256/n256_f2000n5_all.scalars-bous','r');
%fid = endianopen(['~/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/' ...
%                  'n256/n256_f1000n5_all.scalars-bous'],'r');
%fid = endianopen('~/projects/pv/data_analysis/lowforc/low4/Ro1Fr0/n256/n256_f5n2000_all.scalars-bous','r');
%fid = endianopen(['~/projects/pv/data_analysis/lowforc/low4/Ro1Fr0/' ...
%                  'n256/n256_f5n1000_all.scalars-bous'],'r');

%fid = endianopen('~/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/n256/n256high_f2000n5_all.scalars-bous','r');
%fid = endianopen(['~/projects/pv/data_analysis/lowforc/low4/Ro0Fr1/' ...
%                  'n256/n256high_f1000n5_all.scalars-bous'],'r');


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
%clf
plot(time,potens,'b');hold on;
plot(time,potens_qg,'r--');
plot(time,potens_ro0fr1,'r.-');
plot(time,potens_ro1fr0,'k-');
title('Potential enstrophy: (blue: Total); (red--: QG); (red.-: Ro->0,Fr1); (black: Ro1, Fr->0)');


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


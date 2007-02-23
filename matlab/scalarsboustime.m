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


%fid=endianopen('~/Desktop/qg64_sto16_all.scalars-bous','r');

%fid=endianopen('~/projects/pv/data_analysis/lowforc/low3/qg64/qg64_low3_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low3/qg64/qg0194.0000.scalars-bous','r');
%fid=endianopen(['~/projects/pv/data_analysis/lowforc/low4/qg64/qg64_low4_all.scalars-bous','r');
fid=endianopen('~/tmp/new/qg64_all.scalars-bous','r');

%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/iso12/qg64_low4_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/noforc/qg64all_noforc.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg256/qg256_all8.0.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg256/qg256_all100.0.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/iso23w/reg_nu/qg64_iso23w_all.scalars-bous','r');

%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/iso23w/hyper_nu/qg64hyper_all15.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_4/hyper_nu/qg64hyper_all10.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/iso12w/qg64_iso12w_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/bous100/qg64_sto16_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/sto_high_16/bous200/qg64_200all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg256/bous1000/hyper_nu2.5/qg256hyper_all.scalars-bous','r');
fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg256/bous1000/hyper_nu15/qg256hyper_all13.scalars-bous','r');

%fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/iso23w/hyper_nu/qg64hyper_all15.scalars-bous','r');

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
pv=ints(5,:);
potens=ints(7,:);
potens_diss=ints(8,:);

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
clf
hold on
plot(time,potens,'g');hold on;
title('Potential enstrophy');
hold off

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


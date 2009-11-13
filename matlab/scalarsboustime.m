%########################################################################
%#  read in Boussenesque scalars   *.scalars-bous
%########################################################################
%
%
clear all;


% defaults:
nx=0;
fcor=0;
f_k=4;

fid2=-1;



%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg640/qg640_b3000_all.scalars-bous','r');

%fid = endianopen('~/INCITE_runs/SW02_tests/bous128_Ro21Fr0.21_all.scalars-bous','r')


fid=endianopen('~/projects/INCITE_runs/Intrepid/qg/n640_bous3000_all.scalars-bous','r');
fid=endianopen('~/projects/INCITE_runs/Intrepid/Ro1Fr0/n640_fcor14bous3000_all.scalars-bous','r');
%fid=endianopen('~/projects/INCITE_runs/Intrepid/Ro0Fr1/n640_fcor3000bous14_all.scalars-bous','r');

%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/n256_Ro1Fr0.01_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/n512_Ro1Fr0.01_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/n1024_Ro1Fr0.01_all.scalars-bous','r')

fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_Ro1Fr0.002_all.scalars-bous','r')

%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_bous/shift_force/n1600_d0.2_Ro0.05_all.scalars-bous','r')

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
subplot(1,1,1);
%clf
plot(time,potens,'b','Linewidth',2);hold on;
%plot(time,potens_qg,'r--','Markersize',6);
%plot(time,potens_ro0fr1,'ro','Markersize',6);
plot(time,potens_ro1fr0,'k*','Markersize',6);
title('Potential enstrophy: (blue: Total); (red -- QG); (red o: Ro->0,Fr1); (black *: Ro1, Fr->0)');
%legend('$Q$','$Q_{qg}$','$Q_{q\sim f \partial_z \rho}$', '$Q_{q\sim N \omega_3}$');
legend('$Q$','$|N \omega_3|^2$')

%check that dissipation matches dQ/dt for decaying case
if(0)
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
end

print -depsc bous-scalars.ps


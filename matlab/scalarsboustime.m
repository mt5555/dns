%########################################################################
%#  read in Boussenesque scalars   *.scalars-bous
%########################################################################
%
%
clear all;


% defaults:
nx=0;
fcor=0;
bous=0;
f_k=4;

fid2=-1;



%fid=endianopen('~/research.old/projects/pv/data_analysis/lowforc/low4/qg/qg640/qg640_b3000_all.scalars-bous','r');

%fid = endianopen('~/INCITE_runs/SW02_tests/bous128_Ro21Fr0.21_all.scalars-bous','r')


%fid=endianopen('~/projects/INCITE_runs/Intrepid/qg/n640_bous3000_all.scalars-bous','r');
%fid=endianopen('~/projects/INCITE_runs/Intrepid/Ro1Fr0/n640_fcor14bous3000_all.scalars-bous','r');
%fid=endianopen('~/projects/INCITE_runs/Intrepid/Ro0Fr1/n640_fcor3000bous14_all.scalars-bous','r');
%nx=640;ny=640;nz=640;

%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/n256_Ro1Fr0.01_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/n512_Ro1Fr0.01_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.01/n1024_Ro1Fr0.01_all.scalars-bous','r')

%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_nu.2e-4/n1024_Ro1Fr0.002_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_nu.1e-4/n1024_Ro1Fr0.002_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_nu.7e-5/n1024_Ro1Fr0.002_all.scalars-bous','r')
fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.002/n1024_nu.5e-5/n1024_Ro1Fr0.002_nu.5e-5_all.scalars-bous','r')

%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.001/n1024_nu.5e-5/n1024_Ro1Fr0.001_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.001/n1024_nu.5e-5/n1024_Ro1Fr0.001_0003.6000.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro1Fr0.001/n1024_nu.5e-5/n1024_Ro1Fr0.001_0003.7000.scalars-bous','r')
%nx=1024;ny=1024;nz=1024;

%fid=endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.002Fr1/n1024_nu.7e-5/n1024_Ro0.002Fr1_all.scalars-bous','r')

%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr1/n1024_nu.5e-5/n1024_Ro0.001Fr1_all.scalars-bous','r')
%nx=1024;ny=1024;nz=1024;


%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.002Fr0.002/n1024_Ro0.002Fr0.002_all.scalars-bous','r')
%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.002Fr0.002/n1024_Ro0.002Fr0.002_old_all.scalars-bous','r')
%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.002Fr0.002/n1024_Ro0.002Fr0.002_new_all.scalars-bous','r')
%nx=1024;ny=1024;nz=1024;

%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr0.001/n1024_Ro0.001Fr0.001_all_old.scalars-bous','r')
%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr0.001/n1024_Ro0.001Fr0.001_all.scalars-bous','r')
%fid =endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr0.001/n1024_Ro0.001Fr0.001_0000.6000.scalars-bous','r')
%fid=endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr0.001/n1024_Ro0.001Fr0.001_0.5delt.3all.scalars-bous','r')
%fid=endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr0.001/n1024_Ro0.001Fr0.001_0.5delt.2all.scalars-bous','r')
%fid=endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/Ro0.001Fr0.001/n1024_Ro0.001Fr0.001_0.5delt.1_all.scalars-bous','r')
%nx=1024;ny=1024;nz=1024;

%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_bous/shift_force/n1600_d0.2_Ro0.05_all.scalars-bous','r')

%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/LOWRES/sto_high_t4/n512_d0.25_Ro0.05_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/LOWRES/sto_high_t4/n512_d0.25_Ro0.05_-newall.scalars-bous','r')
%nx=512;ny=512;nz=128;

%fid =endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n1024_d0.25_Ro0.05_nodamp/n1024_d0.25_Ro0.05_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n1024_d0.25_Ro0.05_nodamp/n1024_d0.25_Ro0.05_-newall.scalars-bous','r')
fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n1024_d0.25_Ro0.002_nodamp/n1024_d0.25_Ro0.002_all.scalars-bous','r')
%nx=1024;ny=1024;nz=256;

%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.05_nodamp/n2048_d0.25_Ro0.05_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.05_nodamp/n2048_d0.25_Ro0.05_-newall.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.01_nodamp/n2048_d0.25_Ro0.01_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.01_nodamp/n2048_d0.25_Ro0.01_new2_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.005_nodamp/n2048_d0.25_Ro0.005_all.scalars-bous','r')
fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n2048_d0.25_Ro0.002_nodamp/n2048_d0.25_Ro0.002_all.scalars-bous','r')
%nx=2048;ny=2048;nz=512;


%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n512_d1.0_Ro0.05_nodamp/n512_d1.0_Ro0.05_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n512_d1.0_Ro0.05Fr1_nodamp/n512_d1.0_Ro0.05Fr1_all.scalars-bous','r')
%fid = endianopen('~/projects/INCITE_runs/Intrepid/lowaspect_bous/n512_d1.0_Ro1Fr0.05_nodamp/n512_d1.0_Ro1Fr0.05_all.scalars-bous','r')
%nx=512;ny=512;nz=512


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
if length(ints(:,1)) >= 16
potens_Lz = ints(16,:);
else potens_Lz = 0;
end
tote = ke + pe;

time_2=[];


figure(1)
%clf
hold on
plot(time,ke,'rx')
plot(time,pe,'gx')
plot(time,tote,'kx')
title('KE: red  PE:  green  Etot: black');
%hold off

figure(2)
%clf
hold on
plot(time,pv,'gx')
title('PV');
%hold off

figure(5)
size=1; % nx*ny*nz;
subplot(1,1,1);
set(gca,'fontsize',16);
%clf
tscale=1;
%if (max(fcor,bous) ~=0) then
%    tscale = max(fcor,bous);
plot(time,potens/size,'b-','Linewidth',2);hold on;
plot(time,potens_qg/size,'r--','Markersize',6);
plot(time,potens_ro0fr1/size,'ro','Markersize',6);
plot(time,potens_ro1fr0/size,'k*','Markersize',6);
if (length(potens_Lz) > 1)
plot(time,potens_Lz/size,'rx','Markersize',10);
legend('Q', 'Q_{qg}','0.5 |f \partial_z \rho |^2', '0.5 |N \omega_3|^2','Q_{Lz}')
%title('Potential enstrophy: (blue: Total); (red -- QG); (red o: Ro->0,Fr1); (black *: Ro1, Fr->0)');
%legend('$Q$','$Q_{qg}$','$Q_{q\sim f \partial_z \rho}$', '$Q_{q\sim N \omega_3}$');
else
legend('Q', 'Q_{qg}','0.5 |f \partial_z \rho |', '0.5 |N \omega_3|^2')
end
xlabel('time t')

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


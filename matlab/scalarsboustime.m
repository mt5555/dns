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


%fid=endianopen('~/projects/pv/data_analysis/lowforc/low3/qg64/qg64_low3_all.scalars-bous','r');
%fid=endianopen('~/projects/pv/data_analysis/lowforc/low3/qg64/qg0194.0000.scalars-bous','r');
fid=endianopen('~/projects/pv/data_analysis/lowforc/low4/qg64/qg64_low4_all.scalars-bous','r');
f_k= 24;



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
pv=ints(5,:)
potens=ints(7,:)
potens_diss=ints(8,:)

tote = ke + pe

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
plot(time(1:length(time)-1),diff(potens)./diff(time),'b');
title('potential enstrophy dissipation rate dQ/dt')

print -depsc bous-scalars.ps


%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
clear all;
apply_fix=0;

% defaults:
nx=0;
fcor=0;
f_k=0;

fid2=-1;


fid=endianopen('/home/wingate/Projects/KH/Boussinesq/n21/all.scalars-bous','r');
f_k= 24;


nscalars=0
nscalars_e=0

while (1) 
  [ni]=fread(fid,1,'float64')
  nints=ni
  time=fread(fid,1,'float64')
  data1=fread(fid,[nints],'float64')
  ints=[data1]
end;  
  nints=ni;
  
  nints;

disp(sprintf('nints=%i  total scalars read=%i',nints,nscalars))
fclose(fid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the scalars computed every time step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
l=size(ints); 
l=l(2); 

ke=ints(1);
pe=ints(2);
pv=ints(5);
potens=ints(7);


tote = ke + pe;

time_2=[];



figure(1)
clf
hold on
plot(time,ke,'r')
plot(time,pe,'g')
plot(time,tote,'k')
title('KE: red  PE:  green  Etot: black');
hold off

figure(2)
clf
hold on
plot(time,pv,'g')
title('PV');
hold off

figure(3)
clf
hold on
plot(time,potens,'g')
title('Potential enstrophy');
hold off


print -depsc bous-scalars.ps


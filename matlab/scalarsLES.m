%
%########################################################################
%#  plot of DNS structure function output file
%########################################################################
%
%
clear all;
apply_fix=0;

range=0:50;
%range = 1
%range=0:.5:2;

i=1;
%fidarray(i)=fopen('../src/impulse/kh0.scalars','r','l');  i=i+1;
fidarray(i)=fopen('../src/impulse/kh10.scalars','r','l'); i=i+1;
fidarray(i)=fopen('../src/impulse/kh17.scalars','r','l'); i=i+1;
fidarray(i)=fopen('../src/impulse/kh18.scalars','r','l'); i=i+1;
fidarray(i)=fopen('../src/impulse/kh19.scalars','r','l'); i=i+1;
fidarray(i)=fopen('../src/impulse/kh20.scalars','r','l'); 

figure(5)
clf
figure(6)
clf

for index=1:length(fidarray);

fid=fidarray(index);


nscalars=0;
nscalars_e=0;
while (1) 
  [ni,count]=fread(fid,1,'float64');
  if (count~=1) break; end;
  nints=ni;
  ns=fread(fid,1,'float64');
  mu=fread(fid,1,'float64');
  alpha=fread(fid,1,'float64');
  
  [nints,ns,nscalars]
  data1=fread(fid,[nints,ns],'float64');
  data2=fread(fid,[nints,ns],'float64');

  if (nscalars==0) 
    ints=data1;
    maxs=data2;
  else
    ints=[ints,data1];
    maxs=[maxs,data2];
  end
  nscalars=nscalars+ns;

  % now read the "expensive" integrals, which are not computed every time step
  ns_e = fread(fid,1,'float64');
  time = fread(fid,1,'float64');
  data1 = fread(fid,[ns_e,1],'float64');
  data1=[time;data1];
  if (nscalars_e==0) 
    ints_e= data1;
  else
    ints_e= [ints_e,data1];
  end
  nscalars_e=nscalars_e+1;
end

disp(sprintf('nints=%i  total scalars read=%i',nints,nscalars))
fclose(fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the scalars computed every time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=size(ints);
l=l(2);

ke_diss_d=-mu*ints(2,:);
ke_diss_f=ints(3,:);        % < u,F >   F=f, except for alpha model F=f'
vor_z=ints(4,:);
hel=ints(5,:);
ke=ints(6,:);   %  ke 
% ints(7,:)   enstrophy
% ints(8,:)   < u,div(tau)' >  (alpha model only)
% ints(9,:)   < u,f>           (alpha model only)
% ints(1,:)  < u_xx,u_xx> >   (used for E_alpha dissapation term)

maxU=maxs(1,:);  % at time_after
maxV=maxs(2,:);  % at time_after
maxW=maxs(3,:);  % at time_after
%  maxs(4,:)     % max used for CFL, at time_after
maxvor=maxs(5,:);
time_after=maxs(6,:);
time=maxs(7,:);

Ea = ints(6,:) + .5*alpha^2 *ints(2,:); % at time

time_2 = .5*(time(2:l)+time(1:l-1));
ke_diss_tot=(ke(2:l)-ke(1:l-1))./(time(2:l)-time(1:l-1));
Ea_diss_tot=(Ea(2:l)-Ea(1:l-1))./(time(2:l)-time(1:l-1));



disp(sprintf('max vor_z = %e',max(vor_z)));

figure(5)
hold on
plot(time,ke-ke(1))
axis([0 .75 -.005 .025])
title('KE');

figure(6)
hold on
plot(time,Ea-Ea(1))
axis([0 .75 -.005 .025])
title('E_\alpha');


end

hold off

figure(5)
print -depsc scalars.ps
figure(6)
print -depsc scalars2.ps









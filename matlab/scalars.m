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
fid2=-1;

%fid=fopen('iso12_256_200.scalars','r','b'); 
%fid=fopen('../src/impulse/kh230000.0000.scalars','r','l'); 
%fid=fopen('../src/kh/khN.scalars','r','l'); 
%fid=endianopen('/ccs/scratch/taylorm/dns/iso12_512.scalars','r'); 
%fid2=endianopen('/ccs/scratch/taylorm/dns/iso12_512b.scalars','r'); 
fid=fopen('/ccs/scratch/taylorm/decay/decay2048.scalars','r','l'); 
%fid=fopen('../src/sht/rung0000.0000.scalars','r','l'); 




nscalars=0;
nscalars_e=0;
while (1) 
  [ni,count]=fread(fid,1,'float64');
  if (count~=1) 
     if (fid2<0) 
        break; 
     end
     fclose(fid);
     fid=fid2;
     fid2=-1;
     [ni,count]=fread(fid,1,'float64');
     if (count~=1) 
       break;
     end
  end;
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
end

disp(sprintf('nints=%i  total scalars read=%i',nints,nscalars))
fclose(fid);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the scalars computed every time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=size(ints);
l=l(2);

ke_diss_d=ints(10,:);
ke_diss_f=ints(3,:);        % < u,F >   F=f, except for alpha model F=f'
vor_z=ints(4,:);
hel=ints(5,:);
ke=ints(6,:);   %  ke 
ens = ints(7,:);   %  enstrophy
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
clf
hold on
plot(time,ke,'r')
%plot(time_2,ke_diss_tot,'r')
plot(time,ke_diss_f+ke_diss_d,'k')
%plot(time,ke_diss_f,'b')
plot(time,ke_diss_d,'b')
plot(time,hel,'g')
title('KE: red,    d(KE)/dt: blue,    hel: green');
hold off
xlabel('time')
print -dpsc ke.ps
print -djpeg -r72 ke.jpg

figure(6)
plot(time,sqrt(ens),'r')
title('average vorticity magnitude')
xlabel('time')
print -djpeg -r72 vor.jpg


lambda=sqrt(  5*(2*ints(6,:))./ints(2,:)  );
R_l = lambda.*sqrt(2*ints(6,:)/3)/mu;

figure(2); subplot(1,1,1)
plot(time,R_l)
title('R_\lambda')
xlabel('time')
print -djpeg -r72 rl.jpg


% Kolm. micro scale
eta = (mu^3 ./ abs(ke_diss_d)).^(.25);

figure(3); subplot(1,1,1)
plot(time,eta* 1024*2*pi*2*sqrt(2)/3 )
title('k_{nmax} \eta')
xlabel('time')
print -djpeg -r72 kmaxeta.jpg



% averge eta to a number
eta = eta(length(eta)/2:length(eta));
eta = sum(eta)/length(eta);
disp(sprintf('eta (average over last half of data) = %f ',eta));

% averge R_l to a number
R_l = R_l(length(R_l)/2:length(R_l));
R_l = sum(R_l)/length(R_l);
disp(sprintf('R_l (average over last half of data) = %f ',R_l));

disp(sprintf('1/250 in units of eta:  %f',(1/250)/eta));
disp(sprintf('1/500 in units of eta:  %f',(1/500)/eta));
disp(sprintf('1/512 in units of eta:            %f   %f', (1/512)/eta,eta*2*pi*512*sqrt(2)/3 )) ;
disp(sprintf('1/2048 in units of eta:  %f',(1/2048)/eta));

tturn=-2*ke./ke_diss_d;
tturn = tturn(length(tturn)/2:length(tturn));
tturn = sum(tturn)/length(tturn);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn));

epsilon=-ke_diss_d;
epsilon = epsilon(length(epsilon)/2:length(epsilon));
epsilon = sum(epsilon)/length(epsilon);
disp(sprintf('epsilon (averaged over last haf of data) = %f ',epsilon));





print -depsc scalars.ps












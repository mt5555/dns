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

%fid=endianopen('/ccs/scratch/taylorm/dns/iso12/iso12_256.scalars','r'); 
%nx=256;

%fid=endianopen('/scratch2/taylorm/sk256/sk2560000.0000.scalars','r'); 
%nx=256;

%fid=endianopen('/ccs/scratch/taylorm/dns/sc512A.scalars','r'); 
%fid2=endianopen('/ccs/scratch/taylorm/dns/iso12_512b.scalars','r'); 
%fid=fopen('../src/sht/rung0000.0000.scalars','r','l'); 

fid=fopen('/ccs/scratch/taylorm/decay/decay2048.scalars','r','l'); 
nx=2048;

%fid=endianopen('/ccs/taylorm/dns/src/temp0000.0000.scalars','r');
%nx=512;

%fid=endianopen('/ccs/scratch/taylorm/dns/sc1024A/sc1024A.scalars','r');
%nx=1024;

%fid=fopen('../src/sk128_alpha25/sk128_alpha250000.0000.scalars','r');
%nx=128;

%fid=endianopen('/scratch2/taylorm/tmix256B/tmix256B.scalars','r');
%nx=256;

%fid=endianopen('../../helicity_data/helical_forced/sk128_hq_0000.0000.scalars','r');
%nx = 128;

%fid=endianopen('../src/hel128_hpi4/hel128_hpi4_0000.0000.scalars','r');
%nx = 128;

%fid=endianopen('/home2/skurien/rotation/Test0000.0000.scalars','r');
%nx = 128;

%fid=endianopen('/home2/skurien/dns/src/sk128_alpha00/v5e-4/sk128_alpha000000.0000.scalars','r');
%nx = 128;

%fid=endianopen('/home2/skurien/rotation/Test/Test0000.0000.scalars','r');
%nx = 128;

%fid=endianopen('/home2/skurien/rotation/Test2/Test0002.scalars','r');
%nx = 128;

%fid=endianopen('/home2/skurien/dns/src/sk128_alpha40/sk128_v5e-4_alpha400000.0000.scalars','r');
%nx = 128;

%fid=endianopen('/home2/skurien/rotation/Rot1/Rot10000.0000.scalars','r');
%nx = 128;

%fid=endianopen('/home/taylorm/ccs/dns/src/rot3d/rot3d_sto0000.0000.scalars','r');
%nx=128;

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
  
[nints,ns,nscalars];
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

if (ns>=11) 
   h_diss=ints(11,:);   
end
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

ke_v=ke + alpha^4*ints(1,:)/2 + 2*alpha^2*ints(2,:)/2;

%Euse=ke; 
Euse=Ea;
%Euse=ke_v;

epsilon=-(  ke_diss_d-mu*alpha^2*ints(1,:) );
lambda=sqrt( mu*(2*Euse/3) / (epsilon/15)  );
R_l = lambda.*sqrt(2*Euse/3)/mu;
R_large = 1 * sqrt(2*Euse) / mu;
eta = (mu^3 ./ abs(epsilon)).^(.25);

epsilon_ke=-ke_diss_d;
lambda_ke=sqrt( mu*(2*ke/3) / (epsilon_ke/15)  );
R_l_ke=lambda_ke.*sqrt(2*ke/3)/mu;


disp(sprintf('max vor_z = %e',max(vor_z)));

figure(5)
clf
hold on
plot(time,ke,'k')
plot(time,ke_diss_f+ke_diss_d,'k')
plot(time_2,ke_diss_tot,'r')
plot(time,ke_diss_f,'b')
plot(time,ke_diss_d,'b')
plot(time,hel,'g')
title('KE: black,   \epsilon: blue  d(KE)/dt: red,    hel: green');
hold off
xlabel('time')

figure(6)
clf
semilogy(time,ke,'r'); hold on
%plot(time_2,-ke_diss_tot,'b.')
plot(time,-ke_diss_d,'b')
title('KE: red,    d(KE)/dt: blue');
hold off
xlabel('time')
print -dpsc ke.ps
print -djpeg -r72 ke.jpg
%fout=fopen('ke.out','wt');
%for i=1:length(time)
%n   fprintf(fout,'%.14f, %.14f, %.14f\n',time(i),ke(i),-ke_diss_d(i));
%end
%fclose(fout)


figure(7)
plot(time,maxvor,'r'); hold on;
%plot(.4188,2500,'o')
%plot(.4328,2500,'o')
%plot(.4603,2500,'o')
plot(.4894,2500,'.')
%plot(.5551,2500,'o')
%plot(.6034,2500,'o')
%plot(.8149,2500,'o')
plot(time,50000*ke,'k');
hold off;
%axis([0,1,0,5000]);
title('maximum vorticity component')



xlabel('time')
print -djpeg -r72 vor.jpg


figure(2);  hold on;
plot(time,R_l,'b'); hold on;
plot(time,R_l_ke,'r'); hold on;
title('R_\lambda');
legend('R_{\lambda}', 'R_{\lambda}(total KE)')
xlabel('time')
print -djpeg -r72 rl.jpg

figure(3);
plot(time,lambda)
     title('\lambda')
     xlabel('time')
print -djpeg -r72 lambda.jpg


figure(4); subplot(1,1,1)
plot(time,eta* nx*pi*2*sqrt(2)/3 )
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

% averge R_l_ke to a number
R_l_ke = R_l_ke(length(R_l_ke)/2:length(R_l_ke));
R_l_ke = sum(R_l_ke)/length(R_l_ke);
disp(sprintf('R_l_ke (average over last half of data) = %f ',R_l_ke));

disp(sprintf('1/250 in units of eta:  %f',(1/250)/eta));
disp(sprintf('1/500 in units of eta:  %f',(1/500)/eta));
disp(sprintf('1/512 in units of eta:            %f   %f', (1/512)/eta,eta*2*pi*512*sqrt(2)/3 )) ;
disp(sprintf('1/nx in units of eta:  %f',(1/nx)/eta));

tturn=-2*ke./ke_diss_d;
tturn = tturn(length(tturn)/2:length(tturn));
tturn = sum(tturn)/length(tturn);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn));


epsilon = epsilon(length(epsilon)/2:length(epsilon));
epsilon = sum(epsilon)/length(epsilon);
disp(sprintf('epsilon_E (averaged over last haf of data) = %f ',epsilon));

epsilon_ke = epsilon_ke(length(epsilon_ke)/2:length(epsilon_ke));
epsilon_ke = sum(epsilon_ke)/length(epsilon_ke);
disp(sprintf('epsilon_ke (averaged over last haf of data) = %f ',epsilon_ke));

print -depsc scalars.ps


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
LZ=1; %default aspect ratio is 1

fid2=-1;


%fid = endianopen('/netscratch/skurien/dns_data/sc1024A/sc1024A.scalars','r')
%nx=1024; f_k=0;fcor=0


%fid = endianopen('/netscratch/skurien/dns_data/sc2048A/sc2048A.scalars','r')
%nx=2048; f_k=0;fcor=0


%fid = endianopen('~/INCITE_runs/SW02_tests/bous128_Ro21Fr0.21/bous128_Ro21Fr0.21_all.scalars','r')
nx=128; f_k=24;fcor=1.07;bous=107.08;

%fid = endianopen('~/INCITE_runs/SW02_tests/bous128_Ro2.1Fr0.21/bous128_Ro2.1Fr0.21_all.scalars','r')
%nx=128; f_k=24;fcor=10.7;bous=107.08;

%fid = endianopen('~/INCITE_runs/SW02_tests/bous128_Fr0.21/bous128_Fr0.21_all.scalars','r')
%nx=128; f_k=24;fcor=0;bous=107.08;

%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/RSSlowres0000.0000.scalars','r')
%nx=400; f_k=4;fcor=136.2;bous=136.2;

%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/lowres10000.0000.scalars','r')
%nx=200; f_k=4;fcor=21.68;bous=108.4;LZ=0.2;

%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/lowres20000.0000.scalars','r')
%nx=200; f_k=4;fcor=136;bous=681;LZ=0.2;

%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/l400_d0.1_0000.0000.scalars','r')
%nx=400; f_k=4;fcor=136;bous=681;LZ=0.2;

%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/lowres30000.0000.scalars','r')
%nx=200; f_k=4;fcor=27.2;bous=136;LZ=0.2;

%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/lowres/l3_400d0.1_0000.0000.scalars','r')
%nx=400; f_k=4;fcor=27.2;bous=136;LZ=0.2;



%fid = endianopen('~/INCITE_runs/RemSukSmi09_tests/highres/RSShighres_all.scalars','r')
%nx=500; f_k=4;fcor=136.2;bous=13; LZ=0.2

%fid = endianopen('~/INCITE_runs/Intrepid/RSS09_tests/aspect/aspect_newUd1_0000.0000.scalars','r')
%fid = endianopen('~/INCITE_runs/Intrepid/aspect/aspect0.2d100000.0000.scalars','r')
%nx=400; f_k=4;fcor=100;bous=500;LZ=0.2;

%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_NS/n256_d1_0000.0000.scalars','r')
%nx=256; f_k=4;fcor=0;bous=0;LZ=1.0;
%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_NS/n256_d0.5_0000.0000.scalars','r')
%nx=512; f_k=4;fcor=0;bous=0;LZ=0.5;

%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_NS/n256_d1_0000.0000.scalars','r')
%nx=256; f_k=4;fcor=0;bous=0;LZ=1.0;
%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_NS/n256_d0.5_0000.0000.scalars','r')
%nx=512; f_k=4;fcor=0;bous=0;LZ=0.5;
%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_NS/n256_d0.25_all.scalars','r')
%nx=1024; f_k=4;fcor=0;bous=0;LZ=0.25;

%fid = endianopen('~/INCITE_runs/Intrepid/lowaspect_bous/shift_force/n1600_d0.2_Ro0.05_all.scalars', 'r')
%nx=1600; f_k=10;fcor=185;bous=925;LZ=0.2;



%fid = endianopen('/home/kurien/INCITE_runs/Intrepid/bous_NSvisc/n256_Ro1Fr0.01_all.scalars','r')
%nx=256; f_k = 4; fcor=8.58; bous = 858; LZ = 1.0;


%fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/n512_Ro1Fr0.01_all.scalars','r')
%nx=512; f_k = 4; fcor=8.58; bous = 858; LZ = 1.0;

fid = endianopen('~/projects/INCITE_runs/Intrepid/bous_NSvisc/n1024_Ro1Fr0.01_all.scalars','r')
nx=1024; f_k = 4; fcor=8.58; bous = 858; LZ = 1.0;


fid = endianopen('~/projects/INCITE_runs/Intrepid/qg/n640_bous3000_all.scalars','r')
nx=640; f_k=4;fcor=3000;bous=3000;LZ=1.0;

fid = endianopen('~/projects/INCITE_runs/Intrepid/Ro1Fr0/n640_fcor14bous3000_all.scalars','r')
nx=640; f_k=4;fcor=14;bous=3000; LZ=1.0;
%fid = endianopen('~/INCITE_runs/Intrepid/Ro0Fr1/n640_fcor3000bous14_all.scalars','r')
%nx=640; f_k=4;fcor=3000;bous=14; LZ=1.0;


nscalars=0;
nscalars_e=0;
time = 0.0;
%while (1) 
while (time < 500)
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
    ni1=size(ints); ni1=ni1(1);
    ni2=size(data1); ni2=ni2(1);
    if (ni2>ni1)
      % the input data is larger - this means that at some point during
      % the run,  we added attitional scalars.  pad out old data with NaN's
      ints=[ints;NaN*ints(1:ni2-ni1,:)];
      maxs=[maxs;NaN*maxs(1:ni2-ni1,:)];
    end
    ints=[ints,data1];
    maxs=[maxs,data2];
  end
  nscalars=nscalars+ns;

  % now read the "expensive" integrals, which are not computed every time step
  ns_e = fread(fid,1,'float64');
  time = fread(fid,1,'float64');
  data1 = fread(fid,[ns_e,1],'float64');
end

disp(sprintf('nints=%i  total scalars read=%i',nints,nscalars))
fclose(fid);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the scalars computed every time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=size(ints);
l=l(2);

if (nints>=11) 
   h_diss=ints(11,:);   
end
ke_diss_d=ints(10,:);
ke_diss_f=ints(3,:);        % < u,F >   F=f, except for alpha model F=f'
disp(sprintf('mean ke_diss_f = %f', mean(ke_diss_f)))

vor_z=ints(4,:);
hel=ints(5,:);
ke=ints(6,:);      %  ke 
ens = ints(7,:);   %  enstrophy
u3 = ints(12,:);   %  problem dependent.  for ns_xpencil.F90 = u^3
% ints(8,:)   < u,div(tau)' >  (alpha model only)
% ints(9,:)   < u,f>           (alpha model only)
% ints(1,:)  < u_xx,u_xx> >   (used for E_alpha dissapation term)

maxU=maxs(1,:);  % at time_after
maxV=maxs(2,:);  % at time_after
maxW=maxs(3,:);  % at time_after
maxUcfl=maxs(4,:);     % max used for CFL, at time_after
maxvor=maxs(5,:);
time_after=maxs(6,:);
time=maxs(7,:);
if (nints >=15)
pe = ints(15,:);
pe_diss = ints(16,:);
end


% $$$ % time  ke  
% $$$ for i=1:length(time)
% $$$    disp(sprintf('time=%7.4f   KE=%18.14f   ens=%18.12f',time(i),ke(i),ens(i)))
% $$$ end
% $$$ for i=1:length(time_after)
% $$$    disp(sprintf('time=%7.4f   maxU=%18.14f  %18.14f  %18.14f',time_after(i),maxU(i),maxV(i),maxW(i)))
% $$$ end%
% $$$ return


Ea = ints(6,:) + .5*alpha^2 *ints(2,:); % at time

time_2=[];
ke_diss_tot=[];
Ea_diss_tot=[];
h_diss_tot=[];
ke_diss_tot2=[];
j=0;
for i=2:l
  % skip repeats, and skip data with large KE increase
  % since it is probably due to decaying run re-adjustment:
  if ((time(i)-time(i-1))>1e-5 & (ke(i)-ke(i-1))<.01  )
    j=j+1;
    time_2(j) = .5*(time(i)+time(i-1));
    ke_diss_tot(j)=(ke(i)-ke(i-1))./(time(i)-time(i-1));
    pe_diss_tot(j)=(pe(i)-pe(i-1))./(time(i)-time(i-1));
    Ea_diss_tot(j)=(Ea(i)-Ea(i-1))./(time(i)-time(i-1));
    h_diss_tot(j)=(hel(i)-hel(i-1))./(time(i)-time(i-1));
    ke_diss_tot2(j)=.5*(ke_diss_d(i)+ke_diss_d(i-1)) + ...
                    .5*(ke_diss_f(i)+ke_diss_f(i-1));
  end
end

ke_v=ke + alpha^4*ints(1,:)/2 + 2*alpha^2*ints(2,:)/2;

%Euse=ke; 
Euse=Ea;
%Euse=ke_v;

epsilon=-(  ke_diss_d-mu*alpha^2*ints(1,:) );
lambda=sqrt( mu*(2*Euse/3) ./ (epsilon/15)  );
eta = (mu^3 ./ abs(epsilon)).^(.25);

epsilon_ke=-ke_diss_d;
epsilon_pe=-pe_diss;
lambda_ke=sqrt( mu*(2*ke/3) ./ (epsilon_ke/15)  );

if (mu>0)
  R_l = lambda.*sqrt(2*Euse/3)/mu;
  R_large = 1 * sqrt(2*Euse) / mu;
  R_l_ke=lambda_ke.*sqrt(2*ke/3)/mu;
end



disp(sprintf('max vor_z = %e',max(vor_z)));

if (mu>0)
  figure(2);  
  plot(time,R_l,'b'); hold on;
  plot(time,R_l_ke,'r'); hold on;
  hold off;
  title('R_\lambda');
  legend('R_{\lambda}', 'R_{\lambda}(total KE)')
  xlabel('time')
  print -djpeg -r72 rl.jpg
end


%figure(3);
% $$$ plot(time,lambda)
% $$$ title('\lambda')
% $$$ xlabel('time')
% $$$ print -djpeg -r72 lambda.jpg





if (nx>0) 
  figure(4); subplot(1,1,1)
  plot(time,eta* nx*pi*2*sqrt(2)/3 )
  title('k_{nmax} \eta')
  xlabel('time');hold on;
  %print -djpeg -r72 kmaxeta.jpg
end

%plot energies as a function of nonlinear times (SW02 paper)
figure(5); 
ke_diss_f_ave = ke_diss_f(ceil(length(ke_diss_f)/2):length(ke_diss_f));
ke_diss_f_ave = sum(ke_diss_f_ave)/length(ke_diss_f_ave);
tn = (ke_diss_f_ave*(2*pi*f_k)^2)^(1/3); %(fcor + bous)/2/pi; 
tl = (1*(4)^2)^(1/3); %smith's time
en = ((ke_diss_f_ave/(2*pi*f_k))^(-2/3));
el = ((1/4))^(-2/3);   %smith's energy 
timen = time*tn/tl;
ken = ke*en/el;
pen = pe*en/el;
plot(timen,ken,'k')
hold on
plot(timen,pen,'ro')
plot(timen,ken+pen,'b')
%plot(time_2,ke_diss_tot,'r')
%plot(time,hel,'g')
title('KE: black, PE: red, Etot: blue');
%hold off
xlabel('time t(\epsilon_f k_f^2)^{1/3}')
ylabel('E(\epsilon_f/k_f)^{-2/3}')



% look at dissipations seperatly
figure(6)
%clf
hold on
plot(time_2,ke_diss_tot,'k')
plot(time,ke_diss_f,'r')
plot(time,ke_diss_d,'g')
%plot(time,ke_diss_f+ke_diss_d,'b')
%plot(time_2,ke_diss_tot2-ke_diss_tot,'g')

title('F: red  D: green  d(KE)/dt: black');
hold off
xlabel('time')
%axis([.2 .4 -6 6]); return


% look at dissipations seperatly
figure(7)
%clf
plot(time_2,ke_diss_tot2-ke_diss_tot,'k')
hold on
title(' (F+D) - d(KE)/dt');
hold off
xlabel('time')




figure(8)
plot(time,maxvor,'r'); hold on;
%plot(.4188,2500,'o')
%plot(.4328,2500,'o')
%plot(.4603,2500,'o')
%plot(.4894,2500,'.')
%plot(.5551,2500,'o')
%plot(.6034,2500,'o')
%plot(.8149,2500,'o')
%plot(time,50000*ke,'k');
%hold off;
%axis([0,1,0,5000]);
title('maximum vorticity component')
xlabel('time');
%print -djpeg -r72 vor.jpg

%plot helicity, enstrophy and ratio H/W
figure(10)
subplot(3,1,1);
plot(timen,hel,'r');grid on;
ylabel('H');
subplot(3,1,2);
plot(timen,ens,'b');grid on;
ylabel('W');
subplot(3,1,3);
plot(timen,hel./ens,'k');grid on;
ylabel('H/W');




tturn=-2*ke./ke_diss_d;
tturn = tturn(ceil(length(tturn)/2):length(tturn));
tturn = sum(tturn)/length(tturn);
disp(sprintf('eddy turnover time (averaged over last half of data) = %f ',tturn));


maxU = maxU(ceil(length(maxU)/2):length(maxU));
maxU = sum(maxU)/length(maxU);
disp(sprintf('maxU (average over last half of data) = %f ',maxU));
maxV = maxV(ceil(length(maxV)/2):length(maxV));
maxV = sum(maxV)/length(maxV);
disp(sprintf('maxV (average over last half of data) = %f ',maxV));
maxW = maxW(ceil(length(maxW)/2):length(maxW));
maxW = sum(maxW)/length(maxW);
disp(sprintf('maxW (average over last half of data) = %f ',maxW));

maxUcfl=maxUcfl/nx;
maxUcfl = maxUcfl(ceil(length(maxUcfl)/2):length(maxUcfl));
maxUcfl = sum(maxUcfl)/length(maxUcfl);
disp(sprintf('maxUcfl (average over last half of data) = %f ',maxUcfl));

ke = ke(ceil(length(ke)/2):length(ke));
ke = sum(ke)/length(ke);
disp(sprintf('ke (average over last half of data) = %f ',ke));

%pe = pe(ceil(length(pe)/2):length(pe));
%pe = sum(pe)/length(pe);
%disp(sprintf('pe (average over last half of data) = %f ',pe));

%disp(sprintf('Total energy (average over last half of data) = %f ',ke+pe));

if (mu>0) 
  % averge eta to a number
  eta = eta(ceil(length(eta)/2):length(eta));
  eta = sum(eta)/length(eta);
  disp(sprintf('eta (average over last half of data) = %f ',eta));

% averge R_l to a number
  R_l = R_l(ceil(length(R_l)/2):length(R_l));
  R_l = sum(R_l)/length(R_l);
  disp(sprintf('R_l (average over last half of data) = %f ',R_l));
  
  % averge R_l_ke to a number
  R_l_ke = R_l_ke(ceil(length(R_l_ke)/2):length(R_l_ke));
  R_l_ke = sum(R_l_ke)/length(R_l_ke);
  disp(sprintf('R_l_ke (average over last half of data) = %f ',R_l_ke));

  disp(sprintf('1/250 in units of eta:  %f',(1/250)/eta));
  disp(sprintf('1/500 in units of eta:  %f',(1/500)/eta));
  disp(sprintf('1/512 in units of eta:            %f   %f', (1/512)/eta,eta*2*pi*512*sqrt(2)/3 )) ;
  disp(sprintf('1/nx in units of eta:  %f',(1/nx)/eta));
end


if (alpha>0)
  epsilon = epsilon(ceil(length(epsilon)/2):length(epsilon));
  epsilon = sum(epsilon)/length(epsilon);
  disp(sprintf('epsilon_E (averaged over last half of data) = %f ',epsilon));
end

epsilon_ke = epsilon_ke(ceil(length(epsilon_ke)/2):length(epsilon_ke));
epsilon_ke = sum(epsilon_ke)/length(epsilon_ke);
disp(sprintf('epsilon_ke (averaged over last half of data) = %f ',epsilon_ke));

epsilon_pe = epsilon_pe(ceil(length(epsilon_pe)/2):length(epsilon_pe));
epsilon_pe = sum(epsilon_pe)/length(epsilon_pe);
disp(sprintf('epsilon_pe (averaged over last half of data) = %f ',epsilon_pe));


ke_diss_f = ke_diss_f(ceil(length(ke_diss_f)/2):length(ke_diss_f));
ke_diss_f = sum(ke_diss_f_ave)/length(ke_diss_f_ave);
disp(sprintf('epsilon_f (averaged over last half of data) = %f ',ke_diss_f));
Ro = (ke_diss_f * (2*pi*f_k)^2 ).^(1/3) / (fcor);
Fr = (ke_diss_f * (2*pi*f_k)^2 ).^(1/3) / (bous);
disp(sprintf('Ro computed from epsilon_f = %f ',Ro));
disp(sprintf('Fr computed from epsilon_f = %f ',Fr));

U = ( ke_diss_f / (f_k*(1/LZ)*2*pi) ).^(1/3);
L = 1/(2*pi*f_k);
H = LZ*L;
Ro = U/(fcor*L);
Fr = U/(bous*H);
disp(sprintf('Ro computed from LZ & epsilon_f = %f ',Ro));
disp(sprintf('Fr computed from LZ & epsilon_f = %f ',Fr));


ke_diss_tot = ke_diss_tot(ceil(length(ke_diss_tot)/2):length(ke_diss_tot));
ke_diss_tot = sum(ke_diss_tot)/length(ke_diss_tot);
disp(sprintf('d(KE)/dt (averaged over last half of data) = %f ',ke_diss_tot));

pe_diss_tot = pe_diss_tot(ceil(length(pe_diss_tot)/2):length(pe_diss_tot));
pe_diss_tot = sum(pe_diss_tot)/length(pe_diss_tot);
disp(sprintf('d(PE)/dt (averaged over last half of data) = %f ',pe_diss_tot));



print -depsc scalars.ps


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


range=0:50;
%range = 1
%range=0:.5:2;
fid2=-1;


fid=fopen('/home/mataylo/codes/dns_data/decay/decay2048.scalars','r','l'); 
nx=2048;

%fid=endianopen('/nh/nest/u/skurien/projects/helicity_data/helical_forced/hel512_hpi2/diag/skhel512_hpi2.scalars','r');
%nx=512;

%fid=endianopen('/nh/nest/u/taylorm/sc1024A.scalars','r');
%nx=1024;

%fid=endianopen('/home/mataylo/codes/kerr/test/temp0000.0000.scalars','r');


%fid=endianopen('/home/mataylo/all.scalars','r');
%f_k= 16; fcor=40;



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

% $$$ % time  ke  
% $$$ for i=1:length(time)
% $$$    disp(sprintf('time=%7.4f   KE=%18.14f   ens=%18.12f',time(i),ke(i),ens(i)))
% $$$ end
% $$$ for i=1:length(time_after)
% $$$    disp(sprintf('time=%7.4f   maxU=%18.14f  %18.14f  %18.14f',time_after(i),maxU(i),maxV(i),maxW(i)))
% $$$ end
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
lambda=sqrt( mu*(2*Euse/3) / (epsilon/15)  );
eta = (mu^3 ./ abs(epsilon)).^(.25);

epsilon_ke=-ke_diss_d;
lambda_ke=sqrt( mu*(2*ke/3) / (epsilon_ke/15)  );

if (mu>0)
  R_l = lambda.*sqrt(2*Euse/3)/mu;
  R_large = 1 * sqrt(2*Euse) / mu;
  R_l_ke=lambda_ke.*sqrt(2*ke/3)/mu;
end




disp(sprintf('max vor_z = %e',max(vor_z)));

if (mu>0)
  figure(2);  hold on;
  plot(time,R_l,'b'); hold on;
  plot(time,R_l_ke,'r'); hold on;
  title('R_\lambda');
  legend('R_{\lambda}', 'R_{\lambda}(total KE)')
  xlabel('time')
  print -djpeg -r72 rl.jpg
end

figure(3);
plot(time,lambda)
     title('\lambda')
     xlabel('time')
print -djpeg -r72 lambda.jpg


if (nx>0) 
  figure(4); subplot(1,1,1)
  plot(time,eta* nx*pi*2*sqrt(2)/3 )
  title('k_{nmax} \eta')
  xlabel('time')
  print -djpeg -r72 kmaxeta.jpg
end

figure(5)
clf
hold on
plot(time,ke,'k')
%plot(time_2,ke_diss_tot,'r')
%plot(time,hel,'g')
title('KE: black,   hel: green');
hold off
xlabel('time')


% look at dissipations seperatly
figure(6)
clf
hold on
plot(time_2,ke_diss_tot,'k')
plot(time,ke_diss_f,'r')
plot(time,ke_diss_d,'g')
%plot(time,ke_diss_f+ke_diss_d,'b')
%plot(time_2,ke_diss_tot2-ke_diss_tot,'g')

title('F: red  D: green  d(KE)/dt: black');
hold off
xlabel('time')


% look at dissipations seperatly
figure(7)
clf
hold on
plot(time_2,ke_diss_tot2-ke_diss_tot,'k')
title(' (F+D) - d(KE)/dt');
hold off
xlabel('time')




figure(8)
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
xlabel('time');
print -djpeg -r72 vor.jpg









tturn=-2*ke./ke_diss_d;
tturn = tturn(ceil(length(tturn)/2):length(tturn));
tturn = sum(tturn)/length(tturn);
disp(sprintf('eddy turnover time (averaged over last haf of data) = %f ',tturn));


% averge Smith's R0 to a number
if (f_k>0) 
  R0 = (ke_diss_f * (2*pi*f_k)^2 ).^(1/3) / (.5 * fcor);
  R0 = R0(ceil(length(R0)/2):length(R0));
  R0 = sum(R0)/length(R0);
  disp(sprintf('Smith''s R0 (average over last half of data) = %f ',R0));
end


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
  disp(sprintf('epsilon_E (averaged over last haf of data) = %f ',epsilon));
end

epsilon_ke = epsilon_ke(ceil(length(epsilon_ke)/2):length(epsilon_ke));
epsilon_ke = sum(epsilon_ke)/length(epsilon_ke);
disp(sprintf('epsilon_ke (averaged over last haf of data) = %f ',epsilon_ke));

ke_diss_f = ke_diss_f(ceil(length(ke_diss_f)/2):length(ke_diss_f));
ke_diss_f = sum(ke_diss_f)/length(ke_diss_f);
disp(sprintf('epsilon_f (averaged over last haf of data) = %f ',ke_diss_f));

ke_diss_tot = ke_diss_tot(ceil(length(ke_diss_tot)/2):length(ke_diss_tot));
ke_diss_tot = sum(ke_diss_tot)/length(ke_diss_tot);
disp(sprintf('d(KE)/dt (averaged over last haf of data) = %f ',ke_diss_tot));


print -depsc scalars.ps


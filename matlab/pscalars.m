%
%   matlab script to read DNS output files *.scalars-turb
%
%    
%

clear all;

name = '/scratch2/taylorm/tmix256C/tmix256C'
times=[1.0000:.05:1.2];


% load in the u,v,w scalars:
use_pscalars_name=1;
disp('running scalars_turb...')
scalarsturb
disp('reading passive scalar data...')

nt=0;
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.pscalars-turb'];
  disp(fname)
  fid=endianopen(fname,'r')
  if (fid>=0) 

    [ns_e,count] = fread(fid,1,'float64');
    [npassive,count] = fread(fid,1,'float64');
    time = fread(fid,1,'float64');
    mu = fread(fid,1,'float64');

    nt=nt+1;
    if(nt==1) pints_e=zeros([2+ns_e,npassive,1]); end;

    for np=1:npassive
       data1 = fread(fid,[ns_e,1],'float64');
       data1=[time;mu;data1];
       pints_e(:,np,nt)= data1;
    end
    fclose(fid);
    [nt,ns_e]
  end
end


% look at passive scalar number np
np=1;

time_e=squeeze(pints_e(1,np,:))';
mu=squeeze(pints_e(2,np,:))';
schmidt=squeeze(pints_e(3,np,:)/128)';
c2=squeeze(pints_e(4,np,:))';

cx2(1,:)=pints_e(5,np,:);
cx2(2,:)=pints_e(6,np,:);
cx2(3,:)=pints_e(7,np,:);
cx3(1,:)=pints_e(8,np,:);
cx3(2,:)=pints_e(9,np,:);
cx3(3,:)=pints_e(10,np,:);

cxx2(1,:)=pints_e(14,np,:);
cxx2(2,:)=pints_e(15,np,:);
cxx2(3,:)=pints_e(16,np,:);

cu(1,:)=pints_e(23,np,:);
cu(2,:)=pints_e(24,np,:);
cu(3,:)=pints_e(25,np,:);




%
% isotropic relations (may not exactly agree with data from scalars.m)
%
% data needed from scalarsturb.m:
%       u2, ux2, ux3, uxx2
%
%


ke=.5*sum(u2,1);
epsilon=15*mu.*mean(ux2,1);                   % only uses ux^2, vy^2, wz^2
lambda=sqrt( mu.*(2*ke/3) ./ (epsilon/15)  );
R_l = lambda.*sqrt(2*ke/3)./mu;
Rt= R_l*R_l*3/20;
eta = mu^3/epsilon;



disp(sprintf('mu = %f',mu ))
disp(sprintf('ke = %f',ke ))
disp(sprintf('R_l = %f',R_l ))
disp(sprintf('lambda = %f',lambda ))
disp(sprintf('epsilon = %f',epsilon ))

epsilon_c=3*(mu./schmidt).*mean(cx2,1);
lambda_c=sqrt(c2./mean(sx2));
eta_c=eta/sqrt(schmidt);



disp(' ')
disp(sprintf('passive scalar n=%i',np))
disp(sprintf('schmidt = %f',schmidt ))
disp(sprintf('s2 = %f',c2 ))
disp(sprintf('lambda_c = %f',lambda_s ))
disp(sprintf('epsilon_c = %f',epsilon_s ))


%
%data for Ray:
% 
% average over x,y,z directions:


Su = mean(ux3,1)/mean(ux2,1)^1.5  ;
Suc = -mean(cu,1)/sqrt(mean(ux2,1))/mean(cx2,1)
G = mean(u2,1) * mean(uxx2,1)/mean(ux2)^2;
Gc = mean(c2,1) * mean(cxx2,1)/mean(cx2)^2;


ff = Su*sqrt(Rt)*7/3/sqrt(15);
ff = ff + G*7/15;

r=3*(lambda/lambda_c)^2 /5 / schmidt; 
gg = sqrt(5/3)*Suc*sqrt(Rt) + r*Gc;











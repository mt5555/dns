%
%   matlab script to read DNS output files *.scalars-turb
%
%    
%

readdata=0;
%clear all; readdata=1;



name = '/scratch2/taylorm/tmix256C/tmix256C'
times=[1.0000:.01:1.24];



if (readdata)
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
end

% look at passive scalar number np
for np=10:10

time_e=squeeze(pints_e(1,np,:))';
mu=squeeze(pints_e(2,np,:))';
schmidt=squeeze(pints_e(3,np,:))';   % index=1 from fortran data file
c2=squeeze(pints_e(4,np,:))';        % index=2 

cx2(1,:)=pints_e(5,np,:);
cx2(2,:)=pints_e(6,np,:);
cx2(3,:)=pints_e(7,np,:);
cx3(1,:)=pints_e(8,np,:);
cx3(2,:)=pints_e(9,np,:);
cx3(3,:)=pints_e(10,np,:);
% cx4                       % 11,12,13

cxx2(1,:)=pints_e(14,np,:);
cxx2(2,:)=pints_e(15,np,:);
cxx2(3,:)=pints_e(16,np,:);
% cxx3                        17,18,19
% cxx4                        20,21,22

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
Rt= R_l.*R_l*3/20;
eta = mu.^3./epsilon;



disp(sprintf('mu = %f',mu(1) ))
disp(sprintf('ke = %f',ke(1) ))
disp(sprintf('R_l = %f',R_l(1) ))
disp(sprintf('lambda = %f',lambda(1) ))
disp(sprintf('epsilon = %f',epsilon(1) ))

epsilon_c=3*(mu./schmidt).*mean(cx2,1);
lambda_c=sqrt(c2./mean(cx2));
eta_c=eta/sqrt(schmidt);



disp(' ')
disp(sprintf('passive scalar n=%i',np))
disp(sprintf('schmidt = %f',schmidt(1) ))
disp(sprintf('s2 = %f',c2(1) ))
disp(sprintf('lambda_c = %f',lambda_c(1) ))
disp(sprintf('epsilon_c = %f',epsilon_c(1) ))


%
%data for Ray:
% 
% average over x,y,z directions:


Su = mean(ux3,1)./mean(ux2,1).^1.5  ;
Suc = -mean(cu,1)./(sqrt(mean(ux2,1)).*mean(cx2,1));
G = (mean(u2,1) .* mean(uxx2,1))./mean(ux2).^2;
Gc = (c2 .* mean(cxx2,1))./mean(cx2).^2;


ff = Su.*sqrt(Rt)*7/3/sqrt(15);
ff = ff + G*7/15;

r=3* ((lambda./lambda_c).^2) /5 ./ schmidt; 
gg = sqrt(5/3)*Suc.*sqrt(Rt) + r.*Gc;


figure(1); clf;

subplot(5,2,1)
plot(time_e,Su)
title('S_u')

subplot(5,2,2)
plot(time_e,Suc)
title('S_{u\theta}')

subplot(5,2,3)
plot(time_e,G)
title('G')

subplot(5,2,4)
plot(time_e,Gc)
title('G_\theta')



subplot(5,2,5)
plot(time_e,ff)
title('f')

subplot(5,2,6)
plot(time_e,gg)
title('g')


subplot(5,2,7)
plot(time_e,lambda)
title('\lambda')

subplot(5,2,8)
plot(time_e,lambda_c)
title('\lambda_\theta')

subplot(5,2,9)
plot(time_e,r)
title('r')

subplot(5,2,10)
plot(time_e,R_l)
title('R_\lambda')

orient tall
pname=sprintf('pscalars%i',np);
print('-djpeg','-r90',pname); 

end



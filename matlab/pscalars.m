%
%   matlab script to read DNS output files *.scalars-turb
%
% each file contains only IEEE 8 byte floats, in this order:
%
%   number of data for each passive scalar
%   number of passive scalars
%   time
%
% then, repeated for each passive scalar:
%
%   mu                          2  (matlab array starts at 2 since we add time)
%   schmidt                     3 
%   s^2                         4
%   s,i ^2     i=1,3            5
%   s,i ^3     i=1,3            8
%   s,i ^4     i=1,3            11
%   s,ii ^2     i=1,3           14
%   s,ii ^3     i=1,3           17
%   s,ii ^4     i=1,3           20
%   ui,i*s,i    i=1,3           23
%    
%

clear all;

name = '/scratch2/taylorm/tmix256C/tmix256C'
times=[1.0001];

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
    if(nt==1) pints_e=zeros([1+ns_e,npassive,1]); end;

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

time_e=pints_e(1,np,:);
mu=pints_e(2,np,:);
schmidt=pints_e(3,np,:);
s2=pints_e(4,np,:);

sx2(1)=pints_e(5,np,:);
sx2(2)=pints_e(6,np,:);
sx2(3)=pints_e(7,np,:);
sx3(1)=pints_e(8,np,:);
sx3(2)=pints_e(9,np,:);
sx3(3)=pints_e(10,np,:);

sxx2(1)=pints_e(14,np,:);
sxx2(2)=pints_e(15,np,:);
sxx2(3)=pints_e(16,np,:);

su(1)=pints_e(23,np,:);
su(2)=pints_e(24,np,:);
su(3)=pints_e(25,np,:);


%
% isotropic relations (may not exactly agree with data from scalars.m)
%
% data needed from scalarsturb.m:
%       u2, ux2, ux3, uxx2
%
%


ke=.5*sum(u2,1);
epsilon=15*mu*mean(ux2,1);                   % only uses ux^2, vy^2, wz^2
lambda=sqrt( mu*(2*ke/3) / (epsilon/15)  );
R_l = lambda.*sqrt(2*ke/3)/mu;

disp(sprintf('mu = %f',mu ))
disp(sprintf('ke = %f',ke ))
disp(sprintf('R_l = %f',R_l ))
disp(sprintf('lambda = %f',lambda ))
disp(sprintf('epsilon = %f',epsilon ))

ke_s=.5*sum(s2,1);
epsilon_s=15*(mu/schmidt)*sum(s2,1)/3;  
lambda_s=sqrt(mean(s2)/mean(sx2));






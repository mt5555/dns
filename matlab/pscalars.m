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
%   mu                          1
%   schmidt                     2 
%   s^2                         3
%   s,i ^2     i=1,3            4-6
%   s,i ^3     i=1,3            7-9
%   s,i ^4     i=1,3            10-12
%   s,ii ^2     i=1,3           13-15
%   s,ii ^3     i=1,3           16-18
%   s,ii ^4     i=1,3           19-21
%   ui,i*s,i    i=1,3           22-24
%    
%

clear all;

name = '/ccs/scratch/taylorm/decay/decay2048'
times=[0,.0167,.0283,.0524,.0650,.0776, .2131  .2207 .2357]
times=[times,.24:.01:4.0];

%name = '/ccs/scratch/taylorm/dns/iso12_512'
%times=7;

%name = '/ccs/scratch/taylorm/dns/sc512A'
%times=3;

%name = '/home/skurien/dns/src/sk128_alpha0';
%times=0.5

name = '/scratch2/taylorm/tmix256B/tmix256B'
times=[0:.5:3]


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

    nt=nt+1;
    if(nt==1) pints_e=zeros([1+ns_e,npassive,1]); end;

    for np=1:npassive
       data1 = fread(fid,[ns_e,1],'float64');
       data1=[time;data1];
       pints_e(:,np,nt)= data1;
    end
    fclose(fid);
    [nt,ns_e]
  end
end


% look at passive scalar number np
np=1;

mu=pints_e(1,np,:);
schmidt=pints_e(2,np,:);
time_e=pints_e(3,np,:);
s2(1)=pints_e(4,np,:);
s2(2)=pints_e(5,np,:);
s2(3)=pints_e(6,np,:);

sx2(1)=pints_e(7,np,:);
sx2(2)=pints_e(8,np,:);
sx2(3)=pints_e(9,np,:);

sxx2(1)=pints_e(13,np,:);
sxx2(2)=pints_e(14,np,:);
sxx2(3)=pints_e(15,np,:);

su(1)=pints_e(22,np,:);
su(2)=pints_e(23,np,:);
su(3)=pints_e(24,np,:);


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


ke_s=.5*sum(s2,1);
epsilon_s=15*(mu/schmidt)*sum(s2,1)/3;  
lambda_s=sqrt(mean(s2)/mean(sx2));






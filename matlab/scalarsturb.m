%
%   matlab script to read DNS output files *.scalars-turb
%
% each file contains only IEEE 8 byte floats, in this order:
%
%   number of scalars in file (not including time). should be 23.0
%   time
%   u_i,i^2    i=1,3
%   u_i,i^3    i=1,3
%   u_i,i^4    i=1,3
%   wSw
%   u_i^2      i=1,3 
%   S2
%   w_i^2      i=1,3
%   w_i^3      i=1,3
%   w_i^4      i=1,3
%

if (exist('use_pscalars_name') & use_pscalars_name)
   % name set in pscalars.m
else

   %name = '/home/mataylo/data/dns/decay/decay2048'
   %times=[0,.0167,.0283,.0524,.0650,.0776, .2131  .2207 .2357]
   %times=[times,.24:.01:4.0];

   %name = '/ccs/scratch/taylorm/dns/iso12_512'
   %times=7;

   %name = '/ccs/scratch/taylorm/dns/sc512A'
   %times=3;

   %name = '/home/skurien/dns/src/sk128_alpha0';
   %times=0.5

   %name='/scratch2/taylorm/tmix256D-noscalars/tmix256D-noscalars';
   %times=0:.1:.5;

   %name = '/home/skurien/dns/src/sk128_alpha15/sk128_alpha15';
   %times=[1,2,3]

   name = '~/INCITE_runs/Intrepid/lowaspect_NS/n256_d1_';
   times=[0:.1:2.0];

end

nt=0;
for t=times
  tstr=sprintf('%10.4f',t+10000);
  fname=[name,tstr(2:10),'.scalars-turb'];
  disp(fname)
  fid=fopen(fname,'r','l');
  if (fid>=0) 

    [ns_e,count] = fread(fid,1,'float64');
    
    time = fread(fid,1,'float64');
    data1 = fread(fid,[ns_e,1],'float64');
    data1=[time;data1];
    if (nt==0) 
      ints_e= data1;
    else
      ni1=size(ints_e); ni1=ni1(1);
      ni2=size(data1); ni2=ni2(1);
      if (ni2>ni1)
        % the input data is larger - this means that at some point during
        % the run,  we added attitional scalars.  pad out old data with NaN's
        ints_e=[ints_e;NaN*ints_e(1:ni2-ni1,:)];
      end
      if (ni1>ni2)
        % the input data is smaller - this means that at some point during
        % the run,  we removed attitional scalars.  pad out input NaN's
        data1=[data1;NaN*ni2:ni1-1];
      end
      ints_e= [ints_e,data1];
    end
    nt=nt+1;
    [nt,ns_e]
    fclose(fid);
  else
    disp('skipping file...')
  end
end




time_e=ints_e(1,:);
ux2=zeros([3,nt]);
ux3=ux2;
ux4=ux2;
for i=1:3
   ux2(i,:)=ints_e(i+1,:);    % < u1,1 ^ 2 >
   ux3(i,:)=ints_e(i+4,:);    % < u1,1 ^ 3 >
   ux4(i,:)=ints_e(i+7,:);    % < u1,1 ^ 4 >
end
Sww=ints_e(11,:);
if (ns_e>=13) 
  for i=1:3
    u2(i,:)=ints_e(i+11,:);              % < u^2 >
  end
end
if (ns_e>=14) 
    S2=ints_e(15,:);                     % <S2>
end
if (ns_e>=23) 
for i=1:3
   vor2(i,:)=ints_e(i+15,:);    % < w^2 > 
   vor3(i,:)=ints_e(i+18,:);    % < w^3 >
   vor4(i,:)=ints_e(i+21,:);    % < w^4 >
end
end
% S4sum = ints_e(25,:)
% S2w2  = ints_e(26,:)

if (ns_e>=28) 
for i=1:3
   uxx2(i,:)=ints_e(i+26,:);    % < uxx^2 > 
end
end
if (ns_e>=37)
  for i=1:3; u1(i,:)=ints_e(i+29,:);   end;   % < u1 >
  for i=1:3; u3(i,:)=ints_e(i+32,:);   end;   % < u3 >
  for i=1:3; u4(i,:)=ints_e(i+35,:);   end;   % < u4 >
end


% zero crossing
if (ns_e>=43)
  for i=1:3; nu(i,:)=ints_e(i+38,:); end;  % N0 for (u,v,w), long. direction
  for i=1:3; nux(i,:)=ints_e(i+38,:); end; % N0 for ux,vy,wz, long. direction
end


figure(6)
clf
%plot(time_e,Sww)

plot(time_e,ux3(1,:)./ux2(1,:).^(3/2),'k')
hold on
plot(time_e,ux3(2,:)./ux2(2,:).^(3/2),'b')
plot(time_e,ux3(3,:)./ux2(3,:).^(3/2),'g')

%plot(time_e,ux2(1,:),'k')
%plot(time_e,ux2(2,:),'b')
%plot(time_e,ux2(3,:),'g')


%plot(time_e,ux4(1,:)./ux2(1,:).^2,'r')
%plot(time_e,ux4(2,:)./ux2(2,:).^2,'b')
%plot(time_e,ux4(3,:)./ux2(3,:).^2,'g')

%plot(time_e,vor4(1,:)./vor2(1,:).^2,'r')
%hold on
%plot(time_e,vor4(2,:)./vor2(2,:).^2,'b')
%plot(time_e,vor4(3,:)./vor2(3,:).^2,'g')

ax=axis;
axis([ax(1),ax(2),-.6,-.3]);
print -djpeg -r72 skew.jpg

hold off















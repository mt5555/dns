%
%   matlab script to read DNS output files *.scalars-turb
%
%    
%

readdata=0;
%clear all; readdata=1;


if (readdata==0 & ~exist('cxx2'))
   % force reading of new data anyway
   clear all; readdata=1;
end



name = '/scratch2/taylorm/tmix256C/tmix256C'
times=[.75:.01:3.5];


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
  
  pdfname=[name,tstr(2:10),'.spdf']
  fid=fopen(pdfname,'r','l');
  if (fid>=0)
    time=fread(fid,1,'float64');
    npmax=fread(fid,1,'float64');         
    np=1;
    for p=1:npmax
      [n_del,delta,bin_size,n_bin,n_call,bins,pdf]=read1pdf(fid);
      vtot = sum(pdf);
      Vpdata25(p,nt) = sum(pdf.*(bins>=.25 & bins <=.75));
      Vpdata25(p,nt)=Vpdata25(p,nt)/vtot;
      Vpdata20(p,nt) = sum(pdf.*(bins>=.20 & bins <=.80));
      Vpdata20(p,nt)=Vpdata20(p,nt)/vtot;
      Vpdata15(p,nt) = sum(pdf.*(bins>=.15 & bins <=.85));
      Vpdata15(p,nt)=Vpdata15(p,nt)/vtot;
      c1pdf(p,nt)=sum(pdf.*bins);
      c2pdf(p,nt)=sum(pdf.*(bins-c1pdf(p,nt)).^2);
      c4pdf(p,nt)=sum(pdf.*(bins-c1pdf(p,nt)).^4);
    end
  else
    Vp=zeros([10,nt]);
  end
end
end

% look at passive scalar number np
for np=1:8;

Vp25=Vpdata25(np,:);
Vp20=Vpdata20(np,:);
Vp15=Vpdata15(np,:);
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
cxx3(1,:)=pints_e(17,np,:);
cxx3(2,:)=pints_e(18,np,:);
cxx3(3,:)=pints_e(19,np,:);
cxx4(1,:)=pints_e(20,np,:);
cxx4(2,:)=pints_e(21,np,:);
cxx4(3,:)=pints_e(22,np,:);


cu(1,:)=pints_e(23,np,:);
cu(2,:)=pints_e(24,np,:);
cu(3,:)=pints_e(25,np,:);
c1=pints_e(26,np,:);
c1=squeeze(c1)';
c3=pints_e(27,np,:);
c3=squeeze(c3)';
c4=pints_e(28,np,:);
c4=squeeze(c4)';


n0(1,:)=pints_e(29,np,:);
n0_8(1,:)=pints_e(30,np,:);   % .3 after mean removed
n0_9(1,:)=pints_e(31,np,:);   % .4 after mean removed  

n0(2,:)=pints_e(32,np,:);
n0_8(2,:)=pints_e(33,np,:);
n0_9(2,:)=pints_e(34,np,:);

n0(3,:)=pints_e(35,np,:);
n0_8(3,:)=pints_e(36,np,:);
n0_9(3,:)=pints_e(37,np,:);

n0=mean(n0,1);
n0_8=mean(n0_8,1);
n0_9=mean(n0_9,1);

i=37;
for n=1:3
n0x(1,n,:)=pints_e(i+1,np,:);   % x direction N0 for c,n
n0x(2,n,:)=pints_e(i+2,np,:);   % y direction
n0x(3,n,:)=pints_e(i+3,np,:);   % z direciton
i=i+3;
end

%n0x_t=(n0x(1,1,:)+n0x(2,2,:)+n0x(3,3,:))/3;
n0x_l=(n0x(1,2,:)+n0x(1,3,:)+n0x(2,1,:)+n0x(2,3,:)+n0x(3,1,:)+n0x(3,2,:))/6;
%n0x_t=squeeze(n0x_t);
n0x_l=squeeze(n0x_l);

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
eta = (mu.^3./epsilon).^(.25);



disp(sprintf('mu = %f',mu(1) ))
disp(sprintf('ke = %f',ke(1) ))
disp(sprintf('R_l = %f',R_l(1) ))
disp(sprintf('lambda = %f',lambda(1) ))
disp(sprintf('epsilon = %f',epsilon(1) ))

epsilon_c=3*(mu./schmidt).*mean(cx2,1);
lambda_c=sqrt(c2./mean(cx2,1));
lambda_x_c=sqrt(mean(cx2,1)./mean(cxx2,1));

% mean=.5, so zero crossing of .9 is a=.4
lambda_c_9 = lambda_c.*exp(.16./(2*c2));
lambda_c_8 = lambda_c.*exp(.09./(2*c2));


%lambda_c=sqrt(  3*(mu./schmidt).*c2./epsilon_c  );
% eta_c = Bachelor scale
% lambda_c = taylor microscale
if(schmidt(1)<=1) 
  eta_c=eta./(schmidt.^(.75));
else
  eta_c=eta./sqrt(schmidt);
end





disp(' ')
disp(sprintf('passive scalar n=%i',np))
disp(sprintf('schmidt = %f',schmidt(1) ))
disp(sprintf('s2 = %f',c2(1) ))
disp(sprintf('epsilon_c = %f',epsilon_c(1) ))
disp(sprintf('lambda_c = %f',lambda_c(1) ))
disp(sprintf('1/2*pi*N0 = %f %f %f',1./(2*pi*n0(:,1)) ))



%
%data for Ray:
% 
% average over x,y,z directions:


Su = mean(ux3,1)./mean(ux2,1).^1.5  ;
Suc = mean(cu,1)./(sqrt(mean(ux2,1)).*mean(cx2,1));
G = (mean(u2,1) .* mean(uxx2,1))./mean(ux2).^2;
Gc = (c2 .* mean(cxx2,1) )./(mean(cx2,1).^2);



ff = Su.*sqrt(Rt)*7/3/sqrt(15);
ff = ff + G*7/15;

%r=3* ((lambda./lambda_c).^2) /5 ./ schmidt; 
r = 2*ke.*epsilon_c./(c2.*epsilon);
gg = sqrt(5/3)*Suc.*sqrt(Rt) + (5/9)*r.*Gc;


figure(1); clf;


subplot(4,2,1)
plot(time_e,Su,time_e,Suc)
title('S_u,  S_{uc}')

subplot(4,2,2)
plot(time_e,ff,time_e,gg)
title('f,  g')
%ax=axis;
%axis([ax(1),ax(2),0,2]);


subplot(4,2,3)
plot(time_e,G)
title('G')

subplot(4,2,4)
if (schmidt(1)<=1) 
  plot(time_e,Gc,time_e,lambda_c.^2 ./ (4*eta_c.^2))
  title('G_c Sc^{-1/2},    \lambda_c^2 / (4 \eta_c^2 )')
else
  plot(time_e,Gc,time_e,lambda_c.^2 ./ (4*eta_c.^2))
  title('G_c ,    \lambda_c^2 / (4 \eta_c^2)')
end


subplot(4,2,5)
plot(time_e,r)
title('r')
ax=axis;
axis([ax(1),ax(2),0,4]);

subplot(4,2,6)
plot(time_e,R_l)
title('R_\lambda')



subplot(4,2,7)
plot(time_e,lambda,time_e,lambda_c)
title('\lambda,  \lambda_c')

subplot(4,2,8)
plot(time_e,eta,time_e,eta_c)
title('\eta,   \eta_c ')


figure(1)
orient tall

% force all axis to be the same.
% (otherwise, when printing, the plots where we changed y axis above
% will have a different x-axis.  wierd
% Note: might instead try some of the print "render" options - see
% notes in 'help print' about screen and printed output not matching
for i=1:8
  subplot(4,2,i)
  ax=axis;
  axis([ax(1),ax(2),ax(3),ax(4)]);
end


pname=sprintf('pscalars%i',np);
print('-djpeg','-r90',pname); 


%
% some more scalars
%
figure(2); clf;


subplot(4,2,1)
plot(time_e,c2)
title('<c^2>')


subplot(4,2,2)
%plot(time_e,0*gg)
title(' ')

subplot(4,2,3)
plot(time_e,c4./c2.^2,time_e,mean(u4)./mean(u2).^2)
title('Kurtosis (c,u)')

subplot(4,2,4)
%plot(time_e,0*lambda_c)
title(' ')


subplot(4,2,5)
plot(time_e,pi*mean(nu).*lambda,time_e,pi*n0.*lambda_c,...
     time_e,pi*n0_8.*lambda_c_8,time_e,pi*n0_9.*lambda_c_9)
title('\pi N_u \lambda,    \pi N_c \lambda_c (.5,.8,.9)')


subplot(4,2,6)
plot(time_e,2*pi*n0x_l.*eta_c',time_e,2*eta_c./lambda_x_c,time_e,...
     pi*n0x_l.*lambda_x_c')
title('\pi N_{\nabla}_c 2 \eta_c, 2 \eta_c / \lambda_{\nabla}_c,   \pi N_{\nabla}_c \lambda_{\nabla}_c  ')


subplot(4,2,7)
plot(time_e,1-Vp25,time_e,1-Vp20,time_e,1-Vp15,time_e,4*c2)
title('m_{% pf}, 4<c^2>')


ym = (1 - eta_c./lambda_c).^3;
yp = (1 + eta_c./lambda_c).^3; 

subplot(4,2,8)
plot(time_e,ym./yp)
title('% pure fluid model 0')


orient tall
pname=sprintf('pscalars%i_b',np);
print('-djpeg','-r90',pname); 

end



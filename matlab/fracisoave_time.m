%
% Calculation of fractional structure function scaling exponents from time-averaged and angle-averaged data
% S.Kurien
%

mu=0;
ke=0;
nx=1;
delx_over_eta=1;
eta = 1/(nx*delx_over_eta);
ext='.isostrf';
ndir = 73;

name='/nh/xraid/skurien/dns_data/sc1024a/sc1024A'
nx=1024; delx_over_eta=2.95; epsilon=3.57;teddy = 1.05;

times=[1.4:0.1:1.5, 1.7:0.1:2.4];

ndir_use = 0;

teddy=1;


xx=(1:.5:(nx./2.5)) / nx;        %interpolation points
xx_plot = xx*nx*delx_over_eta;    %units of r/neta  
lex = round(length(xx_plot)/4);
lenx = length(xx_plot);

yyl_tangave = zeros([length(xx),9]);
yyt_tangave = zeros([length(xx),9]);

ndir_use=73;
%
% use only 49 directions:
if (ndir_use>0) 
ndir=ndir_use; 
end



if (ndir==3)
  w=ones([1,3])/3;
else
  equalw=0;
  if (equalw) 
    % put this in to use equally weighted:
    w=ones([1,ndir])/ndir;
    disp(sprintf('Using equall weights for spherical integration'))
  else
    % get the weights:
    wname=sprintf('../src/voronoi/isoave.weights%i',ndir);
    disp(sprintf('Reading Voronio weights from file:  %s',wname))
    w=textread(wname,'%f');
    % take every other weight
    w=2*w(1:2:length(w));
  end
end
if (abs(1-sum(w))>1e-7) 
  disp('error: weights do not sum to 1')
  return;
end

  close all;

times_plot=[];
for t=times
sprintf('%f',t);
tstr=sprintf('%10.4f',t+10000);
fname=[name,tstr(2:10)];
disp([fname,ext]);
     [nx,ndelta,ndir,r_val,ke,eps_l,mu,D_ll,D_lll,tmp,tmp,tmp,tmp,...
	tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,Dl,Dt,tmp] ...
= readisostr( [fname,ext] );
 

eta_l = (mu^3 / eps_l)^.25;
delx_over_eta_l=(1/nx)/eta_l;

r_val(:,:) = r_val(:,:)/nx;

ndir = ndir_use;

for ord = 1:9

yyl_angave = zeros([length(xx),9]);
yyt_angave = zeros([length(xx),9]);

for dir=1:ndir;

x=r_val(:,dir);                % units of box length
x_plot=x*nx*delx_over_eta;  % units of r/eta

yl=Dl(:,dir,ord);
yt=Dt(:,dir,ord);

yyl_angave(:,ord) = yyl_angave(:,ord) + [w(dir)*spline(x,yl,xx)]';
yyt_angave(:,ord) = yyt_angave(:,ord) + [w(dir)*spline(x,yt,xx)]';

end

yyl_tangave(:,ord) = yyl_tangave(:,ord) + yyl_angave(:,ord) ;
yyt_tangave(:,ord) = yyt_tangave(:,ord) + yyt_angave(:,ord) ;
end

end

yyl_tangave(:,:) = yyl_tangave(:,:)/length(times);
yyt_tangave(:,:) = yyt_tangave(:,:)/length(times);

for i= 1:9;
figure(i);
loglog(xx_plot(3:lenx),yyl_tangave(3:lenx,i),'r.-','LineWidth',0.1); hold on;
loglog(xx_plot(3:lenx),yyt_tangave(3:lenx,i),'k.-','LineWidth',0.1); hold on;
set(gca,'fontsize',14)
lo_str = ['\langle |\Delta u_L(r)|^{',sprintf('%3.1f',0.1*i),'}\rangle'];
tr_str = ['\langle |\Delta u_T(r)|^{',sprintf('%3.1f',0.1*i),'}\rangle'];

legend(lo_str, tr_str,4); 
xlabel('r/\eta','FontSize',14);

%ESS
%loglog(abs(yyl_tangave(:,2)),yyl_tangave(:,i),'r.-','LineWidth',1.0); hold on;
%loglog(abs(yyl_tangave(:,2)),yyt_tangave(:,i),'k.-','LineWidth',1.0); hold on;
%xlabel('|<\delta u_r ^3 >|','FontSize',16);
end

% LONGITUDINAL 
for i= 1:9;

% plots local slopes, inflection points etc.

figure(10)
der = localslp(xx_plot',yyl_tangave(:,i),1);hold on;
grid on;
set(gca,'fontsize',14)
title('Local slopes - longitudinal');
figure(15)
ind = inflect(xx_plot',yyl_tangave(:,i),3);hold on;
set(gca,'fontsize',14)
xx_plot(ind);
grid on;
title('Inflection points - longitudinal');
p=0.1*i;
sprintf('p=%d;',p)
     exp = mean(der(ind));
sprintf('Exps = %f; ',exp)
     experr = std(der(ind));
sprintf('Error on exps = %f',std(der(ind)))

     figure(17)
     k41 = p/3;
     plot(p,mean(der(ind)),'o',p,k41,'x');hold on;
     title('Comparison of exps with K41 - longitudinal');grid on;
%plot(k41,mean(der(ind)),'ro'); hold on;

figure(19)
     plot(p,(mean(der(ind))-k41)/k41,'ko');hold on;
     title('Relative scaling exponents - longitudinal');
     grid on;
     figure(21)
     errorbar(p,mean(der(ind)),experr,'o');hold on;
     title('Comparison of long and trans exps');
     grid on;
     

% PLOT THE STRUCTURE FUNCTIONS AND SCALING LINE

     y_const= 1.2*(yyl_tangave(lex,i))/(xx_plot(lex)).^exp;     
     
     figure(i)
     loglog(xx_plot, y_const*xx_plot.^exp,'r');hold on;
set(gca,'fontsize',14)
     text(xx_plot(1)-1,1.2*(yyl_tangave(4,i)),['r^{',sprintf('%3.2f',exp),'}']);


% TRANSVERSE

% plots local slopes, inflection points etc.


figure(14)
der = localslp(xx_plot',yyt_tangave(:,i),1);hold on;
title('Local slopes - transverse')
grid on;
figure(16)
ind = inflect(xx_plot',yyt_tangave(:,i),3);hold on;
xx_plot(ind);
title('Inflection points - transverse');
grid on;
sprintf('p=%d;',p)
     exp = mean(der(ind));
     sprintf('Exps = %f; ',exp)
     experr = std(der(ind));
     sprintf('Error on exps = %f',experr)

     

     figure(18)
     k41 = p/3;
     plot(p,mean(der(ind)),'o',p,k41,'x');hold on;
     title('Comparison of exps with K41 - transverse'); grid on;
%plot(k41,mean(der(ind)),'ro'); hold on;

figure(20)
     plot(p,(mean(der(ind))-k41)/k41,'ko');hold on;
     title('Relative scaling exponents -  transverse');
     grid on;
     figure(21)
     errorbar(p,mean(der(ind)),experr,'x');hold on;
     title('comparison of long and trans exps');
     grid on;


%PLOT STRUCTURE FUNCTIONS AND SCALING EXPONENTS

y_const= 1.2*(yyt_tangave(lex,i))/(xx_plot(lex)).^exp;

     figure(i)
     loglog(xx_plot, y_const*xx_plot.^exp,'k');hold on;
     set(gca,'fontsize',14)
     text(xx_plot(1)-1,1.2*(yyt_tangave(4,i)),['r^{',sprintf('%3.2f',exp),'}']);
     sname = sprintf('scalings_%i_dir%d',p,ndir);
%     print('-depsc', sname);

     figure(21)
     axis([0 1 0 0.5]) 
     cname = sprintf('exponents_dir%d',ndir);
%     print('-depsc', cname);

end
